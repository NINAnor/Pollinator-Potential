################################################################################
# ESTIMATE INDICATOR VALUES
#
# For each indicator and Nature region, calculate:
#
#   1. Global mean
#   2. Global SD
#   3. Proportion of area >= 0.6
#   4. Bootstrap mean
#   5. Bootstrap SD
#   6. 2.5% bootstrap quantile
#   7. 97.5% bootstrap quantile
#
# An additional national estimate is calculated across the entire study area.
#
# The representative value is estimated using an intercept-only beta regression:
#
#       indicator ~ 1
#
# The model intercept is estimated on the logit scale and transformed back
# to the original 0-1 scale using the inverse-logit transformation.
################################################################################


################################################################################
# 1. Packages
################################################################################

library(terra)
library(sf)
library(dplyr)
library(purrr)
library(tibble)
library(betareg)

library(sf)

## ========================================================
## Path to mountain-area geodatabase
## ========================================================

library(terra)




################################################################################
# 2. Load indicator raster
################################################################################

if (!file.exists(
  file.path(
    plotFolder,
    "region_map.shp"
  )
)) {
  
  allRasts <- terra::rast(
    file.path(
      plotFolder,
      "all_indicator_raster.tif"
    )
  )
  
  ## ========================================================
  ## Mask indicator rasters to exclude mountain areas
  ## ========================================================
  
  masked_file <- file.path(
    plotFolder,
    "all_indicator_raster_no_mountains.tif"
  )
  
  
  ## ========================================================
  ## Mask mountain areas from all indicator rasters
  ## ========================================================
  
  masked_file <- file.path(
    plotFolder,
    "all_indicator_raster_no_mountains.tif"
  )
  
  
  ## --------------------------------------------------------
  ## Read existing result if available
  ## --------------------------------------------------------
  
  if (file.exists(masked_file)) {
    
    allRasts_masked <- terra::rast(
      masked_file
    )
    
  } else {
    
    input_file <- file.path(
      "R:/GeoSpatialData",
      "LandCover",
      "Norway_MountainAreas",
      "Original",
      "Fjell_norge_publisering.gdb"
    )
    
    
    ## ========================================================
    ## Read 50 m mountain-area grid
    ## ========================================================
    
    fjell_grid50 <- sf::st_read(
      dsn = input_file,
      layer = "fjell_grid50",
      quiet = FALSE
    )
    
    library(terra)
    
    
    ## Convert mountain polygons to SpatVector
    fjell_vect <- terra::vect(
      fjell_grid50
    )
    
    
    ## ------------------------------------------------------
    ## Create a raster mask using the grid of allRasts
    ##
    ## 1 = mountain
    ## 0/NA = outside mountain
    ## ------------------------------------------------------
    
    mountain_mask <- terra::rasterize(
      fjell_vect,
      allRasts[[1]],
      field = 1,
      background = 0
    )
    
    
    ## ------------------------------------------------------
    ## Set mountain cells to NA
    ##
    ## !mountain_mask gives:
    ##   TRUE  = outside mountains
    ##   FALSE = mountains
    ## ------------------------------------------------------
    
    allRasts_masked <- terra::mask(
      allRasts,
      mountain_mask,
      maskvalues = 1,
      updatevalue = NA
    )
    
    
    ## ------------------------------------------------------
    ## Save
    ## ------------------------------------------------------
    
    terra::writeRaster(
      allRasts_masked,
      masked_file,
      overwrite = TRUE
    )
  }
  
  ##############################################################################
  # Extract the six indicator layers
  #
  # region_id is kept separately because it is used to identify Nature regions,
  # but it is NOT used in the beta regression.
  ##############################################################################
  
  indicator_names <- c(
    "Regional indicator",
    "Calibrated regional indicator",
    "Regional indicator without interaction",
    "Global indicator",
    "Calibrated global indicator",
    "Global indicator without interaction"
  )
  
  indicator_layer_names <- c(
    "regional_clamped",
    "calibrated_regional_clamped",
    "regional_clamped_without_interaction",
    "global_clamped",
    "calibrated_global_clamped",
    "global_clamped_without_interaction"
  )
  
  
  ##############################################################################
  # Create indicator raster stack
  ##############################################################################
  
  indicator_rasters <- allRasts_masked[[indicator_layer_names]]
  
  names(indicator_rasters) <- indicator_names
  
  
  ##############################################################################
  # Extract Nature-region ID raster
  ##############################################################################
  
  region_id <- allRasts_masked[["region_id"]]
  
  
  ##############################################################################
  # Region lookup table
  ##############################################################################
  insect_indicator_region_lookup <- readRDS(paste0(plotFolder, "/insect_indicator_region_lookup.rds"))
  region_lookup <- insect_indicator_region_lookup %>%
    dplyr::select(
      region_id,
      nature_region
    )
  
  
  ##############################################################################
  # Bootstrap settings
  ##############################################################################
  
  B <- 50
  
  sample_size <- 2000
  
  bootstrap_seed <- 1000
  
  
  ##############################################################################
  # Pre-build beta-regression formulas
  #
  # All models are intercept-only:
  #
  #       indicator ~ 1
  #
  # Backticks are required because the indicator names contain spaces.
  ##############################################################################
  
  indicator_formulas <- setNames(
    lapply(
      indicator_names,
      function(indicator) {
        
        stats::as.formula(
          paste0(
            "`",
            indicator,
            "` ~ 1"
          )
        )
      }
    ),
    indicator_names
  )
  
  
  ################################################################################
  # FUNCTION 1: Calculate bootstrap estimates
  ################################################################################
  #
  # This function takes a raster containing only the six indicators and repeatedly
  # samples cells and fits an intercept-only beta regression.
  #
  # The intercept is transformed from the logit scale to the original 0-1 scale.
  ################################################################################
  
  bootstrap_indicator_values <- function(
    indicator_rasters,
    B = 100,
    sample_size = 10000,
    seed = 1000
  ) {
    
    purrr::map_dfr(
      seq_len(B),
      function(i) {
        
        message(
          "  Bootstrap ",
          i,
          "/",
          B
        )
        
        
        ##########################################################################
        # Give every bootstrap replicate a unique seed
        ##########################################################################
        
        set.seed(
          seed + i
        )
        
        
        ##########################################################################
        # Sample raster cells
        #
        # Only the six indicator layers are sampled.
        #
        # xy coordinates are not required for the beta regression, so they are
        # omitted to reduce the amount of data returned.
        ##########################################################################
        
        indicator_sample <- terra::spatSample(
          indicator_rasters,
          size = sample_size,
          method = "random",
          na.rm = TRUE,
          as.df = TRUE,
          xy = TRUE
        ) %>%
          na.omit()
        
        
        ##########################################################################
        # Fit one beta regression for each indicator
        ##########################################################################
        
        estimates <- purrr::map_dbl(
          indicator_names,
          function(indicator) {
            
            model <- betareg::betareg(
              indicator_formulas[[indicator]],
              data = indicator_sample,
              link = "logit"
            )
            
            
            ######################################################################
            # Extract intercept and transform from logit to 0-1 scale
            ######################################################################
            
            intercept_logit <- coef(
              model,
              model = "mean"
            )[["(Intercept)"]]
            
            stats::plogis(
              intercept_logit
            )
          }
        )
        
        
        ##########################################################################
        # Return bootstrap estimates
        ##########################################################################
        
        tibble::tibble(
          replicate = i,
          indicator = indicator_names,
          estimate = estimates
        )
      }
    )
  }
  
  
  ################################################################################
  # FUNCTION 2: Summarise bootstrap distribution
  ################################################################################
  
  summarise_bootstrap <- function(
    boot_results
  ) {
    
    boot_results %>%
      dplyr::group_by(
        indicator
      ) %>%
      dplyr::summarise(
        
        bootstrap_mean = mean(
          estimate,
          na.rm = TRUE
        ),
        
        bootstrap_sd = sd(
          estimate,
          na.rm = TRUE
        ),
        
        lower_quantile = quantile(
          estimate,
          probs = 0.025,
          na.rm = TRUE
        ),
        
        upper_quantile = quantile(
          estimate,
          probs = 0.975,
          na.rm = TRUE
        ),
        
        .groups = "drop"
      )
  }
  
  
  ################################################################################
  # FUNCTION 3: Calculate statistics for one spatial area
  ################################################################################
  #
  # This function:
  #
  #   1. Crops the indicator raster to a polygon
  #   2. Masks it to the exact polygon boundary
  #   3. Calculates mean, SD and proportion >= 0.6
  #   4. Runs the beta-regression bootstrap
  #   5. Combines everything into one result
  #
  ################################################################################
  
  calculate_indicator_statistics <- function(
    indicator_rasters,
    polygon = NULL,
    label = "national",
    B = 100,
    sample_size = 10000,
    seed = 1000
  ) {
    
    message(
      "\nProcessing: ",
      label
    )
    
    
    ##############################################################################
    # Crop and mask to polygon when a polygon is supplied
    ##############################################################################
    
    if (!is.null(polygon)) {
      
      indicator_vals <- terra::crop(
        indicator_rasters,
        polygon
      ) %>%
        terra::mask(
          polygon
        )
      
    } else {
      
      indicator_vals <- indicator_rasters
      
    }
    
    
    ##############################################################################
    # Calculate raster mean and SD
    ##############################################################################
    
    raster_stats <- terra::global(
      indicator_vals,
      fun = c(
        "mean",
        "sd"
      ),
      na.rm = TRUE
    )
    
    
    ##############################################################################
    # Calculate proportion of valid cells >= 0.6
    ##############################################################################
    
    proportion_ge_06 <- terra::global(
      indicator_vals >= 0.6,
      fun = "mean",
      na.rm = TRUE
    )[, 1]
    
    
    ##############################################################################
    # Bootstrap
    ##############################################################################
    
    boot_results <- bootstrap_indicator_values(
      indicator_rasters = indicator_vals,
      B = B,
      sample_size = sample_size,
      seed = seed
    )
    
    
    ##############################################################################
    # Summarise bootstrap
    ##############################################################################
    
    bootstrap_summary <- summarise_bootstrap(
      boot_results
    )
    
    
    ##############################################################################
    # Combine raster and bootstrap statistics
    ##############################################################################
    
    results <- tibble::tibble(
      
      nature_region = label,
      
      indicator = indicator_names,
      
      global_mean = raster_stats$mean,
      
      global_sd = raster_stats$sd,
      
      proportion_ge_0.6 = proportion_ge_06
      
    ) %>%
      
      dplyr::left_join(
        bootstrap_summary,
        by = "indicator"
      ) %>%
      
      dplyr::mutate(
        percentage_ge_0.6 =
          proportion_ge_0.6 * 100
      )
    
    
    return(results)
  }
  
  
  ################################################################################
  # 8. Calculate national estimates
  ################################################################################
  national_rasters <- terra::aggregate(
                                        indicator_rasters,
                                          fact = 20,
                                          fun = mean,
                                         na.rm = TRUE)
  
  national_results <- calculate_indicator_statistics(
    
    indicator_rasters =  national_rasters,
    
    polygon = NULL,
    
    label = "national",
    
    B = B,
    
    sample_size = sample_size,
    
    seed = bootstrap_seed
  )
  
  
  ################################################################################
  # 9. Calculate Nature-region estimates
  ################################################################################
  
  regions <- unique(
    nature_regions$nature_region
  )
  
  
  region_results <- purrr::map_dfr(
    
    seq_along(regions),
    
    function(i) {
      
      region <- regions[i]
      
      
      ##########################################################################
      # Extract Nature-region polygon
      ##########################################################################
      
      region_polygon <- nature_regions %>%
        dplyr::filter(
          nature_region == region
        ) %>%
        terra::vect()
      
      
      ##########################################################################
      # Use a different seed for every region
      #
      # This prevents the same bootstrap random-number sequence being reused
      # across regions.
      ##########################################################################
      
      calculate_indicator_statistics(
        
        indicator_rasters = indicator_rasters,
        
        polygon = region_polygon,
        
        label = region,
        
        B = B,
        
        sample_size = sample_size,
        
        seed = bootstrap_seed + i * 10000
      )
    }
  )
  
  
  ################################################################################
  # 10. Combine national and regional results
  ################################################################################
  
  allIndicatorsEstimated <- dplyr::bind_rows(
    
    national_results,
    
    region_results
    
  )
  
  
  ################################################################################
  # 11. Save indicator estimates
  ################################################################################
  
  saveRDS(
    allIndicatorsEstimated,
    file = file.path(
      plotFolder,
      "region_results.rds"
    )
  )
  
  
  ################################################################################
  # 12. Create Nature-region map
  ################################################################################
  
  region_polygons <- nature_regions %>%
    sf::st_as_sf() %>%
    dplyr::rename(
      geometry = SHAPE
    ) %>%
    dplyr::select(
      nature_region,
      geometry
    )
  
  
  ################################################################################
  # Add indicator results to Nature-region polygons
  ################################################################################
  
  region_map <- region_polygons %>%
    dplyr::left_join(
      region_results,
      by = "nature_region"
    )
  
  
  ################################################################################
  # Save Nature-region map
  ################################################################################
  
  sf::st_write(
    region_map,
    file.path(
      plotFolder,
      "region_map.shp"
    ),
    delete_dsn = TRUE,
    quiet = TRUE
  )
  
}


################################################################################
# 13. Load existing region map
################################################################################

region_map <- sf::st_read(
  file.path(
    plotFolder,
    "region_map.shp"
  ),
  quiet = TRUE
)


################################################################################
# END
################################################################################

# p <- ggplot() +
#   
#   # Actual nature-region polygons
#   geom_sf(
#     data = region_map,
#     aes(fill = ntr_rgn),
#     colour = "black",
#     linewidth = 0.4
#   ) +
#   
#   # Values inside each region
#   geom_sf_text(
#     data = region_labels,
#     aes(label = label),
#     size = 3.5
#   ) +
#   
#   # Six indicator maps
#   facet_wrap(
#     ~ indictr,
#     ncol = 3
#   ) +
#   
#   scale_fill_brewer(
#     palette = "Set2",
#     name = "Nature region"
#   ) +
#   
#   labs(
#     x = NULL,
#     y = NULL
#   ) +
#   
#   theme_minimal(base_size = 12) +
#   theme(
#     axis.text = element_blank(),
#     axis.ticks = element_blank(),
#     panel.grid = element_blank(),
#     strip.text = element_text(
#       face = "bold",
#       size = 11
#     ),
#     legend.position = "bottom",
#     panel.spacing = unit(1, "lines")
#   )
# 
# ggsave(filename = file.path(plotFolder,
#                             "indicator_map.png"),
#        plot = p,
#        units = "in",
#        width = 9,
#        height = 10)
# 
# ggsave(filename = file.path(plotFolder,
#                             "indicator_map.pdf"),
#        plot = p,
#        units = "in",
#        width = 9,
#        height = 10)
