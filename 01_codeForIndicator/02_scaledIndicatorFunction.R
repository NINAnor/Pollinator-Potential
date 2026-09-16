calculate_scaled_indicator <- function(
    richness,
    insects,
    study_region,
    classificationRegions,
    referencePolygons,
    crs,
    resultFolder,
    observed_richness,
    allCalibrationRaster,
    plant_expected_richness,
    calibration_vars,
    prediction_year,
    clc_vars,
    cores = 1,
    overwrite = TRUE,
    rerun = FALSE
) {
  
  ## ========================================================
  ## 0. Basic checks
  ## ========================================================
  
  if (!dir.exists(resultFolder)) {
    dir.create(
      resultFolder,
      recursive = TRUE
    )
  }
  
  if (!inherits(richness, "SpatRaster")) {
    stop("richness must be a terra SpatRaster.")
  }
  
  if (!inherits(insects, "SpatRaster")) {
    stop("insects must be a terra SpatRaster.")
  }
  
  if (!inherits(study_region, "SpatRaster")) {
    stop("study_region must be a terra SpatRaster.")
  }
  
  if (!inherits(referencePolygons, "sf")) {
    stop("referencePolygons must be an sf object.")
  }
  
  if (!inherits(classificationRegions, "sf")) {
    stop("classificationRegions must be an sf object.")
  }
  
  
  ## ========================================================
  ## 1. Transform spatial objects
  ## ========================================================
  
  classificationRegions <- sf::st_transform(
    classificationRegions,
    crs
  )
  
  referencePolygons <- sf::st_transform(
    referencePolygons,
    crs
  )
  
  
  ## ========================================================
  ## 2a. Calculate interaction-supported richness
  ## ========================================================
  
  richness_sum <- terra::app(
    richness,
    fun = sum,
    na.rm = TRUE,
    cores = cores,
    filename = file.path(
      resultFolder,
      "insect_richness_sum.tif"
    ),
    overwrite = overwrite
  )
  
  names(richness_sum) <- "insect_richness"
  
  ## ========================================================
  ## 2a. Calculate caliberated interaction-supported richness
  ## ======================================================== 
  if(rerun){
  insect_calibration <- fit_predict_insect_calibration(
    
    observed_richness =
      observed_richness,
    
    allCalibrationRaster =
      allCalibrationRaster,
  
    
    richness = richness_sum,
    
    nature_regions =
      classificationRegions,
    
    calibration_vars =
      calibration_vars,
    
    clc_vars =
      clc_vars,
    
    output_folder =
      plotFolder,
    
    study_region =
      study_region,
    
    prediction_year = 
      prediction_year,
    
    cores = 1,
    
    overwrite = TRUE,
    
    seed = 12345,
    
    checks = TRUE
  ) 
  
  calibrated_richness <- insect_calibration$predicted_insect_richness
  } else {
    calibrated_richness <- terra::rast(file.path(
      resultFolder,
      paste0("predicted_insect_richness_", prediction_year,".tif")
    ))   
  }
  
  # calibrated_richness <- terra::project(
  #   calibrated_richness,
  #   richness_sum,
  #   method = "bilinear"
  # )
  
  names(calibrated_richness) <-
    "calibrated_insect_richness"
  
  ## ========================================================
  ## 3. Calculate richness without interactions
  ## ========================================================
  
  richness_sum_without_interaction <- terra::app(
    insects,
    fun = sum,
    na.rm = TRUE,
    cores = cores,
    filename = file.path(
      resultFolder,
      "insect_richness_without_interaction_sum.tif"
    ),
    overwrite = overwrite
  )
  
  names(
    richness_sum_without_interaction
  ) <- "insect_richness_without_interaction"
  
  
  ## ========================================================
  ## 4. Prepare reference polygons
  ## ========================================================
  
  referencePolygons_vect <- terra::vect(
    referencePolygons
  )
  
  referencePolygons_vect <- terra::project(
    referencePolygons_vect,
    terra::crs(richness_sum)
  )
  
  
  ## ========================================================
  ## 5. Extract both richness measures
  ## ========================================================
  
  richness_sum <- terra::project(
    richness_sum,
    calibrated_richness,
    method = "bilinear"
  )
  
  richness_sum_without_interaction <- terra::project(
    richness_sum_without_interaction,
    calibrated_richness,
    method = "bilinear"
  )
  
  reference_rasters <- c(
    richness_sum,
    richness_sum_without_interaction,
    calibrated_richness
  )
  
  reference_polygon_values <- terra::extract(
    reference_rasters,
    referencePolygons_vect,
    fun = mean,
    na.rm = TRUE,
    exact = TRUE
  )
  
  
  ## ========================================================
  ## 6. Classify reference polygons by nature region
  ## ========================================================
  
  referencePolygons_classified <- sf::st_join(
    referencePolygons,
    classificationRegions["nature_region"],
    join = sf::st_within,
    left = TRUE
  )
  
  
  ## ========================================================
  ## 7. Check for unclassified reference polygons
  ## ========================================================
  
  n_unclassified <- sum(
    is.na(
      referencePolygons_classified$nature_region
    )
  )
  
  if (n_unclassified > 0) {
    
    warning(
      n_unclassified,
      " reference polygons were not assigned ",
      "to a nature region."
    )
  }
  
  
  ## ========================================================
  ## 8. Identify available reference attributes
  ## ========================================================
  
  reference_cols <- intersect(
    c(
      "naturtype",
      "naturtypeKode",
      "natrtyp",
      "ntrtypK"
    ),
    names(referencePolygons_classified)
  )
  
  
  ## ========================================================
  ## 9. Build reference-level data
  ## ========================================================
  
  reference_values_raw <- referencePolygons_classified %>%
    sf::st_drop_geometry() %>%
    dplyr::select(
      dplyr::any_of(
        reference_cols
      )
    ) %>%
    dplyr::mutate(
      
      reference_id =
        dplyr::row_number(),
      
      nature_region =
        referencePolygons_classified$nature_region,
      
      insect_richness =
        reference_polygon_values$insect_richness,
      
      insect_richness_without_interaction =
        reference_polygon_values$
        insect_richness_without_interaction,
      
      calibrated_richness = reference_polygon_values$calibrated_insect_richness
    )
  
  
  ## ========================================================
  ## 10. Remove invalid reference values
  ## ========================================================
  
  reference_values_raw <- reference_values_raw %>%
    dplyr::filter(
      !is.na(nature_region),
      !is.na(insect_richness),
      !is.na(
        insect_richness_without_interaction
      ),
      !is.na(
        calibrated_richness
      )
    ) %>%
    dplyr::filter(ntrtypK %in% c("ntyp_D02_02", "ntyp_D02"))
  
  
  ## ========================================================
  ## 11. Calculate regional reference statistics
  ## ========================================================
  
  reference_regional_richness <- reference_values_raw %>%
    dplyr::group_by(
      nature_region
    ) %>%
    dplyr::summarise(
      
      ## --------------------------------------------
      ## Interaction-supported richness
      ## --------------------------------------------
      
      richness = mean(
        insect_richness,
        na.rm = TRUE
      ),
      
      min_richness = min(
        insect_richness,
        na.rm = TRUE
      ),
      
      max_richness = max(
        insect_richness,
        na.rm = TRUE
      ),
      
      ## --------------------------------------------
      ## Richness without interaction
      ## --------------------------------------------
      
      richness_without_interaction =
        mean(
          insect_richness_without_interaction,
          na.rm = TRUE
        ),
      
      min_richness_without_interaction =
        min(
          insect_richness_without_interaction,
          na.rm = TRUE
        ),
      
      max_richness_without_interaction =
        max(
          insect_richness_without_interaction,
          na.rm = TRUE
        ),
      
      ## --------------------------------------------
      ## Caliberated richness
      ## --------------------------------------------
      
      cal_richness =
        mean(
          calibrated_richness,
          na.rm = TRUE
        ),
      
      min_cal_richness =
        min(
          calibrated_richness,
          na.rm = TRUE
        ),
      
      max_cal_richness =
        max(
          calibrated_richness,
          na.rm = TRUE
        ),
      
      
      ## --------------------------------------------
      ## Number of reference polygons
      ## --------------------------------------------
      
      n_reference = dplyr::n(),
      
      .groups = "drop"
    )
  
  
  ## ========================================================
  ## 12. Calculate GLOBAL reference ranges
  ## ========================================================
  
  ## Interaction-supported
  
  # global_ref_min <- min(
  #   reference_values_raw$insect_richness,
  #   na.rm = TRUE
  # )
  
  global_ref_min <- 0
  global_ref_max <- max(
    reference_values_raw$insect_richness,
    na.rm = TRUE
  )
  
  
  ## Without interaction
  
  # global_ref_min_without_interaction <- min(
  #   reference_values_raw$
  #     insect_richness_without_interaction,
  #   na.rm = TRUE
  # )
  
  global_ref_min_without_interaction <- 0
  
  global_ref_max_without_interaction <- max(
    reference_values_raw$
      insect_richness_without_interaction,
    na.rm = TRUE
  )
  
  # With caliberation
  
  # global_ref_min_caliberation <- min(
  #   reference_values_raw$calibrated_richness,
  #   na.rm = TRUE
  # )
  
  global_ref_min_caliberation <- 0
  
  global_ref_max_caliberation <- max(
    reference_values_raw$calibrated_richness,
    na.rm = TRUE
  )
  
  ## ========================================================
  ## 13. Validate global reference ranges
  ## ========================================================
  
  if (
    global_ref_max <=
    global_ref_min
  ) {
    
    stop(
      "Global interaction-supported ",
      "reference maximum must be greater ",
      "than the reference minimum."
    )
  }
  
  
  if (
    global_ref_max_without_interaction <=
    global_ref_min_without_interaction
  ) {
    
    stop(
      "Global non-interaction ",
      "reference maximum must be greater ",
      "than the reference minimum."
    )
  }
  
  
  ## ========================================================
  ## 14. Convert classification regions to terra
  ## ========================================================
  
  classificationRegions_vect <- terra::vect(
    classificationRegions
  )
  
  classificationRegions_vect <- terra::project(
    classificationRegions_vect,
    terra::crs(richness_sum)
  )
  
  
  ## ========================================================
  ## 15. Assign numeric region IDs
  ## ========================================================
  
  classificationRegions_vect$region_id <-
    seq_len(
      nrow(
        classificationRegions_vect
      )
    )
  
  
  ## ========================================================
  ## 16. Rasterize nature regions
  ## ========================================================
  
  region_raster <- terra::rasterize(
    classificationRegions_vect,
    richness_sum,
    field = "region_id"
  )
  
  names(region_raster) <- "region_id"
  
  
  ## ========================================================
  ## 17. Create regional reference lookup
  ## ========================================================
  
  region_lookup <- data.frame(
    
    region_id =
      classificationRegions_vect$region_id,
    
    nature_region =
      classificationRegions_vect$nature_region
    
  ) %>%
    
    dplyr::left_join(
      reference_regional_richness %>%
        dplyr::select(
          
          nature_region,
          
          ## Interaction-supported
          min_richness,
          max_richness,
          
          ## Without interaction
          min_richness_without_interaction,
          max_richness_without_interaction,
          
          #Caliberation
          min_cal_richness,
          max_cal_richness 
        ),
      by = "nature_region"
    )
  
  
  ## ========================================================
  ## 18. Check regional reference values
  ## ========================================================
  
  missing_regions <- region_lookup %>%
    
    dplyr::filter(
      
      is.na(min_richness) |
        is.na(max_richness) |
        
        is.na(
          min_richness_without_interaction
        ) |
        
        is.na(
          max_richness_without_interaction
        )
      
    ) %>%
    
    dplyr::pull(
      nature_region
    )
  
  
  if (
    length(missing_regions) > 0
  ) {
    
    stop(
      "No reference min/max values found for: ",
      paste(
        unique(missing_regions),
        collapse = ", "
      )
    )
  }
  
  
  ## ========================================================
  ## 19. Check regional reference ranges
  ## ========================================================
  
  invalid_regions <- region_lookup %>%
    
    dplyr::filter(
      
      max_richness <=
        min_richness |
        
        max_richness_without_interaction <=
        min_richness_without_interaction |
        
        max_cal_richness <= min_cal_richness 
      
    ) %>%
    
    dplyr::pull(
      nature_region
    )
  
  
  if (
    length(invalid_regions) > 0
  ) {
    
    stop(
      "Reference maximum must be greater ",
      "than reference minimum for: ",
      paste(
        unique(invalid_regions),
        collapse = ", "
      )
    )
  }
  
  
  ## ========================================================
  ## 20. Combine richness layers and region ID
  ## ========================================================
  
  richness_sum_region <- c(
    
    richness_sum,
    
    richness_sum_without_interaction,
    
    calibrated_richness,
    
    region_raster
    
  )
  
  names(richness_sum_region) <- c(
    
    "insect_richness",
    
    "insect_richness_without_interaction",
    
    "calibrated_richness",
    
    "region_id"
    
  )
  
  
  ## ========================================================
  ## 21. Reference lookup vectors
  ## ========================================================
  
  ## Interaction-supported
  
  #ref_min <- region_lookup$min_richness
  ref_min <- rep(0, length(region_lookup$min_richness))
  
  ref_max <- region_lookup$max_richness
  
  
  ## Without interaction
  
  ref_min_without_interaction <-
  rep(0,  length(region_lookup$min_richness_without_interaction))
  
  ref_max_without_interaction <-
    region_lookup$
    max_richness_without_interaction
  
  ref_cal_min <- rep(0, length(region_lookup$min_cal_richness))
  ref_cal_max <- region_lookup$max_cal_richness
  
  
  ## ========================================================
  ## 22. Regional scaling
  ## ========================================================
  
  regional_scaled <- terra::lapp(
    
    richness_sum_region,
    
    fun = function(
    richness,
    richness_without_interaction,
    calibrated_richness,
    region_id
    ) {
      
      ## --------------------------------------------
      ## Interaction-supported indicator
      ## --------------------------------------------
      
      min_value <-
        ref_min[region_id]
      
      max_value <-
        ref_max[region_id]
      
      regional_indicator <-
        (
          richness - min_value
        ) /
        (
          max_value - min_value
        )
      
      ## --------------------------------------------
      ## Calibrated
      ## --------------------------------------------
      
      min_value <-
        ref_cal_min[region_id]
      
      max_value <-
        ref_cal_max[region_id]
      
      cal_regional_indicator <-
        (
          calibrated_richness - min_value
        ) /
        (
          max_value - min_value
        )
      
      
      ## --------------------------------------------
      ## Without-interaction indicator
      ## --------------------------------------------
      
      min_value_no_interaction <-
        ref_min_without_interaction[
          region_id
        ]
      
      max_value_no_interaction <-
        ref_max_without_interaction[
          region_id
        ]
      
      regional_indicator_no_interaction <-
        (
          richness_without_interaction -
            min_value_no_interaction
        ) /
        (
          max_value_no_interaction -
            min_value_no_interaction
        )
      
      
      ## --------------------------------------------
      ## Return both layers
      ## --------------------------------------------
      
      c(
        regional_indicator,
        cal_regional_indicator,
        regional_indicator_no_interaction
      )
    },
    
    cores = cores,
    
    filename = file.path(
      resultFolder,
      "scaled_insect_indicators_regional.tif"
    ),
    
    overwrite = overwrite
  )
  
  
  names(regional_scaled) <- c(
    
    "regional_scaled",
    "caliberated_regional_scaled",
    
    "regional_scaled_without_interaction"
    
  )
  
  
  ## ========================================================
  ## 23. Mask regional indicators
  ## ========================================================
  
  regional_masked <- terra::mask(
    regional_scaled,
    study_region,
    #maskvalues = 0,
    updatevalue = NA
  )
  
  
  names(regional_masked) <- c(
    
    "regional_masked",
    
    "caliberated_regional_masked",
    
    "regional_masked_without_interaction"
    
  )
  
  
  ## ========================================================
  ## 24. Clamp regional indicators
  ## ========================================================
  
  regional_clamped <- terra::clamp(
    regional_masked,
    lower = 0,
    upper = 1,
    values = TRUE
  )
  
  
  names(regional_clamped) <- c(
    
    "regional_clamped",
    "cal_regional_clamped",
    
    "regional_clamped_without_interaction"
    
  )
  
  
  ## ========================================================
  ## 25. Global scaling
  ## ========================================================
  
  global_scaled <- terra::lapp(
    
    c(
      richness_sum,
      calibrated_richness,
      richness_sum_without_interaction
    ),
    
    fun = function(
    richness,
    calibrated_richness,
    richness_without_interaction
    ) {
      
      ## --------------------------------------------
      ## Interaction-supported
      ## --------------------------------------------
      
      global_indicator <-
        (
          richness -
            global_ref_min
        ) /
        (
          global_ref_max -
            global_ref_min
        )
      
      ## --------------------------------------------
      ## Caliberated 
      ## --------------------------------------------
      
      global_cal_indicator <-
        (
          calibrated_richness -
            global_ref_min_caliberation
        ) /
        (
          global_ref_max_caliberation -
            global_ref_min_caliberation
        )
      
      
      ## --------------------------------------------
      ## Without interaction
      ## --------------------------------------------
      
      global_indicator_no_interaction <-
        (
          richness_without_interaction -
            global_ref_min_without_interaction
        ) /
        (
          global_ref_max_without_interaction -
            global_ref_min_without_interaction
        )
      
      
      ## --------------------------------------------
      ## Return both
      ## --------------------------------------------
      
      c(
        global_indicator,
        global_cal_indicator,
        global_indicator_no_interaction
      )
    },
    
    cores = cores,
    
    filename = file.path(
      resultFolder,
      "scaled_insect_indicators_global.tif"
    ),
    
    overwrite = overwrite
  )
  
  
  names(global_scaled) <- c(
    
    "global_scaled",
    "global_cal_scaled",
    
    "global_scaled_without_interaction"
    
  )
  
  
  ## ========================================================
  ## 26. Mask global indicators
  ## ========================================================
  
  global_masked <- terra::mask(
    global_scaled,
    study_region#,
    #maskvalues = 0
  )
  
  
  names(global_masked) <- c(
    
    "global_masked",
    "global_cal_masked",
    "global_masked_without_interaction"
    
  )
  
  
  ## ========================================================
  ## 27. Clamp global indicators
  ## ========================================================
  
  global_clamped <- terra::clamp(
    global_masked,
    lower = 0,
    upper = 1,
    values = TRUE
  )
  
  
  names(global_clamped) <- c(
    
    "global_clamped",
    "global_cal_clamped",
    "global_clamped_without_interaction"
    
  )
  
  
  ## ========================================================
  ## 28. Combine all raster outputs
  ## ========================================================
  
  allRasts <- c(
    
    ## --------------------------------------------
    ## Raw richness
    ## --------------------------------------------
    
    richness_sum,
    calibrated_richness,
    richness_sum_without_interaction,
    
    
    ## --------------------------------------------
    ## Region ID
    ## --------------------------------------------
    
    region_raster,
    
    
    ## --------------------------------------------
    ## Regional indicators
    ## --------------------------------------------
    
    regional_clamped,
    
    regional_masked,
    
    regional_scaled,
    
    
    ## --------------------------------------------
    ## Global indicators
    ## --------------------------------------------
    
    global_clamped,
    
    global_masked,
    
    global_scaled
    
  )
  
  
  names(allRasts) <- c(
    
    ## Raw richness
    
    "richness",
    "calibrated_richness",
    
    "richness_without_interaction",
    
    
    ## Region
    
    "region_id",
    
    
    ## Regional
    
    "regional_clamped",
    "calibrated_regional_clamped",
    
    "regional_clamped_without_interaction",
    
    "regional_masked",
    "calibrated_regional_masked",
    
    "regional_masked_without_interaction",
    
    "regional_scaled",
    "calibrated_regional_scaled",
    
    "regional_scaled_without_interaction",
    
    
    ## Global
    
    "global_clamped",
    "calibrated_global_clamped",
    
    "global_clamped_without_interaction",
    
    "global_masked",
    "calibrated_global_masked",
    "global_masked_without_interaction",
    
    "global_scaled",
    "calibrated_global_scaled",
    
    "global_scaled_without_interaction"
    
  )
  
  
  terra::writeRaster(allRasts,
              filename = file.path(
                resultFolder,
                "all_indicator_raster.tif"
              ),
              overwrite = overwrite)
  
  ## ========================================================
  ## 29. Return results
  ## ========================================================
  
  list(
    
    ## ------------------------------------------------------
    ## All rasters
    ## ------------------------------------------------------
    
    allRasts =
      allRasts,
    
    
    ## ------------------------------------------------------
    ## Raw richness
    ## ------------------------------------------------------
    
    richness_sum =
      richness_sum,
    
    richness_sum_without_interaction =
      richness_sum_without_interaction,
    
    calibrated_richness = 
      calibrated_richness,
    
    
    ## ------------------------------------------------------
    ## Region raster
    ## ------------------------------------------------------
    
    region_raster =
      region_raster,
    
    richness_sum_region =
      richness_sum_region,
    
    
    ## ------------------------------------------------------
    ## Regional indicators
    ## ------------------------------------------------------
    
    scaled_richness =
      regional_scaled,
    
    scaled_richness_masked =
      regional_masked,
    
    scaled_richness_clamped =
      regional_clamped,
    
    
    ## ------------------------------------------------------
    ## Global indicators
    ## ------------------------------------------------------
    
    scaled_richness_global =
      global_scaled,
    
    scaled_richness_global_masked =
      global_masked,
    
    scaled_richness_global_clamped =
      global_clamped,
    
    
    ## ------------------------------------------------------
    ## Reference information
    ## ------------------------------------------------------
    
    reference_values =
      reference_regional_richness,
    
    reference_values_raw =
      reference_values_raw,
    
    reference_polygons =
      referencePolygons_classified,
    
    region_lookup =
      region_lookup,
    
    
    ## ------------------------------------------------------
    ## Global reference levels
    ## ------------------------------------------------------
    
    global_ref_min =
      global_ref_min,
    
    global_ref_max =
      global_ref_max,
    
    global_ref_min_without_interaction =
      global_ref_min_without_interaction,
    
    global_ref_max_without_interaction =
      global_ref_max_without_interaction
    
  )
}