################################################################################
# Create open-lowland mask
#
# This script creates a binary raster identifying areas classified as
# "open lowland".
#
# Two alternative methods can be used:
#
#   Method 1: CLC+ land-cover classification
#       useCLCForLowlandBoundary = TRUE
#
#       Open lowland is defined using CLC+ 10-m land-cover classes 6 and 7.
#
#   Method 2: Percentage-based threshold
#       useCLCForLowlandBoundary = FALSE
#
#       Open lowland is defined as areas where the percentage of open lowland
#       exceeds the value specified by `openlowland_threshold`.
#
# In both cases, the final output is a binary SpatRaster:
#
#       1  = open-lowland area
#       NA = outside the open-lowland area
#
################################################################################

# ==============================================================================
# METHOD 1: Define open lowland using CLC+ land-cover classes
# ==============================================================================

if (useCLCplusForLowlandBoundary) {
  
  
  # ============================================================================
  # 1. Define the output file
  # ============================================================================
  #
  # The mask is saved at 50-m resolution so that it can be used directly with
  # the other 50-m calibration/prediction rasters.
  # ============================================================================
  
  openlowland_mask_file <- file.path(
    resultFolder,
    "polygons",
    "openlowland_mask_CLC_50m.tif"
  )
  
  
  # ============================================================================
  # 2. Create the 50-m open-lowland mask if it does not already exist
  # ============================================================================
  
  if (!file.exists(openlowland_mask_file)) {
    
    # --------------------------------------------------------------------------
    # 2.1 Load the CLC+ 10-m land-cover raster
    # --------------------------------------------------------------------------
    #
    # CLC+ is available at 10-m resolution.
    #
    # `project_dir` points to the location of the CLC+ data.
    # --------------------------------------------------------------------------
    
    project_dir <- "C:/Users/KWAKU~1.ADJ/ONEDRI~1/222608~1"
    
    CLCplus10m <- terra::rast(
      file.path(
        project_dir,
        "CLC+/104035/Results/CLCplus_2018_010m/",
        "CLCplus_2018_010m/CLCplus_2018_010m.tif"
      )
    )
    
    
    # --------------------------------------------------------------------------
    # 2.2 Convert CLC+ classes into a binary open-lowland mask
    # --------------------------------------------------------------------------
    #
    # CLC+ classes 6 and 7 are being used here to represent open lowland.
    #
    # The resulting raster contains:
    #
    #       1  = CLC+ class 6 or 7 (open lowland)
    #       NA = all other land-cover classes
    #
    # `ifel()` is used instead of assigning values directly so that all
    # non-open-lowland cells become NA.
    # --------------------------------------------------------------------------
    
    openlowland_binary <- terra::ifel(
      CLCplus10m %in% c(6, 7),
      1,
      NA
    )
    
    
    # --------------------------------------------------------------------------
    # 2.3 Project the 10-m binary mask onto the 50-m calibration grid
    # --------------------------------------------------------------------------
    #
    # `allCalibrationRaster[[1]]` is used as the target grid.
    #
    # This ensures that the resulting mask has the same:
    #
    #       - CRS
    #       - extent
    #       - resolution
    #       - cell alignment
    #
    # as the calibration raster.
    #
    # Nearest-neighbour ("near") resampling is used because this is a
    # categorical/binary raster. We do NOT want interpolation to create
    # fractional values such as 0.4 or 0.7.
    #
    # The result is saved directly to disk.
    # --------------------------------------------------------------------------
    
    openlowland_50m <- terra::project(
      openlowland_binary,
      allCalibrationRaster[[1]],
      method = "near",
      filename = openlowland_mask_file,
      overwrite = TRUE
    )
  }
  
  
  # ============================================================================
  # 3. Load the CLC+-based 50-m open-lowland mask
  # ============================================================================
  
  openlowland_binary <- terra::rast(
    openlowland_mask_file
  )
}


if(useCLCForLowlandBoundary){
  project_dir <- "C:/Users/KWAKU~1.ADJ/ONEDRI~1/222608~1"
  
  CLC100m <- terra::rast(
    file.path(
      project_dir,
      "CLC100/u2018_clc2018_v2020_20u1_raster100m/u2018_clc2018_v2020_20u1_raster100m/DATA/U2018_CLC2018_V2020_20u1.tif"
    )
  )
  
  ## ========================================================
  ## Extract CLC classes representing open lowland
  ## ========================================================
  
  open_lowland_classes <- c(
    18,  # Pastures
    26,  # Natural grasslands
    27   # Moors and heathland
  )
  
  open_lowland <- terra::ifel(
    CLC100m %in% open_lowland_classes,
    1,
    NA
  )
  
  open_lowland_crop <- terra::crop(open_lowland, 
                                   terra::project(allCalibrationRaster[[1]], crs(open_lowland))) %>%
    terra::mask(terra::project(allCalibrationRaster[[1]], crs(open_lowland)))
  
  
  
  openlowland_50m <- terra::project(
    open_lowland_crop,
    allCalibrationRaster[[1]],
    method = "near",
    filename = openlowland_mask_file,
    overwrite = TRUE
  )
  
 
  

  
  }

# ==============================================================================
# METHOD 2: Define open lowland using a percentage threshold
# ==============================================================================

if (!useCLCForLowlandBoundary && !useCLCplusForLowlandBoundary) {
  
  source("01_codeForEnvAgency/04_process_ASO_data.R")
  # ============================================================================
  # 4. Define the output file
  # ============================================================================
  #
  # The threshold is included in the filename so that different thresholds
  # produce different mask files.
  #
  # For example:
  #
  #     openlowland_mask_50pct.tif
  #
  # would represent a 50% open-lowland threshold.
  # ============================================================================
  
  openlowland_mask_file <- file.path(
    "C:/terra_tmp",
    paste0(
      "openlowland_mask_",
      openlowland_threshold,
      "pct.tif"
    )
  )
  
  
  # ============================================================================
  # 5. Create the threshold-based mask if it does not already exist
  # ============================================================================
  
  if (!file.exists(openlowland_mask_file)) {
    
    
    # --------------------------------------------------------------------------
    # 5.1 Load the open-lowland percentage raster
    # --------------------------------------------------------------------------
    #
    # This raster contains the percentage of open-lowland habitat for each
    # raster cell.
    # --------------------------------------------------------------------------
    
    openlowland_boundary <- terra::rast(
      file.path(
        resultFolder,
        "openlow_pct_NO.tif"
      )
    )
    
    
    # --------------------------------------------------------------------------
    # 5.2 Project the percentage raster onto the prediction grid
    # --------------------------------------------------------------------------
    #
    # `plants` is used as the target raster, ensuring that the open-lowland
    # percentage data are aligned with the plant/insect prediction grid.
    #
    # Bilinear interpolation is appropriate here because the input is a
    # continuous percentage (0–100% or 0–1, depending on the source raster).
    # --------------------------------------------------------------------------
    ## --------------------------------------------------------
    ## Extract open-lowland percentage at unique ASO locations
    ## --------------------------------------------------------
    aso_locations <- aso_plant_richness %>%
      dplyr::distinct(geometry)
    
    thresholdVals <- terra::extract(
      openlowland_boundary,
      terra::vect(aso_locations)
    )
    
    thresholdVal <- quantile(thresholdVals$lavland[thresholdVals$lavland > 0], probs = 0.25)
    
    openlowland_projected <- terra::project(
      openlowland_boundary,
      allCalibrationRaster[[1]],
      method = "bilinear"
    )
    
    
    # --------------------------------------------------------------------------
    # 5.3 Apply the open-lowland threshold
    # --------------------------------------------------------------------------
    #
    # A cell is considered open lowland when its percentage of open lowland is
    # greater than or equal to `openlowland_threshold`.
    #
    # The result is a logical raster:
    #
    #       TRUE  = meets the threshold
    #       FALSE = does not meet the threshold
    # --------------------------------------------------------------------------
    
    openlowland_mask <- (
      openlowland_projected >= thresholdVal
    )
    
    
    # --------------------------------------------------------------------------
    # 5.4 Save the logical mask
    # --------------------------------------------------------------------------
    
    terra::writeRaster(
      openlowland_mask,
      openlowland_mask_file,
      overwrite = FALSE
    )
    
    
  } else {
    
    
    # --------------------------------------------------------------------------
    # 5.5 Load the existing threshold-based mask
    # --------------------------------------------------------------------------
    
    openlowland_mask <- terra::rast(
      openlowland_mask_file
    )
  }
  
  
  # ============================================================================
  # 6. Convert the logical mask to a binary raster
  # ============================================================================
  #
  # The final mask uses:
  #
  #       1  = open lowland
  #       NA = not open lowland
  #
  # This format is convenient for masking subsequent spatial analyses.
  # ============================================================================
  
  openlowland_binary <- terra::ifel(
    openlowland_mask == 1,
    1,
    NA
  )
}


# ==============================================================================
# 7. Plot and save the open-lowland mask
# ==============================================================================
#
# The same mask is exported in both PDF and PNG format.
#
# The plot is useful for visually checking whether the resulting open-lowland
# boundary is spatially sensible.
# ==============================================================================


# ------------------------------------------------------------------------------
# 7.1 Define the output plot filename
# ------------------------------------------------------------------------------

plot_file <- file.path(
  resultFolder,
  "figures",
  "openlowland_mask"
)


# ------------------------------------------------------------------------------
# 7.2 Save PDF version
# ------------------------------------------------------------------------------

pdf(
  paste0(plot_file, ".pdf"),
  width = 8,
  height = 7
)

terra::plot(
  openlowland_binary,
  col = "black",
  legend = FALSE,
  axes = TRUE,
  main = "Open-lowland area"
)

dev.off()


# ------------------------------------------------------------------------------
# 7.3 Save high-resolution PNG version
# ------------------------------------------------------------------------------

png(
  filename = paste0(plot_file, ".png"),
  width = 8,
  height = 7,
  units = "in",
  res = 600
)

terra::plot(
  openlowland_binary,
  col = "black",
  alpha = 1,
  legend = FALSE,
  axes = TRUE,
  main = "Open-lowland area"
)

dev.off()
