################################################################################
# Prepare calibration predictor rasters at 50-m resolution
#
# This script:
#
#   1. Loads the land-cover and elevation predictor rasters.
#   2. Creates a common 50-m raster grid.
#   3. Resamples the calibration predictors to this 50-m grid.
#   4. Resamples the elevation-difference raster to the same grid.
#   5. Reprojects the elevation raster so that it exactly matches the
#      calibration predictor raster.
#   6. Scales the elevation predictor.
#   7. Combines all predictors into a single multi-layer raster.
#   8. Saves the final predictor stack as a GeoTIFF.
#
# Output:
#
#   calibration_preds_50m_scaled.tif
#
# The final raster contains all calibration predictors at 50-m resolution,
# including the scaled elevation predictor.
################################################################################


# ==============================================================================
# 0. Define the output file
# ==============================================================================

calibration_pred_file <- file.path(
  resultFolder,
  "polygons",
  "calibration_preds_50m_scaled.tif"
)


# ==============================================================================
# 1. Create the calibration predictor raster if it does not already exist
# ==============================================================================

if (!file.exists(calibration_pred_file)) {
  
  
  # ============================================================================
  # 2. Load the environmental predictor rasters
  # ============================================================================
  #
  # 00A_get_CLC_raster.R
  #     -> loads the Corine Land Cover / land-cover predictors
  #
  # 00C_get_elevation_raster.R
  #     -> creates/loads the elevation-difference raster
  #
  # These scripts are expected to create the objects:
  #
  #     allCalibrationRaster
  #     elevation_difference
  # ============================================================================
  
  source("01_codeForEnvAgency/00A_get_CLC_raster.R")
  
  source("01_codeForEnvAgency/00C_get_elevation_raster.R")
  
  
  # ============================================================================
  # 3. Create a 50-m raster template
  # ============================================================================
  #
  # The template defines the common spatial grid that will be used for the
  # calibration predictors.
  #
  # The template inherits:
  #
  #     extent -> from allCalibrationRaster
  #     CRS    -> from allCalibrationRaster
  #
  # but sets the resolution explicitly to 50 × 50 m.
  # ============================================================================
  
  template_50m <- terra::rast(
    ext = terra::ext(allCalibrationRaster),
    resolution = 50,
    crs = terra::crs(allCalibrationRaster)
  )
  
  
  # ============================================================================
  # 4. Resample the calibration predictors to 50-m resolution
  # ============================================================================
  #
  # All layers in allCalibrationRaster are transformed onto the new 50-m grid.
  #
  # Bilinear interpolation is used here.
  #
  # IMPORTANT:
  # Bilinear interpolation is appropriate for continuous variables such as
  # temperature, precipitation, elevation, etc.
  #
  # If allCalibrationRaster contains categorical variables (e.g. land-cover
  # classes), nearest-neighbour ("near") resampling should be used for those
  # layers instead. Bilinear interpolation can create values that do not
  # represent valid categories.
  # ============================================================================
  
  allCalibrationRaster_50m <- terra::resample(
    allCalibrationRaster,
    template_50m,
    method = "bilinear"
  )
  
  
  # Optional visual check of the resampled predictors
  plot(allCalibrationRaster_50m)
  
  
  # ============================================================================
  # 5. Resample the elevation-difference raster to the 50-m template
  # ============================================================================
  #
  # First, the elevation-difference raster is resampled to the 50-m grid.
  #
  # Bilinear interpolation is appropriate because elevation difference is a
  # continuous numerical variable.
  # ============================================================================
  
  elevation_difference_50m <- terra::resample(
    elevation_difference,
    template_50m,
    method = "bilinear"
  )
  
  
  # ============================================================================
  # 6. Reproject elevation difference to match the calibration raster
  # ============================================================================
  #
  # The object name "_25833" suggests that the calibration raster uses
  # EPSG:25833 (ETRS89 / UTM zone 33N).
  #
  # Rather than relying on the CRS code in the object name, we use
  # allCalibrationRaster_50m itself as the target raster.
  #
  # This ensures that the elevation raster has:
  #
  #     - the same CRS
  #     - the same extent
  #     - the same resolution
  #     - the same cell alignment
  #
  # as the other calibration predictors.
  # ============================================================================
  
  elevation_difference_50m_25833 <- terra::project(
    elevation_difference_50m,
    allCalibrationRaster_50m,
    method = "bilinear"
  )
  
  
  # ============================================================================
  # 7. Scale the elevation-difference predictor
  # ============================================================================
  #
  # Scaling converts the elevation predictor to a standardized variable:
  #
  #       scaled value = (x - mean(x)) / SD(x)
  #
  # Therefore:
  #
  #       mean ≈ 0
  #       SD   ≈ 1
  #
  # This is useful when predictors have very different numerical ranges and
  # is particularly useful for regression / spatial statistical models.
  #
  # NOTE:
  # The mean and SD used here are calculated from the raster being scaled.
  # If this predictor is later used for prediction, the same scaling
  # parameters should be applied to the prediction raster rather than
  # recalculating mean and SD from the prediction data.
  # ============================================================================
  
  elevation_mean <- terra::global(
    elevation_difference_50m_25833,
    fun = "mean",
    na.rm = TRUE
  )[[1]]
  
  elevation_sd <- terra::global(
    elevation_difference_50m_25833,
    fun = "sd",
    na.rm = TRUE
  )[[1]]
  
  
  elevation_difference_50m_scaled <- (
    elevation_difference_50m_25833 - elevation_mean
  ) / elevation_sd
  
  
  # ============================================================================
  # 8. Combine all calibration predictors
  # ============================================================================
  #
  # The calibration predictors and the scaled elevation-difference raster are
  # combined into one multi-layer SpatRaster.
  #
  # The resulting object contains:
  #
  #     allCalibrationRaster_50m
  #     +
  #     elevation_difference_50m_scaled
  #
  # as separate raster layers.
  # ============================================================================
  
  all_calibration_pred <- c(
    allCalibrationRaster_50m,
    elevation_difference_50m_scaled
  )
  
  
  # ============================================================================
  # 9. Save the complete predictor stack
  # ============================================================================
  #
  # FLT4S = 32-bit floating point
  # LZW   = lossless GeoTIFF compression
  # ============================================================================
  
  terra::writeRaster(
    all_calibration_pred,
    filename = calibration_pred_file,
    overwrite = TRUE,
    wopt = list(
      datatype = "FLT4S",
      gdal = "COMPRESS=LZW"
    )
  )
}


# ==============================================================================
# 10. Load the final calibration predictor raster
# ==============================================================================

allCalibrationRaster <- terra::rast(
  calibration_pred_file
)


# ==============================================================================
# 11. Rename the elevation layer
# ==============================================================================
#
# The elevation-difference layer is the 13th layer in the predictor stack.
#
# Rename it to "elevation" so that it can be referred to explicitly in
# downstream modelling code.
# ==============================================================================

names(allCalibrationRaster)[13] <- "elevation"
