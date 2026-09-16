################################################################################
# Calculate local elevation difference
#
# This script:
#
#   1. Creates a merged 50-m Digital Terrain Model (DTM) from the original
#      Kartverket DTM50 tiles, if the merged raster does not already exist.
#
#   2. Aggregates the 50-m DTM to 1-km cells using the mean elevation.
#
#   3. Resamples the 1-km mean elevation back onto the original 50-m grid.
#
#   4. Calculates the difference between:
#         - the elevation at each 50-m cell, and
#         - the mean elevation of its corresponding 1-km cell.
#
#      Positive values = location is higher than the surrounding 1-km mean.
#      Negative values = location is lower than the surrounding 1-km mean.
#
#   5. Saves the resulting elevation-difference raster.
#
################################################################################


# ==============================================================================
# 0. Define output file paths
# ==============================================================================

merged_dtm_file <- file.path(
  resultFolder,
  "DTM50_UTM33_merged.tif"
)

elevation_difference_file <- file.path(
  resultFolder,
  "elevation_UTM33_difference_50m_1km.tif"
)


# ==============================================================================
# 1. Create the merged 50-m DTM if it does not already exist
# ==============================================================================

if (!file.exists(elevation_difference_file)) {
  
  # --------------------------------------------------------------------------
  # 1.1 Build the merged DTM from the original Kartverket tiles
  #
  # The merge is only performed if the physical GeoTIFF does not already exist.
  # --------------------------------------------------------------------------
  
  if (!file.exists(merged_dtm_file)) {
    
    # Location of the original DTM50 tiles
    nin_gdb <- paste0(
      "R:/GeoSpatialData/Elevation/",
      "Norway_DEM_50m_Kartverket/Original/",
      "DTM50_UTM33_20210416"
    )
    
    
    # ------------------------------------------------------------------------
    # Find all GeoTIFF tiles
    # ------------------------------------------------------------------------
    
    dtm_tiles <- list.files(
      nin_gdb,
      pattern = "\\.tif$",
      full.names = TRUE,
      ignore.case = TRUE
    )
    
    # Check how many tiles were found
    length(dtm_tiles)
    
    
    # ------------------------------------------------------------------------
    # Create a temporary directory for the VRT
    #
    # A VRT (Virtual Raster) references the original tiles without creating
    # another full copy of the data.
    # ------------------------------------------------------------------------
    
    tmp_dtm <- file.path(
      "C:/terra_tmp",
      "dtm50_UTM33_tiles"
    )
    
    dir.create(
      tmp_dtm,
      recursive = TRUE,
      showWarnings = FALSE
    )
    
    
    # Path for the virtual mosaic
    dtm50_vrt <- file.path(
      tmp_dtm,
      "DTM50_UTM33_mosaic.vrt"
    )
    
    
    # ------------------------------------------------------------------------
    # Create a virtual mosaic of all DTM tiles
    # ------------------------------------------------------------------------
    
    dtm50 <- terra::vrt(
      dtm_tiles,
      filename = dtm50_vrt,
      overwrite = TRUE
    )
    
    dtm50
    
    
    # ------------------------------------------------------------------------
    # Write the virtual mosaic to a physical GeoTIFF
    #
    # This creates:
    #
    #   DTM50_UTM33_merged.tif
    #
    # FLT4S = 32-bit floating-point values
    # LZW   = lossless compression
    # ------------------------------------------------------------------------
    
    terra::writeRaster(
      dtm50,
      filename = merged_dtm_file,
      overwrite = TRUE,
      wopt = list(
        datatype = "FLT4S",
        gdal = "COMPRESS=LZW"
      )
    )
  }
  
  
  # =============================================================================
  # 2. Load the merged 50-m DTM
  # =============================================================================
  
  dtm_50m <- terra::rast(
    merged_dtm_file
  )
  
  
  # =============================================================================
  # 3. Calculate mean elevation at 1-km resolution
  # =============================================================================
  #
  # The original DTM has a 50-m resolution.
  #
  # Therefore:
  #
  #     1,000 m / 50 m = 20 cells
  #
  # So an aggregation factor of 20 converts the 50-m raster into a 1-km
  # raster.
  #
  # Each 1-km cell contains:
  #
  #     20 × 20 = 400
  #
  # original 50-m cells.
  #
  # The value of each 1-km cell is the mean elevation of those 400 cells.
  # =============================================================================
  
  dtm_1km <- terra::aggregate(
    dtm_50m,
    fact = 20,
    fun = mean,
    na.rm = TRUE
  )
  
  
  # =============================================================================
  # 4. Resample the 1-km mean elevation back onto the 50-m grid
  # =============================================================================
  #
  # The 1-km raster now contains the mean elevation for each 1-km area.
  #
  # We resample it back onto the original 50-m grid so that every 50-m cell
  # receives the mean elevation of the 1-km cell in which it occurs.
  #
  # "near" (nearest neighbour) is appropriate here because we are assigning
  # categorical blocks of 1-km mean values back to the finer grid rather than
  # interpolating a continuous surface.
  # =============================================================================
  
  dtm_1km_50m <- terra::resample(
    dtm_1km,
    dtm_50m,
    method = "near"
  )
  
  
  # =============================================================================
  # 5. Calculate the elevation difference
  # =============================================================================
  #
  # For every 50-m cell:
  #
  #     elevation difference =
  #         local 50-m elevation
  #         -
  #         mean elevation of corresponding 1-km cell
  #
  # Interpretation:
  #
  #   > 0  = cell is higher than the 1-km mean
  #   < 0  = cell is lower than the 1-km mean
  #   = 0  = cell is equal to the 1-km mean
  #
  # This therefore represents local topographic deviation from the broader
  # 1-km mean elevation.
  # =============================================================================
  
  elevation_difference <- dtm_50m - dtm_1km_50m
  
  
  # =============================================================================
  # 6. Save the elevation-difference raster
  # =============================================================================
  
  terra::writeRaster(
    elevation_difference,
    filename = elevation_difference_file,
    overwrite = TRUE,
    wopt = list(
      datatype = "FLT4S",
      gdal = "COMPRESS=LZW"
    )
  )
}


# ==============================================================================
# 7. Load the final elevation-difference raster
# ==============================================================================

elevation_difference <- terra::rast(
  elevation_difference_file
)

