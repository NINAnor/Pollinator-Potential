

library(terra)

if(!file.exists(file.path(
  resultFolder,
  "elevation_difference_10m_1km.tif"
))){

if(!file.exits(file.path(
  resultFolder,
  "DTM50_merged.tif"
))){
# ------------------------------------------------------------------
# 1. Find the ZIP files
# ------------------------------------------------------------------
# Path to elevation data
nin_gdb <- "R:/GeoSpatialData/Elevation/Norway_DEM_50m_Kartverket/Original/DTM50 TIFF-format"
zip_files <- list.files(
  nin_gdb,
  pattern = "\\.zip$",
  full.names = TRUE,
  ignore.case = TRUE
)

length(zip_files)

# -------------------------------------------------------------------
# 2. Create a temporary extraction directory
# ___________________________________________________________________
tmp_dtm <- file.path(
  "C:/terra_tmp",
  "dtm50_tiles"
)

dir.create(
  tmp_dtm,
  recursive = TRUE,
  showWarnings = FALSE
)

tmp_dtm

# -------------------------------------------------------------------
# 3. Unzip all files
# ___________________________________________________________________
for (z in zip_files) {
  
  # Name of ZIP without extension
  tile_name <- tools::file_path_sans_ext(
    basename(z)
  )
  
  # Directory for this tile
  tile_dir <- file.path(
    tmp_dtm,
    tile_name
  )
  
  dir.create(
    tile_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  # Extract
  unzip(
    z,
    exdir = tile_dir
  )
}

# ------------------------------------------------
# 4. Find all TIFFs
# ------------------------------------------------

tif_files <- list.files(
  tmp_dtm,
  pattern = "\\.tif(f)?$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

cat("Number of TIFF files:", length(tif_files), "\n")


# ------------------------------------------------
# 5. Create virtual mosaic
# ------------------------------------------------

dtm_vrt <- file.path(
  tmp_dtm,
  "DTM50_mosaic.vrt"
)

dtm <- terra::vrt(
  tif_files,
  filename = dtm_vrt,
  overwrite = TRUE
)


# ------------------------------------------------
# 6. Inspect
# ------------------------------------------------

print(dtm)

cat("CRS:\n")
print(crs(dtm))

cat("Resolution:\n")
print(res(dtm))

cat("Extent:\n")
print(ext(dtm))


# ------------------------------------------------
# 7. Optional: write physical GeoTIFF
# ------------------------------------------------

dtm_merged <- terra::writeRaster(
  dtm,
  filename = file.path(
    resultFolder,
    "DTM50_merged.tif"
  ),
  overwrite = TRUE,
  wopt = list(
    datatype = "FLT4S",
    gdal = c("COMPRESS=LZW")
  )
)

} 

dtm_merged <- terra::rast(file.path(
  resultFolder,
  "DTM50_merged.tif"
))

# 1. Aggregate 50-m DTM to 1-km mean
dtm_1km <- aggregate(
  dtm_merged,
  fact = 100,
  fun = mean,
  na.rm = TRUE
)

# 2. Put the 1-km mean onto the 50-m grid
dtm_1km_50m <- resample(
  dtm_1km,
  dtm_merged,
  method = "near"
)

# 3. Calculate elevation difference
elevation_difference <- dtm_merged - dtm_1km_50m

writeRaster(
  elevation_difference,
  filename = file.path(
    resultFolder,
    "elevation_difference_10m_1km.tif"
  ),
  overwrite = TRUE,
  wopt = list(
    datatype = "FLT4S",
    gdal = "COMPRESS=LZW"
  )
)
}

elevation_difference <- terra::rast(file.path(
  resultFolder,
  "elevation_difference_10m_1km.tif"
))
