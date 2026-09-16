################################################################################
# CLC calibration / downscaling covariates
#
# This script:
#
#   1. Loads CLC100 and CLC+10 1-km class-frequency rasters.
#   2. Loads landscape heterogeneity rasters for both map products.
#   3. Harmonises CLC100 classes to the CLC+ backbone classification.
#   4. Aggregates CLC100 class frequencies to the CLC+ classes.
#   5. Calculates CLC+10 - CLC100 differences for common classes.
#   6. Standardises landscape heterogeneity variables.
#   7. Standardises the CLC class-difference variables.
#   8. Calculates the difference in standardised landscape heterogeneity.
#   9. Combines all modelling covariates.
#  10. Writes the final raster to disk.
#
# Interpretation of CLC difference variables:
#
#   CLC_diffs = CLC+10 - CLC100
#
#   Positive values:
#       CLC+10 has a greater proportion of the class.
#
#   Negative values:
#       CLC100 has a greater proportion of the class.
#
# The differences are subsequently standardised for modelling.
################################################################################


## ============================================================================
## 1. Packages
## ============================================================================

library(terra)

if(!file.exists(file.path(
  resultFolder,
  "polygons",
  "CLC_calibration_downscaling_covariates_scaled.tif"
))){
## ============================================================================
## 2. Project directory
##
## Use the Windows short path to avoid problems with the Norwegian "å" in
## "åpent" when accessed through R on Windows/OneDrive.
## ============================================================================

project_dir <- "C:/Users/KWAKU~1.ADJ/ONEDRI~1/222608~1"


## ============================================================================
## 3. Input files
## ============================================================================

CLC100_file <- file.path(
  project_dir,
  "r CLC100 1km class freq.tif"
)

CLCplus10_file <- file.path(
  project_dir,
  "r CLC10m 1km class freq.tif"
)

CLC100_heterogeneity_file <- file.path(
  project_dir,
  "r CLC100 1km Landscape.heterogeneity.tif"
)

CLCplus_heterogeneity_file <- file.path(
  project_dir,
  "r CLC10Plus 1km Landscape.heterogeneity.tif"
)


## ============================================================================
## 4. Check input files
## ============================================================================

input_files <- c(
  CLC100_file,
  CLCplus10_file,
  CLC100_heterogeneity_file,
  CLCplus_heterogeneity_file
)

missing_files <- input_files[!file.exists(input_files)]

if (length(missing_files) > 0L) {
  stop(
    "The following input files could not be found:\n",
    paste(missing_files, collapse = "\n")
  )
}


## ============================================================================
## 5. Load rasters
## ============================================================================

r.CLC100 <- terra::rast(CLC100_file)

r.CLC10 <- terra::rast(CLCplus10_file)

CLC_heterogeneity <- terra::rast(
  CLC100_heterogeneity_file
)

CLCplus_heterogeneity <- terra::rast(
  CLCplus_heterogeneity_file
)


## ============================================================================
## 6. Name landscape heterogeneity layers
## ============================================================================

names(CLC_heterogeneity) <- "CLC_heterogeneity"

names(CLCplus_heterogeneity) <- "CLCplus_heterogeneity"


## ============================================================================
## 7. Inspect input rasters
## ============================================================================

print(r.CLC100)
print(r.CLC10)

print(CLC_heterogeneity)
print(CLCplus_heterogeneity)


## ============================================================================
## 8. Check raster geometry
##
## CLC100 and CLC+10 must be on the same grid because their class frequencies
## are compared pixel-by-pixel.
## ============================================================================

if (!terra::compareGeom(
  r.CLC100,
  r.CLC10,
  stopOnError = FALSE
)) {
  
  stop(
    "CLC100 and CLC+10 rasters do not have identical geometry. ",
    "Check CRS, extent, resolution and origin before calculating differences."
  )
}


## ============================================================================
## 9. Check landscape heterogeneity geometry
## ============================================================================

if (!terra::compareGeom(
  CLC_heterogeneity,
  CLCplus_heterogeneity,
  stopOnError = FALSE
)) {
  
  stop(
    "CLC100 and CLC+10 landscape heterogeneity rasters do not have ",
    "identical geometry."
  )
}


## ============================================================================
## 10. CLC100 -> CLC+ crosswalk
##
## CLC100 classes are mapped to the corresponding CLC+ backbone class.
##
## CLC+ backbone:
##
##   1  Sealed
##   2  Woody - needle leaved trees
##   3  Woody - Broadleaved deciduous trees
##   4  Woody - Broadleaved evergreen trees
##   5  Low-growing woody plants
##   6  Permanent herbaceous
##   7  Periodically herbaceous
##   8  Lichens and mosses
##   9  Non- and sparsely-vegetated
##   10 Water
##   11 Snow and ice
##
## CLC100 mixed forest (25) is not assigned because there is no defensible
## allocation among the CLC+ forest classes in the current crosswalk.
## ============================================================================

CLC100_to_CLCplus <- c(
  
  ## --------------------------------------------------------------------------
  ## Sealed
  ## --------------------------------------------------------------------------
  
  `1` = 1,
  `2` = 1,
  `3` = 1,
  `4` = 1,
  `5` = 1,
  `6` = 1,
  
  
  ## --------------------------------------------------------------------------
  ## Non- and sparsely-vegetated
  ## --------------------------------------------------------------------------
  
  `7` = 9,
  `8` = 9,
  `9` = 9,
  
  
  ## --------------------------------------------------------------------------
  ## Permanent herbaceous
  ## --------------------------------------------------------------------------
  
  `10` = 6,
  `11` = 6,
  `18` = 6,
  
  
  ## --------------------------------------------------------------------------
  ## Periodically herbaceous
  ## --------------------------------------------------------------------------
  
  `12` = 7,
  `20` = 7,
  
  
  ## --------------------------------------------------------------------------
  ## Agriculture / herbaceous vegetation
  ## --------------------------------------------------------------------------
  
  `21` = 6,
  
  
  ## --------------------------------------------------------------------------
  ## Forest
  ## --------------------------------------------------------------------------
  
  `23` = 3,
  `24` = 2,
  
  ## Mixed forest deliberately left unassigned
  `25` = NA,
  
  
  ## --------------------------------------------------------------------------
  ## Permanent herbaceous
  ## --------------------------------------------------------------------------
  
  `26` = 6,
  
  
  ## --------------------------------------------------------------------------
  ## Low-growing woody plants
  ## --------------------------------------------------------------------------
  
  `27` = 5,
  `29` = 5,
  
  
  ## --------------------------------------------------------------------------
  ## Non- and sparsely-vegetated
  ## --------------------------------------------------------------------------
  
  `30` = 9,
  `31` = 9,
  `32` = 9,
  `33` = 9,
  
  
  ## --------------------------------------------------------------------------
  ## Snow and ice
  ## --------------------------------------------------------------------------
  
  `34` = 11,
  
  
  ## --------------------------------------------------------------------------
  ## Permanent herbaceous
  ## --------------------------------------------------------------------------
  
  `35` = 6,
  `36` = 6,
  
  
  ## --------------------------------------------------------------------------
  ## Non- and sparsely-vegetated
  ## --------------------------------------------------------------------------
  
  `39` = 9,
  
  
  ## --------------------------------------------------------------------------
  ## Water
  ## --------------------------------------------------------------------------
  
  `40` = 10,
  `41` = 10,
  `43` = 10,
  `44` = 10
)


## ============================================================================
## 11. Extract CLC100 class numbers
## ============================================================================

CLC100_classes <- as.integer(
  sub(
    "^CLC100_",
    "",
    names(r.CLC100)
  )
)


## ============================================================================
## 12. Translate CLC100 classes to CLC+ classes
## ============================================================================

CLCplus_classes_from_CLC100 <- unname(
  CLC100_to_CLCplus[
    as.character(CLC100_classes)
  ]
)


## ============================================================================
## 13. Identify CLC100 layers that can be harmonised
## ============================================================================

keep <- !is.na(
  CLCplus_classes_from_CLC100
)

target_classes <- sort(
  unique(
    CLCplus_classes_from_CLC100[keep]
  )
)


## ============================================================================
## 14. Aggregate CLC100 frequencies to the CLC+ backbone
##
## Because the input raster contains class frequencies within each 1-km cell,
## aggregation is performed by summing the frequencies of all CLC100 classes
## belonging to the same CLC+ class.
## ============================================================================

r.CLC100.to.CLCplus <- terra::tapp(
  r.CLC100[[which(keep)]],
  index = CLCplus_classes_from_CLC100[keep],
  fun = sum,
  na.rm = TRUE
)


## ============================================================================
## 15. Name aggregated CLC100 layers
## ============================================================================

names(r.CLC100.to.CLCplus) <- paste0(
  "CLCplus_",
  target_classes
)


## ============================================================================
## 16. Extract CLC+10 class numbers
## ============================================================================

CLC10_classes <- as.integer(
  sub(
    "^CLCPlus10_",
    "",
    names(r.CLC10)
  )
)


## ============================================================================
## 17. Identify common CLC+ classes
##
## Only classes represented in both:
##
##   - CLC+10
##   - harmonised CLC100
##
## can be directly compared.
## ============================================================================

common_classes <- intersect(
  CLC10_classes,
  target_classes
)


if (length(common_classes) == 0L) {
  
  stop(
    "No common CLC+ classes were found between CLC+10 and the ",
    "harmonised CLC100 raster."
  )
}


## ============================================================================
## 18. Match layer positions
## ============================================================================

CLC10_idx <- match(
  common_classes,
  CLC10_classes
)

CLC100_idx <- match(
  common_classes,
  target_classes
)


## ============================================================================
## 19. Calculate CLC class differences
##
## Difference = CLC+10 - CLC100
## ============================================================================

CLC_diffs <- (
  r.CLC10[[CLC10_idx]] -
    r.CLC100.to.CLCplus[[CLC100_idx]]
)


## ============================================================================
## 20. Name raw difference layers
## ============================================================================

names(CLC_diffs) <- paste0(
  "CLC_",
  common_classes,
  "_diff"
)


## ============================================================================
## 21. Function for raster standardisation
##
## Each raster layer is transformed as:
##
##              x - mean(x)
##      z = -------------------
##              sd(x)
##
## This makes each modelling covariate approximately mean = 0 and SD = 1.
## ============================================================================

scale_raster <- function(x) {
  
  ## Calculate mean for each layer
  raster_mean <- terra::global(
    x,
    fun = "mean",
    na.rm = TRUE
  )[, 1]
  
  
  ## Calculate standard deviation for each layer
  raster_sd <- terra::global(
    x,
    fun = "sd",
    na.rm = TRUE
  )[, 1]
  
  
  ## Check for layers with zero variance
  if (any(
    is.na(raster_sd) |
    raster_sd == 0
  )) {
    
    bad_layers <- names(x)[
      is.na(raster_sd) |
        raster_sd == 0
    ]
    
    stop(
      "Cannot standardise raster layer(s) with zero/undefined SD: ",
      paste(bad_layers, collapse = ", ")
    )
  }
  
  
  ## Standardise each layer
  x_scaled <- (x - raster_mean) / raster_sd
  
  
  ## Preserve original layer names
  names(x_scaled) <- names(x)
  
  
  return(x_scaled)
}


## ============================================================================
## 22. Scale CLC100 landscape heterogeneity
## ============================================================================

CLC_heterogeneity_scaled <- scale_raster(
  CLC_heterogeneity
)

names(CLC_heterogeneity_scaled) <- (
  "CLC_heterogeneity_scaled"
)


## ============================================================================
## 23. Scale CLC+10 landscape heterogeneity
## ============================================================================

CLCplus_heterogeneity_scaled <- scale_raster(
  CLCplus_heterogeneity
)

names(CLCplus_heterogeneity_scaled) <- (
  "CLCplus_heterogeneity_scaled"
)


## ============================================================================
## 24. Calculate difference in standardised landscape heterogeneity
##
## This represents the difference between the two map products after
## standardising their respective heterogeneity measures.
## ============================================================================

landscape_heterogeneity_difference <- (
  CLCplus_heterogeneity_scaled -
    CLC_heterogeneity_scaled
)

names(
  landscape_heterogeneity_difference
) <- "landscape_heterogeneity_difference"


## ============================================================================
## 25. Scale CLC class differences
##
## The raw CLC differences remain available as `CLC_diffs`.
##
## For modelling, however, each difference is standardised to mean = 0 and
## SD = 1.
## ============================================================================

CLC_diffs_scaled <- scale_raster(
  CLC_diffs
)


## ============================================================================
## 26. Rename scaled CLC difference layers
## ============================================================================

names(CLC_diffs_scaled) <- paste0(
  names(CLC_diffs),
  "_scaled"
)


## ============================================================================
## 27. Combine all modelling covariates
##
## Final layers:
##
##   1. CLC100 landscape heterogeneity
##   2. CLC+10 landscape heterogeneity
##   3. Difference in landscape heterogeneity
##   4+. Standardised CLC class differences
## ============================================================================

allCalibrationRaster <- c(
  CLC_heterogeneity_scaled,
  CLCplus_heterogeneity_scaled,
  landscape_heterogeneity_difference,
  CLC_diffs_scaled
)

allCalibrationRaster_withoutScale <- c(
  CLC_heterogeneity,
  CLCplus_heterogeneity,
  landscape_heterogeneity_difference,
  CLC_diffs
)


## ============================================================================
## 28. Check final raster
## ============================================================================

print(allCalibrationRaster)

print(
  names(allCalibrationRaster)
)


## ============================================================================
## 29. Check dimensions and geometry
## ============================================================================

cat("\nNumber of modelling covariates:", nlyr(allCalibrationRaster), "\n")

cat(
  "Rows:",
  nrow(allCalibrationRaster),
  "\n"
)

cat(
  "Columns:",
  ncol(allCalibrationRaster),
  "\n"
)

cat(
  "Resolution:",
  paste(res(allCalibrationRaster), collapse = " x "),
  "\n"
)

cat(
  "CRS:\n",
  crs(allCalibrationRaster),
  "\n"
)


## ============================================================================
## 30. Check summary statistics
##
## This provides a quick check that the scaled variables have approximately
## mean = 0 and SD = 1.
## ============================================================================

calibration_summary <- data.frame(
  variable = names(allCalibrationRaster),
  mean = terra::global(
    allCalibrationRaster,
    "mean",
    na.rm = TRUE
  )[, 1],
  sd = terra::global(
    allCalibrationRaster,
    "sd",
    na.rm = TRUE
  )[, 1],
  min = terra::global(
    allCalibrationRaster,
    "min",
    na.rm = TRUE
  )[, 1],
  max = terra::global(
    allCalibrationRaster,
    "max",
    na.rm = TRUE
  )[, 1]
)

print(calibration_summary)


## ============================================================================
## 31. Output file
## ============================================================================

output_file <- file.path(
  resultFolder,
  "polygons",
  "CLC_calibration_downscaling_covariates_scaled.tif"
)



## ============================================================================
## 32. Write final modelling raster
## ============================================================================

target_crs <- "EPSG:25833"

allCalibrationRaster_25833 <- terra::project(
  allCalibrationRaster,
  target_crs,
  method = "bilinear"
)

allCalibrationRaster_withoutScale_25833 <- terra::project(
  allCalibrationRaster_withoutScale,
  target_crs,
  method = "bilinear"
)

terra::writeRaster(
  allCalibrationRaster_25833,
  output_file,
  overwrite = TRUE
)

terra::writeRaster(
  allCalibrationRaster_withoutScale_25833,
  file.path(
    resultFolder,
    "polygons",
    "CLC_calibration_downscaling_covariates.tif"
  ),
  overwrite = TRUE
)


## ============================================================================
## 33. Confirm output
## ============================================================================

cat(
  "\nFinal calibration/downscaling raster written to:\n",
  output_file,
  "\n"
)

cat(
  "\nFile exists:",
  file.exists(output_file),
  "\n"
)


## ============================================================================
## 34. Optional: save the raw CLC differences separately
##
## This is useful because the raw differences retain their direct ecological
## interpretation and can be used for maps/reporting.
## ============================================================================

raw_difference_file <- file.path(
  resultFolder,
  "polygons",
  "CLC_class_differences_raw.tif"
)

terra::writeRaster(
  CLC_diffs,
  raw_difference_file,
  overwrite = TRUE
)


## ============================================================================
## 35. Optional: save the harmonised CLC100 raster
##
## This allows the intermediate harmonisation product to be inspected later.
## ============================================================================

harmonised_CLC100_file <- file.path(
  resultFolder,
  "polygons",
  "CLC100_harmonised_to_CLCplus.tif"
)

terra::writeRaster(
  r.CLC100.to.CLCplus,
  harmonised_CLC100_file,
  overwrite = TRUE
)


## ============================================================================
## 36. Final objects available in the R session
##
## allCalibrationRaster
##     Final standardised modelling covariates.
##
## CLC_diffs
##     Raw CLC+10 - CLC100 class differences.
##
## CLC_diffs_scaled
##     Standardised CLC class differences used for modelling.
##
## r.CLC100.to.CLCplus
##     CLC100 class frequencies harmonised to the CLC+ backbone.
##
## calibration_summary
##     Summary statistics for the final modelling covariates.
## ============================================================================
} else {
  allCalibrationRaster <- terra::rast(file.path(
    resultFolder,
    "polygons",
    "CLC_calibration_downscaling_covariates_scaled.tif"
    #"CLC_calibration_downscaling_covariates.tif"
  ))
}

################################################################################
# END
################################################################################