
# ==============================================================================
# Calculate and calibrate the insect indicator
# ==============================================================================

# 1. Load functions and prepare observed data ----------------------------------

source("01_codeForEnvAgency/02_getCaliberationModel.R")
source("01_codeForEnvAgency/02_calculate_scaled_indicator.R")
source("01_codeForEnvAgency/02_calculate_reference_values.R")
source("01_codeForEnvAgency/04_process_national_insect_monitoring_data.R")

checks <- FALSE

observed_richness <- dplyr::bind_rows(national_monitoring_insect_richness)


# 2. Define calibration predictors ---------------------------------------------

# CLC differences and elevation
clc_vars <- c(
  "CLC_1_diff_scaled",
  "CLC_2_diff_scaled",
  "CLC_3_diff_scaled",
  "CLC_5_diff_scaled",
  "CLC_6_diff_scaled",
  "CLC_7_diff_scaled",
  "CLC_9_diff_scaled",
  "CLC_10_diff_scaled",
  "CLC_11_diff_scaled",
  "elevation"
)

# Other spatial predictors
other_vars <- c(
  "landscape_heterogeneity_difference",
  "insect_richness"
)

calibration_vars <- c(clc_vars, other_vars)


# 3. Calculate reference values and calibrated richness ------------------------

get_reference_vals <- calculate_reference_values(
  richness = weighted_insects_probability %>%
    tidyterra::select(dplyr::all_of(species_names)),
  insects = insects %>%
    tidyterra::select(dplyr::all_of(species_names)),
  study_region = openlowland_binary_with_survey_locs,
  classificationRegions = nature_regions,
  referencePolygons = ref_locs,
  crs = crs,
  resultFolder = plotFolder,
  observed_richness = observed_richness,
  allCalibrationRaster = allCalibrationRaster,
  plant_expected_richness = plant_expected_richness_50m,
  prediction_year = prediction_year,
  calibration_vars = calibration_vars,
  clc_vars = clc_vars,
  cores = 1,
  rerun = rerunCalibrationModel,
  use_max = use_max,
  overwrite = overwrite
)

allRasts <- get_reference_vals$reference_rasters$calibrated_insect_richness

reference_polygon_values <- get_reference_vals$reference_polygon_values %>%
  dplyr::select(
    nature_type,
    condition,
    nature_region,
    nature_type_code,
    calibrated_reference_richness
  ) %>%
  dplyr::filter(condition == "god") %>%
    dplyr::filter(nature_type_code %in% c("ntyp_D02_02", "ntyp_D02"))




# 4. Prepare reference values for indicator scaling ----------------------------

source("01_codeForEnvAgency/02_getReferenceVals.R", echo = TRUE)

reference_values$X_60 <- NULL  # two-anchor (X_0, X_100) scaling
# 5. Calculate the scaled insect indicator ------------------------------------

scaled_richness <- calculate_scaled_indicator(
  classificationRegions = nature_regions,
  overwrite = TRUE,
  rerun = FALSE,
  use_max = FALSE,
  allRasts = allRasts,
  ref_vals = reference_values,
  output_dir = plotFolder,
  value_layer = NULL,
  truncation = TRUE,
  cores = 1
)
