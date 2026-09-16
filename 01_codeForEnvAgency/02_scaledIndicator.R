################################################################################
# Calculate and calibrate the insect indicator
#
# This script:
#
#   1. Loads the calibration-model and indicator-calculation functions.
#
#   2. Loads observed insect richness from the national insect monitoring data.
#
#   3. Defines the environmental and spatial predictors used in the calibration
#      model.
#
#   4. Combines the predictors into the final calibration predictor set.
#
#   5. Calculates the scaled insect indicator by comparing:
#
#          - interaction-supported insect predictions
#          - predicted insect occurrence
#          - observed insect richness
#          - environmental predictors
#          - reference semi-natural grassland locations
#          - nature-region boundaries
#          - the open-lowland study area
#
#   6. Saves the main calibration outputs and lookup tables for later use.
#
################################################################################


# ==============================================================================
# 1. Load the calibration model and indicator functions
# ==============================================================================
#
# `02_getCaliberationModel.R`
#     -> contains code/functions used to fit or retrieve the calibration model.
#
# `02_scaledIndicatorFunction.R`
#     -> contains `calculate_scaled_indicator()`, which performs the main
#        indicator calibration and scaling.
#
# ==============================================================================

source(
  "01_codeForEnvAgency/02_getCaliberationModel.R"
)

source(
  "01_codeForEnvAgency/02_scaledIndicatorFunction.R"
)


# ==============================================================================
# 2. Load observed insect richness
# ==============================================================================
#
# The national insect monitoring dataset is processed by:
#
#     04_process_national_insect_monitoring_data.R
#
# This script is expected to create the object:
#
#     national_monitoring_insect_richness
#
# which contains observed insect richness at the monitoring locations.
#
# `bind_rows()` is used here so that, if the source object contains multiple
# data frames (e.g. different monitoring groups or subsets), they are combined
# into one data frame.
# ==============================================================================

source(
  "01_codeForEnvAgency/04_process_national_insect_monitoring_data.R"
)


observed_richness <- dplyr::bind_rows(
  national_monitoring_insect_richness
)


# ==============================================================================
# 3. Define calibration predictors
# ==============================================================================
#
# The calibration model relates the spatially predicted insect indicator to
# observed insect richness using a set of environmental and spatial predictors.
#
# The predictors are divided into two groups:
#
#     1. CLC/elevation predictors
#     2. Other spatial predictors
#
# Keeping these groups separate makes it easier to modify the model and to
# distinguish land-cover-derived predictors from other environmental variables.
# ==============================================================================


# ==============================================================================
# 3A. CLC and elevation predictors
# ==============================================================================
#
# These variables represent differences in CLC-derived environmental
# characteristics, together with elevation.
#
# The `_scaled` suffix indicates that these predictors have already been
# standardised.
#
# For a standardised predictor:
#
#     scaled value = (x - mean(x)) / SD(x)
#
# This generally gives predictors approximately:
#
#     mean = 0
#     SD   = 1
#
# which can improve numerical stability when fitting regression models.
# ==============================================================================

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
  
  # --------------------------------------------------------------------------
  # Heterogeneity variables are currently excluded from the model.
  #
  # They can be reintroduced by removing the leading `#`.
  # --------------------------------------------------------------------------
  
  # "CLC_heterogeneity_scaled",
  # "CLCplus_heterogeneity_scaled",
  
  "elevation"
)


# ------------------------------------------------------------------------------
# Alternative: unscaled CLC predictors
#
# This block is retained as a reference but is not currently used.
#
# If activated, it would replace the standardised CLC difference variables
# with their unscaled versions.
# ------------------------------------------------------------------------------

# clc_vars <- c(
#   "CLC_1_diff",
#   "CLC_2_diff",
#   "CLC_3_diff",
#   "CLC_5_diff",
#   "CLC_6_diff",
#   "CLC_7_diff",
#   "CLC_9_diff",
#   "CLC_10_diff",
#   "CLC_11_diff",
#   "elevation"
# )


# ==============================================================================
# 3B. Other spatial predictors
# ==============================================================================
#
# These are additional spatial predictors that are not part of the CLC
# downscaling/calibration variables.
# ==============================================================================

other_vars <- c(
  
  # --------------------------------------------------------------------------
  # Difference in landscape heterogeneity.
  # --------------------------------------------------------------------------
  
  "landscape_heterogeneity_difference",
  
  # --------------------------------------------------------------------------
  # Plant expected richness is currently excluded from the calibration model.
  #
  # It is nevertheless passed separately to `calculate_scaled_indicator()`
  # below, so it can still be used by the indicator calculation if required.
  # --------------------------------------------------------------------------
  
  # "plant_expected_richness",
  
  # --------------------------------------------------------------------------
  # Observed/predicted insect richness predictor.
  # --------------------------------------------------------------------------
  
  "insect_richness"
)


# ==============================================================================
# 4. Combine all calibration predictors
# ==============================================================================
#
# `calibration_vars` contains the complete set of predictors that will be
# supplied to the calibration model.
#
# This combines:
#
#     CLC/elevation predictors
#                 +
#     other spatial predictors
#
# ==============================================================================

calibration_vars <- c(
  clc_vars,
  other_vars
)


# ==============================================================================
# 5. Calculate the scaled insect indicator
# ==============================================================================
#
# This is the main analysis step.
#
# `calculate_scaled_indicator()` combines the predicted insect/plant
# interactions with observed monitoring data and environmental information
# to produce the calibrated indicator.
#
# -------------------------------------------------------------------------------
# Main inputs
# -------------------------------------------------------------------------------
#
# richness:
#     Interaction-supported insect probability rasters.
#
#     `weighted_insects_probability` contains the insect predictions after
#     accounting for plant–insect interactions.
#
#
# insects:
#     The original insect prediction rasters.
#
#     These are used as the underlying insect probability/suitability surface.
#
#
# study_region:
#     Binary open-lowland mask defining the area over which the indicator
#     is calculated.
#
#
# classificationRegions:
#     The five broad nature regions used for regional classification:
#
#         South
#         East
#         West
#         Central
#         North
#
#
# referencePolygons:
#     High-quality semi-natural grassland NiN polygons used as reference
#     locations for calibration/scaling.
#
#
# observed_richness:
#     Observed insect richness from national monitoring data.
#
#
# allCalibrationRaster:
#     The spatial raster stack containing the environmental calibration
#     predictors.
#
#
# plant_expected_richness:
#     Expected plant richness at 50-m resolution.
#
#
# prediction_year:
#     The year associated with the spatial prediction.
#
#
# calibration_vars:
#     Complete set of variables included in the calibration model.
#
#
# clc_vars:
#     Subset of calibration variables associated with CLC/elevation.
#
#
# cores:
#     Number of processing cores used for the calculation.
#
#
# rerun:
#     Controls whether the calibration model should be rerun rather than
#     reusing an existing model/result.
#
# ==============================================================================

indicator_results <- calculate_scaled_indicator(
  
  # --------------------------------------------------------------------------
  # Interaction-supported insect predictions
  #
  # Only the insect species represented in `species_names` are retained.
  # --------------------------------------------------------------------------
  
  richness =
    weighted_insects_probability %>%
    tidyterra::select(
      dplyr::all_of(species_names)
    ),
  
  
  # --------------------------------------------------------------------------
  # Original insect predictions
  # --------------------------------------------------------------------------
  
  insects =
    insects %>%
    tidyterra::select(
      dplyr::all_of(species_names)
    ),
  
  
  # --------------------------------------------------------------------------
  # Spatial extent where the indicator is calculated
  #
  # 1 = open lowland
  # NA = outside the study region
  # --------------------------------------------------------------------------
  
  study_region =
    openlowland_binary,
  
  
  # --------------------------------------------------------------------------
  # Broad geographical/nature regions
  # --------------------------------------------------------------------------
  
  classificationRegions =
    nature_regions,
  
  
  # --------------------------------------------------------------------------
  # Reference habitat polygons
  #
  # These represent good-quality semi-natural grassland NiN polygons.
  # --------------------------------------------------------------------------
  
  referencePolygons =
    ref_locs,
  
  
  # --------------------------------------------------------------------------
  # Coordinate reference system
  # --------------------------------------------------------------------------
  
  crs =
    crs,
  
  
  # --------------------------------------------------------------------------
  # Directory in which indicator plots/results are stored
  # --------------------------------------------------------------------------
  
  resultFolder =
    plotFolder,
  
  
  # --------------------------------------------------------------------------
  # Observed insect richness from national monitoring
  # --------------------------------------------------------------------------
  
  observed_richness =
    observed_richness,
  
  
  # --------------------------------------------------------------------------
  # Raster stack containing the calibration predictors
  # --------------------------------------------------------------------------
  
  allCalibrationRaster =
    allCalibrationRaster,
  
  
  # --------------------------------------------------------------------------
  # Expected plant richness at 50-m resolution
  # --------------------------------------------------------------------------
  
  plant_expected_richness =
    plant_expected_richness_50m,
  
  
  # --------------------------------------------------------------------------
  # Year for which the indicator is being predicted
  # --------------------------------------------------------------------------
  
  prediction_year =
    prediction_year,
  
  
  # --------------------------------------------------------------------------
  # Complete set of calibration predictors
  # --------------------------------------------------------------------------
  
  calibration_vars =
    calibration_vars,
  
  
  # --------------------------------------------------------------------------
  # CLC/elevation subset of calibration predictors
  # --------------------------------------------------------------------------
  
  clc_vars =
    clc_vars,
  
  
  # --------------------------------------------------------------------------
  # Use one processing core
  #
  # This is conservative and can reduce problems with raster file locking
  # and temporary files.
  # --------------------------------------------------------------------------
  
  cores = 1,
  
  
  # --------------------------------------------------------------------------
  # Whether to refit/recalculate the calibration model
  # --------------------------------------------------------------------------
  
  rerun =
    rerunCalibrationModel
)


# ==============================================================================
# 6. Save reference values
# ==============================================================================
#
# `reference_values` contains the values calculated from the reference
# locations/polygons.
#
# These values are important for scaling the indicator relative to the
# reference condition.
#
# They are saved separately so they can be reused without recalculating the
# complete indicator.
# ==============================================================================

saveRDS(
  indicator_results$reference_values,
  file = file.path(
    plotFolder,
    "insect_indicator_reference_values.rds"
  )
)


# ==============================================================================
# 7. Save the regional lookup table
# ==============================================================================
#
# `region_lookup` stores the relationship between spatial locations and their
# assigned nature region.
#
# This can be reused later when analysing or plotting the indicator by region.
# ==============================================================================

saveRDS(
  indicator_results$region_lookup,
  file = file.path(
    plotFolder,
    "insect_indicator_region_lookup.rds"
  )
)


# ==============================================================================
# 8. Save the complete indicator results
# ==============================================================================
#
# The complete `indicator_results` object is saved as an RData file.
#
# Unlike the two RDS files above, this preserves the entire list/object
# returned by `calculate_scaled_indicator()`.
#
# It can subsequently be loaded using:
#
#     load("indicator_results.RData")
#
# ==============================================================================

save(
  indicator_results,
  file = file.path(
    plotFolder,
    "indicator_results.RData"
  )
)
