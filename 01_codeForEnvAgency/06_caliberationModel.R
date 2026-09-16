################################################################################
# INSECT RICHNESS CALIBRATION
#
# Objective:
#   Calibrate observed insect richness against CLC calibration/downscaling
#   variables, landscape heterogeneity, and expected plant richness.
#
# Candidate models:
#
#   A. Poisson GLMM - linear effects
#   B. NB2 GLMM    - linear effects
#   C. NB2 GLMM    - quadratic effects
#   D. NB2 GLM     - CLC effects differ among nature regions
#
# Validation:
#   80% training / 20% validation within each nature region
#
# Important:
#   All continuous raster predictors are standardised before model fitting.
#   The same scaling is applied to the corresponding prediction rasters.
################################################################################


################################################################################
# 0. Packages
################################################################################

library(dplyr)
library(tidyr)
library(sf)
library(terra)
library(glmmTMB)
library(ggplot2)
library(purrr)


################################################################################
# 1. Load observed insect richness
################################################################################

source(
  "01_codeForEnvAgency/04_process_national_insect_monitoring_data.R"
)

observed_richness <- bind_rows(
  national_monitoring_insect_richness
)


################################################################################
# 2. Define predictors
################################################################################

# CLC calibration/downscaling variables

clc_vars <- c(
  "CLC_1_diff_scaled",
  "CLC_2_diff_scaled",
  "CLC_3_diff_scaled",
  "CLC_5_diff_scaled",
  "CLC_6_diff_scaled",
  "CLC_7_diff_scaled",
  "CLC_9_diff_scaled",
  "CLC_10_diff_scaled",
  "CLC_11_diff_scaled"
)


# Other spatial predictors

other_vars <- c(
  "landscape_heterogeneity_difference",
  "plant_expected_richness",
  "richness"
)


# All predictors used in the spatial calibration

calibration_vars <- c(
  clc_vars,
  other_vars
)


################################################################################
# 3. Prepare plant expected richness raster
################################################################################

# The CLC calibration raster is in EPSG:25833 with metre units.
#
# plant_expected_richness is already on approximately the same spatial grid,
# but is represented in km units.
#
# Project it onto the exact CLC calibration raster grid.

plant_expected_richness_prediction <- terra::project(
  plant_expected_richness,
  allCalibrationRaster[[1]],
  method = "bilinear"
)%>%
  scale()

names(plant_expected_richness_prediction) <-
  "plant_expected_richness"

insect_expect_richness_prediction <- indicator_results$allRasts %>%
  tidyterra::select(richness,
                    richness_without_interaction) %>%
  terra::project(
    .,
    allCalibrationRaster[[1]],
    method = "bilinear"
  ) %>%
  scale()


################################################################################
# 4. Scale ALL continuous spatial predictors
################################################################################
#
# The CLC calibration variables are already scaled in
# allCalibrationRaster.
#
# Here we standardise plant_expected_richness in exactly the same way that
# the richness raster was standardised:
#
#       (x - mean(x)) / sd(x)
#
# The scaling parameters are calculated from the raster itself.
#
# These same parameters are therefore used for the raster values extracted
# at the calibration locations AND for the final spatial prediction.
################################################################################


scale_raster <- function(x) {
  
  r_mean <- terra::global(
    x,
    fun = "mean",
    na.rm = TRUE
  )[, 1]
  
  r_sd <- terra::global(
    x,
    fun = "sd",
    na.rm = TRUE
  )[, 1]
  
  if (
    any(
      is.na(r_sd) |
      r_sd == 0
    )
  ) {
    
    stop(
      "Cannot scale raster: at least one layer has ",
      "zero or undefined standard deviation."
    )
  }
  
  x_scaled <- (
    x - r_mean
  ) / r_sd
  
  names(x_scaled) <- names(x)
  
  x_scaled
}


# plant_expected_richness_scaled <- scale_raster(
#   plant_expected_richness_prediction
# )

# names(
#   plant_expected_richness_scaled
# ) <- "plant_expected_richness"


################################################################################
# 5. Create the complete spatial predictor raster
################################################################################

# Select only the CLC predictors required by the model.

clc_prediction_raster <- allCalibrationRaster[[
  c(
    clc_vars,
    "landscape_heterogeneity_difference"
  )
]]


# Add scaled expected plant richness.

calibration_prediction_raster <- c(
  clc_prediction_raster,
  plant_expected_richness_prediction,
  insect_expect_richness_prediction %>% tidyterra::select(richness)
)


################################################################################
# 6. Check predictor names
################################################################################

if (
  !identical(
    names(calibration_prediction_raster),
    calibration_vars
  )
) {
  
  stop(
    "Prediction raster layers do not match calibration_vars."
  )
}


print(
  names(calibration_prediction_raster)
)


################################################################################
# 7. Extract spatial predictors at monitoring locations
################################################################################

# One extraction for all calibration predictors.

extracted_predictors <- terra::extract(
  calibration_prediction_raster,
  terra::vect(observed_richness)
) %>%
  as.data.frame() %>%
  select(-ID)


################################################################################
# 8. Extract existing plant richness indicator
################################################################################

# The original richness raster was scaled before extraction.

# scaled_richness_raster <- indicator_results$allRasts[[
#   "richness"
# ]]
# 
# scaled_richness_raster <- scale_raster(
#   scaled_richness_raster
# )
# 
# names(scaled_richness_raster) <- "richness"
# 
# 
# extracted_richness <- terra::extract(
#   scaled_richness_raster,
#   terra::vect(observed_richness)
# ) %>%
#   as.data.frame() %>%
#   select(-ID)


################################################################################
# 9. Create calibration data
################################################################################

insectCalibration <- observed_richness %>%
  sf::st_drop_geometry() %>%
  select(
    insct_r,
    samp_dur,
    weight
  ) %>%
  # bind_cols(
  #   extracted_richness
  # ) %>%
  bind_cols(
    extracted_predictors
  )


################################################################################
# 10. Add nature regions
################################################################################

# Perform the spatial join directly on the original monitoring locations.

nature_region_values <- observed_richness %>%
  select(geometry) %>%
  sf::st_join(
    nature_regions %>%
      select(nature_region),
    join = sf::st_within,
    left = TRUE
  ) %>%
  st_drop_geometry() %>%
  select(nature_region)


insectCalibration <- insectCalibration %>%
  bind_cols(
    nature_region_values
  )


################################################################################
# 11. Prepare modelling data
################################################################################

model_data <- insectCalibration %>%
  select(
    insct_r,
    nature_region,
    all_of(calibration_vars),
    weight
  ) %>%
  mutate(
    nature_region = factor(
      nature_region
    )
  )


################################################################################
# 12. Remove incomplete observations
################################################################################

model_data <- model_data %>%
  filter(
    complete.cases(
      select(
        .,
        insct_r,
        nature_region,
        all_of(calibration_vars),
        weight
      )
    )
  )


################################################################################
# 13. Response diagnostics
################################################################################

cat("\n========================================\n")
cat("OBSERVED INSECT RICHNESS\n")
cat("========================================\n\n")

print(
  summary(model_data$insct_r)
)

cat("\nRange:\n")

print(
  range(
    model_data$insct_r,
    na.rm = TRUE
  )
)

cat("\nObservations by nature region:\n")

print(
  table(model_data$nature_region)
)


################################################################################
# 14. Check predictor scaling
################################################################################

cat("\n========================================\n")
cat("PREDICTOR RANGES\n")
cat("========================================\n\n")

print(
  model_data %>%
    summarise(
      across(
        all_of(calibration_vars),
        list(
          min = ~ min(.x, na.rm = TRUE),
          max = ~ max(.x, na.rm = TRUE),
          mean = ~ mean(.x, na.rm = TRUE),
          sd = ~ sd(.x, na.rm = TRUE)
        )
      )
    )
)


################################################################################
# 15. Create 80/20 train-validation split
################################################################################
#
# The split is stratified within nature region so that every region contributes
# observations to both training and validation.
################################################################################

set.seed(12345)

model_data <- model_data %>%
  group_by(nature_region) %>%
  mutate(
    validation_set = sample(
      c(
        rep(
          "train",
          floor(0.80 * n())
        ),
        rep(
          "validation",
          n() - floor(0.80 * n())
        )
      )
    )
  ) %>%
  ungroup()


################################################################################
# 16. Split data
################################################################################

train_data <- model_data %>%
  filter(
    validation_set == "train"
  ) %>%
  select(
    -validation_set
  )


validation_data <- model_data %>%
  filter(
    validation_set == "validation"
  ) %>%
  select(
    -validation_set
  )


################################################################################
# 17. Check split
################################################################################

cat("\nTraining observations by region:\n")

print(
  table(train_data$nature_region)
)

cat("\nValidation observations by region:\n")

print(
  table(validation_data$nature_region)
)


################################################################################
# 18. Construct model formulas
################################################################################


# ---------------------------------------------------------------------------
# Model A and B
# Linear effects
# ---------------------------------------------------------------------------

linear_formula <- as.formula(
  paste(
    "insct_r ~",
    paste(
      calibration_vars,
      collapse = " + "
    ),
    "+ (1 | nature_region)"
  )
)


# ---------------------------------------------------------------------------
# Model C
# Quadratic effects
# ---------------------------------------------------------------------------

polynomial_terms <- paste0(
  "poly(",
  calibration_vars,
  ", degree = 2, raw = TRUE)"
)


polynomial_formula <- as.formula(
  paste(
    "insct_r ~",
    paste(
      polynomial_terms,
      collapse = " + "
    ),
    "+ (1 | nature_region)"
  )
)


# ---------------------------------------------------------------------------
# Model D
# Nature-specific CLC effects
# ---------------------------------------------------------------------------

region_specific_formula <- as.formula(
  paste(
    "insct_r ~",
    paste(
      c(
        paste0(
          clc_vars,
          " * nature_region"
        ),
        other_vars
      ),
      collapse = " + "
    )
  )
)


################################################################################
# 19. Fit candidate models
################################################################################

cat("\n========================================\n")
cat("FITTING MODELS\n")
cat("========================================\n\n")


models <- list(
  
  "A. Poisson - linear" = glmmTMB(
    formula = linear_formula,
    data = train_data,
    family = poisson(link = "log")
  ),
  
  "B. NB2 - linear" = glmmTMB(
    formula = linear_formula,
    data = train_data,
    family = nbinom2(link = "log")
  ),
  
  "C. NB2 - polynomial" = glmmTMB(
    formula = polynomial_formula,
    data = train_data,
    family = nbinom2(link = "log")
  ),
  
  "D. NB2 - nature-specific CLC" = glmmTMB(
    formula = region_specific_formula,
    data = train_data,
    family = nbinom2(link = "log")
  )
)


################################################################################
# 20. AIC comparison
################################################################################

model_comparison_AIC <- purrr::map_dfr(
  models,
  ~ tibble(
    AIC = AIC(.x),
    BIC = BIC(.x),
    logLik = as.numeric(logLik(.x)),
    df = attr(logLik(.x), "df")
  ),
  .id = "model"
) %>%
  arrange(AIC)

print(model_comparison_AIC)


cat("\n========================================\n")
cat("AIC COMPARISON\n")
cat("========================================\n\n")

print(
  model_comparison_AIC
)


################################################################################
# 21. Predict validation data
################################################################################

validation_predictions <- purrr::map(
  models,
  ~ predict(
    .x,
    newdata = validation_data,
    type = "response",
    allow.new.levels = FALSE
  )
)


################################################################################
# 22. Add predictions to validation data
################################################################################

prediction_names <- c(
  "predicted_poisson",
  "predicted_nb",
  "predicted_polynomial",
  "predicted_region_specific"
)

names(validation_predictions) <- prediction_names

validation_data <- bind_cols(
  validation_data,
  as.data.frame(validation_predictions)
)


################################################################################
# 23. Validation statistics
################################################################################

calculate_validation_stats <- function(
    observed,
    predicted
) {
  
  keep <- complete.cases(
    observed,
    predicted
  )
  
  observed <- observed[keep]
  predicted <- predicted[keep]
  
  correlation <- cor(
    observed,
    predicted
  )
  
  tibble(
    n = length(observed),
    
    RMSE = sqrt(
      mean(
        (observed - predicted)^2
      )
    ),
    
    MAE = mean(
      abs(
        observed - predicted
      )
    ),
    
    bias = mean(
      predicted - observed
    ),
    
    correlation = correlation,
    
    R2 = correlation^2
  )
}


################################################################################
# 24. Compare predictive performance
################################################################################

validation_comparison <- purrr::map2_dfr(
  prediction_names,
  names(models),
  ~ calculate_validation_stats(
    observed = validation_data$insct_r,
    predicted = validation_data[[.x]]
  ) %>%
    mutate(
      model = .y,
      .before = 1
    )
) %>%
  arrange(RMSE)


cat("\n========================================\n")
cat("VALIDATION PERFORMANCE\n")
cat("========================================\n\n")

print(
  validation_comparison
)


################################################################################
# 25. Select best model
################################################################################

best_model_name <- validation_comparison %>%
  slice_min(
    RMSE,
    n = 1,
    with_ties = FALSE
  ) %>%
  pull(model)


insect_model <- models[[best_model_name]]


cat("\n========================================\n")
cat("SELECTED MODEL\n")
cat("========================================\n\n")

cat(
  best_model_name,
  "\n"
)


################################################################################
# 26. Identify selected prediction column
################################################################################

selected_prediction_column <- prediction_names[
  match(
    best_model_name,
    names(models)
  )
]


validation_data <- validation_data %>%
  mutate(
    predicted_insct_r =
      .data[[selected_prediction_column]]
  )


################################################################################
# 27. Validation performance by nature region
################################################################################

validation_by_region <- purrr::map2_dfr(
  prediction_names,
  names(models),
  function(pred_col, model_name) {
    
    validation_data %>%
      group_by(nature_region) %>%
      group_modify(
        function(region_data, region_key) {
          
          calculate_validation_stats(
            observed = region_data$insct_r,
            predicted = region_data[[pred_col]]
          )
        }
      ) %>%
      ungroup() %>%
      mutate(
        model = model_name,
        .before = 1
      )
  }
)


cat("\n========================================\n")
cat("VALIDATION BY NATURE REGION\n")
cat("========================================\n\n")

print(
  validation_by_region
)


################################################################################
# 28. Observed versus predicted richness
################################################################################

validation_long <- validation_data %>%
  select(
    insct_r,
    nature_region,
    all_of(prediction_names)
  ) %>%
  pivot_longer(
    cols = all_of(prediction_names),
    names_to = "prediction",
    values_to = "predicted_insct_r"
  ) %>%
  mutate(
    model = recode(
      prediction,
      
      predicted_poisson =
        "A. Poisson - linear",
      
      predicted_nb =
        "B. NB2 - linear",
      
      predicted_polynomial =
        "C. NB2 - polynomial",
      
      predicted_region_specific =
        "D. NB2 - nature-specific CLC"
    )
  )


ggplot(
  validation_long,
  aes(
    x = predicted_insct_r,
    y = insct_r,
    colour = nature_region
  )
) +
  
  geom_point(
    alpha = 0.5
  ) +
  
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed"
  ) +
  
  geom_smooth(
    method = "lm",
    se = TRUE
  ) +
  
  facet_wrap(
    ~ model,
    scales = "free"
  ) +
  
  labs(
    x = "Predicted insect richness",
    y = "Observed insect richness",
    colour = "Nature region",
    title = "Validation of insect richness calibration models"
  ) +
  
  theme_bw() +
  
  theme(
    legend.position = "bottom"
  )

################################################################################
# 29. Predicting observed richness
################################################################################

nature_regions_raster_crs <- nature_regions %>%
  sf::st_transform(terra::crs(allCalibrationRaster)) %>%
  mutate(
    nature_region = factor(
      nature_region,
      levels = levels(model_data$nature_region)
    ),
    nature_region_id = as.integer(nature_region)
  )

nature_region_raster <- terra::rasterize(
  terra::vect(nature_regions_raster_crs),
  allCalibrationRaster[[1]],
  field = "nature_region_id"
)

names(nature_region_raster) <- "nature_region"

################################################################################
# 7b. Categorical region raster
################################################################################

# nature_regions_pred$region_id <- as.integer(
#   nature_regions_pred$nature_region
# )
# 
# nature_region_raster <- terra::rasterize(
#   terra::vect(nature_regions_pred),
#   prediction_raster[[1]],
#   field = "region_id"
# )
# 
# names(nature_region_raster) <- "nature_region"


# pred_data <- terra::as.data.frame(
#   calibration_prediction_raster,
#   xy = TRUE,
#   na.rm = FALSE
# )
# 
# pred_data$nature_region_id <- terra::values(nature_region_raster)
# 
# pred_data$nature_region <- factor(
#   pred_data$nature_region_id,
#   levels = 1:5,
#   labels = c("Central", "East", "North", "South", "West")
# )
#   
################################################################################
# 8. Predict insect richness
################################################################################

prediction_raster <- c(calibration_prediction_raster,
                       nature_region_raster)



################################################################################
# 10. Prediction function
################################################################################

predict_insect_richness <- function(model, data) {
  
  # Convert raster region IDs (1-5) back to the factor levels
  # used when fitting the model
  region_levels <- levels(model_data$nature_region)
  
  data$nature_region <- factor(
    data$nature_region,
    levels = seq_along(region_levels),
    labels = region_levels
  )
  
  # Hold weight constant for spatial prediction
  data$weight <- 0
  
  predict(
    model,
    newdata = data,
    type = "response",
    allow.new.levels = FALSE
  )
}

################################################################################
# 11. Predict across the raster
################################################################################

predicted_insect_richness <- terra::predict(
  prediction_raster,
  insect_model,
  fun = predict_insect_richness,
  na.rm = TRUE
)

names(predicted_insect_richness) <- "predicted_insect_richness"




