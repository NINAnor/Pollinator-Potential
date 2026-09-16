################################################################################
# Fit insect calibration models and predict spatial insect richness
#
# Workflow
# --------
# 1. Validate inputs
# 2. Project the insect-richness predictor
# 3. Create the nature-region raster
# 4. Build the calibration predictor raster
# 5. Extract predictors at insect monitoring locations
# 6. Build the modelling dataset
# 7. Create a stratified train/validation split
# 8. Scale predictors using TRAINING data only
# 9. Fit candidate models
# 10. Validate and select the best model
# 11. Scale the spatial prediction raster using training parameters
# 12. Predict spatial insect richness
# 13. Return model, predictions and diagnostics
#
# Important:
#   - insct_r is NOT scaled
#   - nature_region is NOT scaled
#   - year is NOT used as a spatial prediction predictor
#   - the year random effect is excluded from spatial predictions
################################################################################


fit_predict_insect_calibration <- function(
    observed_richness,
    allCalibrationRaster,
   # plant_expected_richness,
    richness,
    nature_regions,
    calibration_vars,
    clc_vars,
    study_region,
    output_folder = NULL,
   scaleCovs = FALSE,
   prediction_year = prediction_year,
    cores = 1,
    overwrite = TRUE,
    seed = 12345,
    checks = FALSE
) {
  
  
  ################################################################################
  # 1. INPUT CHECKS
  ################################################################################
  
  stopifnot(
    inherits(observed_richness, "sf"),
    inherits(allCalibrationRaster, "SpatRaster"),
    #inherits(plant_expected_richness, "SpatRaster"),
    inherits(richness, "SpatRaster"),
    inherits(nature_regions, "sf")
  )
  
  if (!"insct_r" %in% names(observed_richness)) {
    
    stop(
      "'observed_richness' must contain the response variable 'insct_r'."
    )
  }
  
  if (!"nature_region" %in% names(nature_regions)) {
    
    stop(
      "'nature_regions' must contain 'nature_region'."
    )
  }
  
  if (is.null(output_folder)) {
    
    validation_dir <- NULL
    
  } else {
    
    dir.create(
      output_folder,
      recursive = TRUE,
      showWarnings = FALSE
    )
    
    validation_dir <- file.path(
      output_folder,
      "validation_output"
    )
    
    dir.create(
      validation_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )
  }
  
  
  ################################################################################
  # 2. PROJECT INSECT-RICHNESS PREDICTOR
  #
  # The richness raster is projected to the same CRS, resolution and grid as
  # the calibration raster.
  #
  # bilinear interpolation is appropriate because insect_richness is continuous.
  ################################################################################
  
  richness_output <- if (!is.null(output_folder)) {
    file.path(
      output_folder,
      "richness_prediction_projected.tif"
    )
  } else {
    NULL
  }
  
  
  if (
    !is.null(richness_output) &&
    file.exists(richness_output) 
  ) {
    
    richness_prediction <- terra::rast(
      richness_output
    )
    
  } else {
    
    richness_prediction <- terra::project(
      richness,
      allCalibrationRaster[[1]],
      method = "bilinear",
      filename = richness_output,
      overwrite = overwrite
    )
  }
  
  
  names(
    richness_prediction
  ) <- "insect_richness"
  
  
  ################################################################################
  # 3. CREATE NATURE-REGION RASTER
  #
  # Region IDs are explicitly defined from the factor levels:
  #
  #   1 = Central
  #   2 = East
  #   3 = North
  #   4 = South
  #   5 = West
  #
  # This raster is aligned directly to the calibration raster.
  ################################################################################
  
  region_levels <- levels(
    factor(
      nature_regions$nature_region
    )
  )
  
  # Make sure the factor order is explicit.
  #
  # This prevents alphabetical/factor-order changes from silently changing
  # the region IDs.
  region_levels <- c(
    "Central",
    "East",
    "North",
    "South",
    "West"
  )
  
  nature_regions_raster_crs <- sf::st_transform(
    nature_regions,
    terra::crs(allCalibrationRaster)
  ) |>
    dplyr::mutate(
      nature_region = factor(
        nature_region,
        levels = region_levels
      ),
      nature_region_id = as.integer(
        nature_region
      )
    )
  
  
  if (
    anyNA(
      nature_regions_raster_crs$nature_region_id
    )
  ) {
    
    stop(
      "Some nature-region values could not be matched to region_levels."
    )
  }
  
  
  nature_region_raster <- terra::rasterize(
    terra::vect(
      nature_regions_raster_crs
    ),
    allCalibrationRaster[[1]],
    field = "nature_region_id"
  )
  
  names(
    nature_region_raster
  ) <- "nature_region"
  
  
  # Attach categorical labels for documentation and downstream use.
  
  region_levels_df <- data.frame(
    ID = seq_along(region_levels),
    nature_region = region_levels
  )
  
  nature_region_raster <- terra::as.factor(
    nature_region_raster
  )
  
  levels(
    nature_region_raster
  ) <- region_levels_df
  
  
  ################################################################################
  # 4. BUILD CALIBRATION PREDICTOR RASTER
  #
  # Only variables actually used by the candidate models are retained.
  ################################################################################
  
  clc_prediction_raster <- allCalibrationRaster[[
    c(
      clc_vars,
      "landscape_heterogeneity_difference"
    )
  ]]
  
  
  calibration_prediction_raster <- c(
    clc_prediction_raster,
    richness_prediction
  )
  
  
  names(
    calibration_prediction_raster
  ) <- c(
    clc_vars,
    "landscape_heterogeneity_difference",
    "insect_richness"
  )
  
  
  calibration_prediction_raster$insect_richness <- 
    scale(calibration_prediction_raster$insect_richness)
  
  ################################################################################
  # 5. CHECK PREDICTOR NAMES
  ################################################################################
  
  if (
    !identical(
      names(calibration_prediction_raster),
      calibration_vars
    )
  ) {
    
    stop(
      "Prediction raster layers do not match calibration_vars.\n",
      "Raster layers: ",
      paste(
        names(calibration_prediction_raster),
        collapse = ", "
      ),
      "\nModel predictors: ",
      paste(
        calibration_vars,
        collapse = ", "
      )
    )
  }
  
  
  ################################################################################
  # 6. EXTRACT SPATIAL PREDICTORS AT MONITORING LOCATIONS
  #
  # terra::extract() is considerably simpler than performing an sf spatial
  # join for every monitoring location.
  ################################################################################
  
  extracted_predictors <- terra::extract(
    calibration_prediction_raster,
    terra::vect(observed_richness)
  ) |>
    dplyr::select(
      -ID
    )
  
  
  # Extract nature-region ID from the same aligned raster.
  
  extracted_region <- terra::extract(
    nature_region_raster,
    terra::vect(observed_richness)
  ) |>
    dplyr::pull(
      nature_region
    )
  
  
  ################################################################################
  # 7. BUILD CALIBRATION DATASET
  ################################################################################
  
  insectCalibration <- observed_richness |>
    sf::st_drop_geometry() |>
    dplyr::select(
      insct_r,
      dplyr::any_of(
        c(
          "samp_dur",
          "year",
          "local"
        )
      )
    ) |>
    dplyr::bind_cols(
      extracted_predictors
    ) |>
    dplyr::mutate(
      nature_region = factor(
        extracted_region#,
        # levels = seq_along(region_levels),
        # labels = region_levels
      )
    )
  
  
  ################################################################################
  # 8. PREPARE MODELLING DATA
  ################################################################################
  
  model_data <- insectCalibration |>
    dplyr::select(
      insct_r,
      nature_region,
      dplyr::all_of(calibration_vars),
      "samp_dur",
      "year",
      "local"
    ) |>
    dplyr::mutate(
      nature_region = factor(
        nature_region,
        levels = region_levels
      ),
      year = factor(
        year
      )
    ) 
  
  
  ################################################################################
  # 9. TRAIN / VALIDATION SPLIT
  #
  # Stratification is performed by nature region so that every region is
  # represented in both training and validation datasets.
  ################################################################################
  
  set.seed(
    seed
  )
  
  model_data <- model_data |>
    dplyr::group_by(
      nature_region
    ) |>
    dplyr::mutate(
      validation_set = sample(
        c(
          rep(
            "train",
            floor(0.80 * dplyr::n())
          ),
          rep(
            "validation",
            dplyr::n() -
              floor(0.80 * dplyr::n())
          )
        )
      )
    ) |>
    dplyr::ungroup()
  
  
  train_data <- model_data |>
    dplyr::filter(
      validation_set == "train"
    ) |>
    dplyr::select(
      -validation_set
    )
  
  
  validation_data <- model_data |>
    dplyr::filter(
      validation_set == "validation"
    ) |>
    dplyr::select(
      -validation_set
    )
  
  
  ################################################################################
  # 10. SCALE CONTINUOUS PREDICTORS
  #
  # IMPORTANT:
  #
  # Scaling parameters are calculated ONLY from the training data.
  #
  # This prevents information from the validation dataset from entering
  # the model-fitting process.
  ################################################################################
  
  if(scaleCovs){
  scale_vars <- calibration_vars
  
  
  calculate_scaling_parameters <- function(
    data,
    variables
  ) {
    
    purrr::map_dfr(
      variables,
      function(variable) {
        
        x <- data[[variable]]
        
        if (!is.numeric(x)) {
          
          stop(
            "Variable '",
            variable,
            "' must be numeric."
          )
        }
        
        mean_x <- mean(
          x,
          na.rm = TRUE
        )
        
        sd_x <- stats::sd(
          x,
          na.rm = TRUE
        )
        
        if (
          !is.finite(sd_x) ||
          sd_x <= 0
        ) {
          
          stop(
            "Variable '",
            variable,
            "' has an invalid standard deviation."
          )
        }
        
        tibble::tibble(
          variable = variable,
          mean = mean_x,
          sd = sd_x
        )
      }
    )
  }
  
  
  apply_scaling_parameters <- function(
    data,
    scaling_parameters
  ) {
    
    for (i in seq_len(
      nrow(scaling_parameters)
    )) {
      
      variable <- scaling_parameters$variable[i]
      
      data[[variable]] <- (
        data[[variable]] -
          scaling_parameters$mean[i]
      ) /
        scaling_parameters$sd[i]
    }
    
    data
  }
  
  
  # Calculate parameters from TRAINING data only.
  
  scaling_parameters <- calculate_scaling_parameters(
    train_data,
    scale_vars
  )
  
  
  # Scale training data.
  
  train_data <- apply_scaling_parameters(
    train_data,
    scaling_parameters
  )
  
  
  # Scale validation data using EXACTLY the same parameters.
  
  validation_data <- apply_scaling_parameters(
    validation_data,
    scaling_parameters
  )
  }
  
  ################################################################################
  # 11. CHECK MODEL DATA
  ################################################################################
  
  stopifnot(
    all(
      is.finite(
        as.matrix(
          train_data[, calibration_vars]
        )
      )
    )
  )
  
  
  ################################################################################
  # 12. OPTIONAL EXPLORATORY DIAGNOSTICS
  ################################################################################
  
  if (checks) {
    
    predictors <- train_data |>
      dplyr::select(
        dplyr::where(is.numeric),
        -insct_r
      ) |>
      names()
    
    plot_data <- train_data |>
      tidyr::pivot_longer(
        cols = dplyr::all_of(predictors),
        names_to = "predictor",
        values_to = "value"
      )
    
    region_colours <- c(
      "Central" = "#0072B2",
      "East"    = "#D55E00",
      "North"   = "#CC79A7",
      "South"   = "#009E73",
      "West"    = "#E69F00"
    )
    
    p <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(
        x = value,
        y = insct_r,
        colour = nature_region
      )
    ) +
      ggplot2::geom_point(
        alpha = 0.5,
        size = 1.5
      ) +
      ggplot2::geom_smooth(
        method = "loess",
        se = TRUE
      ) +
      ggplot2::facet_wrap(
        ~ predictor,
        scales = "free_x"
      ) +
      ggplot2::coord_cartesian(
        ylim = c(0, 200)
      ) +
      ggplot2::scale_colour_manual(
        values = region_colours
      ) +
      ggplot2::theme_bw() +
      ggplot2::labs(
        x = "Scaled predictor",
        y = "Insect richness",
        colour = "Nature region",
        title = "Insect richness versus environmental predictors"
      ) +
      ggplot2::theme(
        legend.position = "bottom",
        strip.text = ggplot2::element_text(
          face = "bold"
        )
      )
    
    print(p)
    
    if (!is.null(validation_dir)) {
      
      ggplot2::ggsave(
        file.path(
          validation_dir,
          "exploration.pdf"
        ),
        p,
        units = "in",
        width = 10,
        height = 10
      )
      
      ggplot2::ggsave(
        file.path(
          validation_dir,
          "exploration.png"
        ),
        p,
        units = "in",
        width = 10,
        height = 10,
        dpi = 300
      )
    }
  }
  
  
  ################################################################################
  # 13. MODEL FORMULAS
  ################################################################################
  
  # Nature region is a FIXED EFFECT.
  #
  # Year is a RANDOM EFFECT.
  #
  # Sampling duration is NOT included because it was not part of the
  # model specification you selected.
  
  random_effects <- "+ nature_region + (1 | year)"
  
  
  linear_formula <- stats::as.formula(
    paste(
      "insct_r ~",
      paste(
        calibration_vars,
        collapse = " + "
      ),
      random_effects
    )
  )
  
  
  polynomial_terms <- paste0(
    "poly(",
    calibration_vars,
    ", degree = 2, raw = TRUE)"
  )
  
  
  polynomial_formula <- stats::as.formula(
    paste(
      "insct_r ~",
      paste(
        polynomial_terms,
        collapse = " + "
      ),
      random_effects
    )
  )
  
  
  region_specific_formula <- stats::as.formula(
    paste(
      "insct_r ~",
      paste(
        paste0(
          calibration_vars,
          " * nature_region"
        ),
        collapse = " + "
      ),
      "+ (1 | year)"
    )
  )
  
  
  gam_terms <- paste0(
    "s(",
    calibration_vars,
    ", k = 3)"
  )
  
  
  gam_formula <- stats::as.formula(
    paste(
      "insct_r ~",
      paste(
        gam_terms,
        collapse = " + "
      ),
      "+ s(nature_region, bs = 're')",
      "+ s(year, bs = 're')"
    )
  )
  
  
  ################################################################################
  # 14. FIT CANDIDATE MODELS
  ################################################################################
  
  models <- list(
    
    "A. Poisson - linear" =
      glmmTMB::glmmTMB(
        linear_formula,
        data = train_data,
        family = poisson(
          link = "log"
        )
      ),
    
    "B. NB2 - linear" =
      glmmTMB::glmmTMB(
        linear_formula,
        data = train_data,
        family = glmmTMB::nbinom2(
          link = "log"
        )
      ),
    
    # "C. NB2 - polynomial" =
    #   glmmTMB::glmmTMB(
    #     polynomial_formula,
    #     data = train_data,
    #     family = glmmTMB::nbinom2(
    #       link = "log"
    #     )
    #   ),
    # 
    # "D. NB2 - nature-specific CLC" =
    #   glmmTMB::glmmTMB(
    #     region_specific_formula,
    #     data = train_data,
    #     family = glmmTMB::nbinom2(
    #       link = "log"
    #     )
    #   ),
    
    "C. GAM - negative binomial" =
      mgcv::gam(
        gam_formula,
        data = train_data,
        family = mgcv::nb(),
        method = "REML"
      )#,
    
    # "F. Poisson - polynomial" =
    #   glmmTMB::glmmTMB(
    #     polynomial_formula,
    #     data = train_data,
    #     family = poisson(
    #       link = "log"
    #     )
    #   )
  )
  
  
  ################################################################################
  # 15. VALIDATION PREDICTIONS
  ################################################################################
  
  # For the glmmTMB models, validation predictions currently include the
  # estimated year random effect because year is known in validation_data.
  #
  # If your goal is to evaluate GENERALISATION to new years, change this to
  # re.form = NA.
  #
  # For consistency with your final spatial prediction, I recommend excluding
  # the year random effect here as well.
  
  validation_predictions <- purrr::map(
    models,
    function(model) {
      
      if (inherits(model, "glmmTMB")) {
        
        stats::predict(
          model,
          newdata = validation_data,
          type = "response",
          #re.form = NA,
          allow.new.levels = FALSE
        )
        
      } else {
        
        stats::predict(
          model,
          newdata = validation_data,
          type = "response"
        )
      }
    }
  )
  
  
  prediction_names <- c(
    "predicted_poisson",
    "predicted_nb",
    # "predicted_polynomial",
    # "predicted_region_specific",
    "predicted_gam"#,
    #"predicted_poisson_polynomial"
  )
  
  names(
    validation_predictions
  ) <- prediction_names
  
  
  validation_data <- dplyr::bind_cols(
    validation_data,
    as.data.frame(
      validation_predictions
    )
  )
  
  
  ################################################################################
  # 16. VALIDATION STATISTICS
  ################################################################################
  
  calculate_validation_stats <- function(
    observed,
    predicted
  ) {
    
    keep <- is.finite(observed) &
      is.finite(predicted)
    
    observed <- observed[keep]
    
    predicted <- predicted[keep]
    
    correlation <- cor(
      observed,
      predicted,
      method = "pearson"
    )
    
    tibble::tibble(
      
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
  
  
  validation_comparison <- purrr::map2_dfr(
    
    prediction_names,
    
    names(models),
    
    ~ calculate_validation_stats(
      observed = validation_data$insct_r,
      predicted = validation_data[[.x]]
    ) |>
      dplyr::mutate(
        model = .y,
        .before = 1
      )
    
  ) |>
    dplyr::arrange(
      RMSE
    )
  
  
  ################################################################################
  # 17. SELECT BEST MODEL
  ################################################################################
  
  best_model_name <- validation_comparison |>
    dplyr::slice_min(
      RMSE,
      n = 1,
      with_ties = FALSE
    ) |>
    dplyr::pull(
      model
    )
  
  
  insect_model <- models[[
    best_model_name
  ]]
  
  
  ################################################################################
  # 18. OPTIONAL VALIDATION PLOT
  ################################################################################
  
  if (checks) {
    
    plot_data <- validation_data |>
      dplyr::select(
        insct_r,
        dplyr::all_of(prediction_names)
      ) |>
      tidyr::pivot_longer(
        cols = dplyr::all_of(prediction_names),
        names_to = "model",
        values_to = "predicted"
      ) |>
      dplyr::mutate(
        model = dplyr::recode(
          model,
          predicted_poisson =
            "Poisson - linear",
          predicted_nb =
            "NB2 - linear",
          # predicted_polynomial =
          #   "NB2 - polynomial",
          # predicted_region_specific =
          #   "NB2 - nature-specific CLC",
          predicted_gam =
            "GAM - negative binomial"#,
          # predicted_poisson_polynomial =
          #   "Poisson - polynomial"
        )
      ) |>
      dplyr::filter(
        is.finite(insct_r),
        is.finite(predicted)
      )
    
    
    p <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(
        x = predicted,
        y = insct_r
      )
    ) +
      
      ggplot2::geom_point(
        alpha = 0.65,
        size = 2
      ) +
      
      ggplot2::geom_abline(
        slope = 1,
        intercept = 0,
        linetype = "dashed",
        linewidth = 0.8
      ) +
      
      ggplot2::geom_smooth(
        method = "lm",
        se = FALSE,
        linewidth = 0.8
      ) +
      
      ggplot2::facet_wrap(
        ~ model,
        scales = "free"
      ) +
      
      ggplot2::labs(
        y = "Observed insect richness",
        x = "Predicted insect richness",
        title = "Observed versus predicted insect richness",
        subtitle = paste(
          "Validation predictions; year random effect included"
        )
      ) +
      
      ggplot2::theme_bw() +
      
      ggplot2::theme(
        strip.text = ggplot2::element_text(
          face = "bold"
        ),
        plot.title = ggplot2::element_text(
          face = "bold"
        ),
        panel.grid.minor = ggplot2::element_blank()
      )
    
    
    print(p)
    
    
    if (!is.null(validation_dir)) {
      
      ggplot2::ggsave(
        file.path(
          validation_dir,
          "validation_plots.pdf"
        ),
        p,
        units = "in",
        width = 10,
        height = 10
      )
      
      ggplot2::ggsave(
        file.path(
          validation_dir,
          "validation_plots.png"
        ),
        p,
        units = "in",
        width = 10,
        height = 10,
        dpi = 300
      )
    }
  }
  
  
  ################################################################################
  # 19. SCALE SPATIAL PREDICTORS
  #
  # The raster is transformed using the means and SDs calculated from the
  # TRAINING DATA.
  ################################################################################
  
  prediction_raster_scaled <- calibration_prediction_raster
  
  # for (i in seq_len(
  #   nrow(scaling_parameters)
  # )) {
  #   
  #   variable <- scaling_parameters$variable[i]
  #   
  #   prediction_raster_scaled[[variable]] <-
  #     (
  #       prediction_raster_scaled[[variable]] -
  #         scaling_parameters$mean[i]
  #     ) /
  #     scaling_parameters$sd[i]
  # }
  
  
  ################################################################################
  # 20. ADD NATURE REGION TO PREDICTION RASTER
  ################################################################################
  
  prediction_raster <- c(
    prediction_raster_scaled,
    nature_region_raster
  ) %>%
    terra::mask(
    .,
    study_region
  )
  
  ################################################################################
  # 21. FINAL PREDICTION FUNCTION
  #
  # The selected model contains:
  #
  #   fixed effects:
  #       calibration predictors
  #       nature_region
  #
  #   random effect:
  #       year
  #
  # re.form = NA removes the year random effect.
  #
  # nature_region remains because it is a FIXED EFFECT.
  ################################################################################
  
  predict_insect_richness <- function(
    model,
    data
  ) {
    
    # terra may pass the region raster as its numeric IDs.
    #
    # Convert:
    #   1 -> Central
    #   2 -> East
    #   3 -> North
    #   4 -> South
    #   5 -> West
    
    data$nature_region <- factor(
      data$nature_region
    )
    
    data$year <- prediction_year
    
    data$year <- factor(
      data$year,
      levels = levels(train_data$year)
    )
    
    
    # --------------------------------------------------------------
    # GAM
    # --------------------------------------------------------------
    
    if (inherits(model, "gam")) {
      
      return(
        mgcv::predict.gam(
          model,
          newdata = data,
          type = "response"#,
          #exclude = "s(year)"
        )
      )
    }
    
    
    # --------------------------------------------------------------
    # glmmTMB
    # --------------------------------------------------------------
    
    if (inherits(model, "glmmTMB")) {
      
      return(
        stats::predict(
          model,
          newdata = data,
          type = "response",
          #re.form = NA,
          allow.new.levels = FALSE
        )
      )
    }
  }
  
  
  ################################################################################
  # 22. TEST PREDICTION BEFORE PROCESSING THE FULL RASTER
  #
  # This is important because the raster contains hundreds of millions of cells.
  ################################################################################
  
  # test_data <- terra::as.data.frame(
  #   prediction_raster,
  #   na.rm = TRUE
  # )
  # 
  # if (nrow(test_data) == 0) {
  #   
  #   stop(
  #     "No valid cells are available in the prediction raster."
  #   )
  # }
  # 
  # 
  # set.seed(
  #   seed
  # )
  # 
  # test_index <- sample(
  #   seq_len(nrow(test_data)),
  #   min(
  #     1000,
  #     nrow(test_data)
  #   )
  # )
  # 
  # test_prediction <- predict_insect_richness(
  #   insect_model,
  #   test_data[test_index, , drop = FALSE]
  # )
  # 
  # 
  # if (
  #   any(
  #     !is.finite(test_prediction)
  #   )
  # ) {
  #   
  #   stop(
  #     "Spatial prediction test produced non-finite values."
  #   )
  # }
  # 
  # 
  # if (checks) {
  #   
  #   cat("\n")
  #   cat("============================================================\n")
  #   cat("Spatial prediction test\n")
  #   cat("============================================================\n")
  #   
  #   print(
  #     summary(test_prediction)
  #   )
  #   
  #   cat("\nRange:\n")
  #   
  #   print(
  #     range(
  #       test_prediction,
  #       na.rm = TRUE
  #     )
  #   )
  # }
  
  
  ################################################################################
  # 23. PREDICT SPATIAL INSECT RICHNESS
  ################################################################################
  
  prediction_output <- if (!is.null(output_folder)) {
    
    file.path(
      output_folder,
      paste0("predicted_insect_richness_", prediction_year,".tif")
    )
    
  } else {
    
    NULL
  }
  
  
  predicted_insect_richness <- terra::predict(
    prediction_raster,
    insect_model,
    fun = predict_insect_richness,
    na.rm = TRUE,
    cores = cores,
    filename = prediction_output,
    overwrite = overwrite
  )
  
  
  names(
    predicted_insect_richness
  ) <- "predicted_insect_richness"
  
  
  ################################################################################
  # 24. RETURN RESULTS
  ################################################################################
  
  return(
    list(
      
      # Final spatial prediction
      predicted_insect_richness =
        predicted_insect_richness,
      
      # Selected model
      insect_model =
        insect_model,
      
      # Name of selected model
      best_model_name =
        best_model_name,
      
      # All candidate models
      models =
        models,
      
      # Validation performance
      validation_comparison =
        validation_comparison,
      
      # Validation observations and predictions
      validation_data =
        validation_data,
      
      # Training data after scaling
      train_data =
        train_data,
      
      # Scaling parameters used for the model
      # scaling_parameters =
      #   scaling_parameters,
      
      # Unscaled predictor raster
      calibration_prediction_raster =
        calibration_prediction_raster,
      
      # Scaled prediction raster
      prediction_raster_scaled =
        prediction_raster_scaled,
      
      # Nature-region raster
      nature_region_raster =
        nature_region_raster,
      
      # Region levels
      region_levels =
        region_levels
    )
  )
}