################################################################################
# Fit insect calibration model and predict calibrated insect richness
################################################################################

fit_predict_insect_calibration <- function(
    observed_richness,
    allCalibrationRaster,
    plant_expected_richness,
    richness,
    nature_regions,
    calibration_vars,
    clc_vars,
    output_folder = NULL,
    cores = 1,
    overwrite = TRUE,
    seed = 12345,
    checks = FALSE,
    splitData = FALSE
) {
  
  ################################################################################
  # 1. Checks
  ################################################################################
  
  if (!inherits(observed_richness, "sf")) {
    stop("observed_richness must be an sf object.")
  }
  
  if (!inherits(allCalibrationRaster, "SpatRaster")) {
    stop("allCalibrationRaster must be a terra SpatRaster.")
  }
  
  if (!inherits(plant_expected_richness, "SpatRaster")) {
    stop("plant_expected_richness must be a terra SpatRaster.")
  }
  
  if (!inherits(richness, "SpatRaster")) {
    stop("richness must be a terra SpatRaster.")
  }
  
  if (!inherits(nature_regions, "sf")) {
    stop("nature_regions must be an sf object.")
  }
  
  if (!"insct_r" %in% names(observed_richness)) {
    stop("observed_richness must contain 'insct_r'.")
  }
  
  if (!"nature_region" %in% names(nature_regions)) {
    stop("nature_regions must contain 'nature_region'.")
  }
  
  
  ################################################################################
  # 2. Prepare plant expected richness
  ################################################################################
  
  # plant_expected_richness_prediction <- terra::project(
  #   plant_expected_richness,
  #   allCalibrationRaster[[1]],
  #   method = "bilinear"
  # )
  
  # plant_expected_richness_prediction <- plant_expected_richness %>%
  #   scale()
  # 
  # names(plant_expected_richness_prediction) <-
  #   "plant_expected_richness"
  
  
  ################################################################################
  # 3. Prepare richness predictor
  ################################################################################
  
  output_file <- file.path(
    output_folder,
    "richness_prediction_projected.tif"
  )
  
  if(!file.exists(output_file)){
    
    richness_prediction <- terra::project(
      richness,
      allCalibrationRaster[[1]],
      method = "bilinear",
      filename = output_file,
      overwrite = TRUE
    )
  } else {
    richness_prediction <- terra::rast(output_file)  
  }
  
  names(richness_prediction) <- "insect_richness"
  
  
  
  ################################################################################
  # 4. Create spatial predictor raster
  ################################################################################
  
  clc_prediction_raster <- allCalibrationRaster[[
    c(
      clc_vars,
      "landscape_heterogeneity_difference"
    )
  ]]
  
  
  calibration_prediction_raster <- c(
    clc_prediction_raster,
    #plant_expected_richness_prediction,
    richness_prediction
  )
  
  
  ################################################################################
  # 5. Check predictor names
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
  # 6. Extract predictors at monitoring locations
  ################################################################################
  
  extracted_predictors <- terra::extract(
    calibration_prediction_raster,
    terra::vect(observed_richness)
  ) %>%
    as.data.frame() %>%
    dplyr::select(-ID)
  
  
  ################################################################################
  # 7. Create calibration data
  ################################################################################
  
  insectCalibration <- observed_richness %>%
    sf::st_drop_geometry() %>%
    dplyr::select(
      insct_r,
      dplyr::any_of(
        c(
          "samp_dur",
          "year",
          "local"
        )
      )
    ) %>%
    dplyr::bind_cols(
      extracted_predictors
    )
  
  
  ################################################################################
  # 8. Add nature regions
  ################################################################################
  
  nature_region_values <- observed_richness %>%
    dplyr::select(geometry) %>%
    sf::st_join(
      nature_regions %>%
        dplyr::select(nature_region),
      join = sf::st_within,
      left = TRUE
    ) %>%
    sf::st_drop_geometry() %>%
    dplyr::select(nature_region)
  
  
  insectCalibration <- insectCalibration %>%
    dplyr::bind_cols(
      nature_region_values
    )
  
  
  ################################################################################
  # 9. Prepare modelling data
  ################################################################################
  
  model_data <- insectCalibration %>%
    dplyr::select(
      insct_r,
      nature_region,
      dplyr::all_of(calibration_vars),
      "samp_dur",
      "year",
      "local"
    ) %>%
    dplyr::mutate(
      nature_region = factor(nature_region)
    ) %>%
    dplyr::filter(
      complete.cases(.)
    ) 
  
  
  ################################################################################
  # 10. Create stratified train-validation split
  ################################################################################
  validation_dir <- file.path(output_folder,
                              "validation_output")
  if(!dir.exists(validation_dir)) dir.create(validation_dir)
  
  
  if(splitData){
    set.seed(seed)
    
    model_data <- model_data %>%
      dplyr::group_by(nature_region) %>%
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
      ) %>%
      dplyr::ungroup()
    
    
    train_data <- model_data %>%
      dplyr::filter(
        validation_set == "train"
      ) %>%
      dplyr::select(
        -validation_set
      )
    
    
    validation_data <- model_data %>%
      dplyr::filter(
        validation_set == "validation"
      ) %>%
      dplyr::select(
        -validation_set
      )
    
    if(checks){
      
      library(ggplot2)
      library(dplyr)
      library(tidyr)
      
      predictors <- train_data %>%
        select(where(is.numeric), -insct_r) %>%
        names()
      
      plot_data <- train_data %>%
        pivot_longer(
          cols = all_of(predictors),
          names_to = "predictor",
          values_to = "value"
        )
      
      ## --------------------------------------------------------
      ## Color-blind-friendly palette
      ## --------------------------------------------------------
      
      region_colours <- c(
        "Central" = "#0072B2",  # blue
        "East"    = "#D55E00",  # vermillion
        "North"   = "#CC79A7",  # reddish purple
        "South"   = "#009E73",  # bluish green
        "West"    = "#E69F00"   # orange
      )
      
      ## --------------------------------------------------------
      ## Plot
      ## --------------------------------------------------------
      
      p <- ggplot(
        plot_data,
        aes(
          x = value,
          y = insct_r,
          colour = nature_region
        )
      ) +
        geom_point(
          alpha = 0.5,
          size = 1.5
        ) +
        geom_smooth(
          method = "loess",
          se = TRUE
        ) +
        facet_wrap(
          ~ predictor,
          scales = "free_x"
        ) +
        coord_cartesian(
          ylim = c(0, 200)
        ) +
        scale_colour_manual(
          values = region_colours
        ) +
        theme_bw() +
        labs(
          x = "Predictor value",
          y = "Insect richness",
          colour = "Nature region",
          title = "Insect richness versus environmental predictors"
        ) +
        theme(
          legend.position = "bottom",
          strip.text = element_text(face = "bold")
        )
      
      print(p)
      
      ggsave(filename = file.path(validation_dir,
                                  "exploration.pdf"),
             plot = p,
             units = "in",
             width = 10,
             height = 10)
      
      ggsave(filename = file.path(validation_dir,
                                  "exploration.png"),
             plot = p,
             units = "in",
             width = 10,
             height = 10)
    }
    
    
    ################################################################################
    # 11. Scale calibration predictors using TRAINING data only
    ################################################################################
    
    ## ------------------------------------------------------------------------
    ## Variables to scale
    ##
    ## insect_richness is included because it is a continuous calibration
    ## predictor and should be treated consistently with the other predictors.
    ##
    ## Response:
    ##   insct_r          -> NOT scaled
    ##
    ## Offset:
    ##   weight           -> NOT scaled
    ##
    ## Random effect:
    ##   nature_region    -> NOT scaled
    ## ------------------------------------------------------------------------
    
    scale_vars <- calibration_vars
    
    
    ## ------------------------------------------------------------------------
    ## Function to calculate scaling parameters from training data
    ## ------------------------------------------------------------------------
    
    calculate_scaling_parameters <- function(
    data,
    variables
    ) {
      
      scaling_parameters <- purrr::map_dfr(
        variables,
        function(variable) {
          
          x <- data[[variable]]
          
          if (!is.numeric(x)) {
            
            stop(
              "Variable '",
              variable,
              "' is not numeric and cannot be scaled."
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
          
          if (!is.finite(sd_x) || sd_x == 0) {
            
            stop(
              "Variable '",
              variable,
              "' has zero or non-finite standard deviation ",
              "in the training data."
            )
          }
          
          tibble::tibble(
            variable = variable,
            mean = mean_x,
            sd = sd_x
          )
        }
      )
      
      scaling_parameters
    }
    
    
    ## ------------------------------------------------------------------------
    ## Function to apply scaling parameters
    ## ------------------------------------------------------------------------
    
    apply_scaling_parameters <- function(
    data,
    scaling_parameters
    ) {
      
      for (i in seq_len(nrow(scaling_parameters))) {
        
        variable <- scaling_parameters$variable[i]
        
        mean_x <- scaling_parameters$mean[i]
        
        sd_x <- scaling_parameters$sd[i]
        
        data[[variable]] <-
          (
            data[[variable]] - mean_x
          ) /
          sd_x
      }
      
      data
    }
    
    
    ## ------------------------------------------------------------------------
    ## Calculate scaling parameters from TRAINING data
    ## ------------------------------------------------------------------------
    
    scaling_parameters <- calculate_scaling_parameters(
      data = train_data,
      variables = scale_vars
    )
    
    
    ## ------------------------------------------------------------------------
    ## Apply training-set scaling to training data
    ## ------------------------------------------------------------------------
    
    train_data <- apply_scaling_parameters(
      data = train_data,
      scaling_parameters = scaling_parameters
    ) %>%
      dplyr::mutate(
        nature_region = factor(nature_region),
        year = factor(year)
      )
    
    
    ## ------------------------------------------------------------------------
    ## Apply EXACT SAME scaling to validation data
    ## ------------------------------------------------------------------------
    
    validation_data <- apply_scaling_parameters(
      data = validation_data,
      scaling_parameters = scaling_parameters
    ) %>%
      dplyr::mutate(
        nature_region = factor(nature_region),
        year = factor(year)
      )
    
    # Sampling duration must be positive
    stopifnot(all(train_data$samp_dur > 0, na.rm = TRUE))
    stopifnot(all(validation_data$samp_dur > 0, na.rm = TRUE))
    
    
    ## ------------------------------------------------------------------------
    ## Optional diagnostic
    ## ------------------------------------------------------------------------
    
    if (exists("checks") && checks) {
      
      cat("\n")
      cat("============================================================\n")
      cat("Scaling parameters from training data\n")
      cat("============================================================\n")
      
      print(scaling_parameters)
    }
    
    
    ################################################################################
    # 12. Model formulas
    ################################################################################
    
    #random_effects <- "+ (1 | nature_region) + (1 | year)"
    # offset_term <- "+ offset(log(samp_dur))"
    
    random_effects <- "+ nature_region + (1 | year)"
    
    ## ------------------------------------------------------------------------
    ## Linear formula
    ## ------------------------------------------------------------------------
    
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
    
    
    ## ------------------------------------------------------------------------
    ## Polynomial terms
    ## ------------------------------------------------------------------------
    
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
    
    
    ## ------------------------------------------------------------------------
    ## Polynomial formula
    ## ------------------------------------------------------------------------
    
    # polynomial_formula <- stats::as.formula(
    #   paste(
    #     "insct_r ~",
    #     paste(
    #       polynomial_terms,
    #       collapse = " + "
    #     ),
    #     "+ (1 | nature_region)"
    #   )
    # )
    
    
    ## ------------------------------------------------------------------------
    ## Nature-region-specific formula
    ## ------------------------------------------------------------------------
    
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
    
    
    ## ------------------------------------------------------------------------
    ## GAM formula
    ##
    ## k = 5 gives each continuous predictor a relatively simple smooth.
    ##
    ## nature_region is represented as a random-effect smooth.
    ## ------------------------------------------------------------------------
    
    gam_terms <- paste0(
      "s(",
      calibration_vars,
      ", k = 5)"
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
    # 13. Fit candidate models
    ################################################################################
    
    ################################################################################
    # 13. Fit candidate models
    ################################################################################
    
    models <- list(
      
      ## ----------------------------------------------------------------------
      ## A. Poisson - linear
      ## ----------------------------------------------------------------------
      
      "A. Poisson - linear" =
        
        glmmTMB::glmmTMB(
          formula = linear_formula,
          data = train_data,
          family = poisson(link = "log")
        ),
      
      
      ## ----------------------------------------------------------------------
      ## B. NB2 - linear
      ## ----------------------------------------------------------------------
      
      "B. NB2 - linear" =
        
        glmmTMB::glmmTMB(
          formula = linear_formula,
          data = train_data,
          family = glmmTMB::nbinom2(link = "log")
        ),
      
      
      ## ----------------------------------------------------------------------
      ## C. NB2 - polynomial
      ## ----------------------------------------------------------------------
      
      "C. NB2 - polynomial" =
        
        glmmTMB::glmmTMB(
          formula = polynomial_formula,
          data = train_data,
          family = glmmTMB::nbinom2(link = "log")
        ),
      
      
      ## ----------------------------------------------------------------------
      ## D. NB2 - nature-region-specific effects
      ## ----------------------------------------------------------------------
      
      "D. NB2 - nature-specific CLC" =
        
        glmmTMB::glmmTMB(
          formula = region_specific_formula,
          data = train_data,
          family = glmmTMB::nbinom2(link = "log")
        ),
      
      
      ## ----------------------------------------------------------------------
      ## E. GAM - negative binomial
      ## ----------------------------------------------------------------------
      
      "E. GAM - negative binomial" =
        
        mgcv::gam(
          formula = gam_formula,
          data = train_data,
          family = mgcv::nb(),
          method = "REML"
        ),
      
      ## ----------------------------------------------------------------------
      ## F. Poisson - polynomiaø
      ## ----------------------------------------------------------------------
      
      "F. Poisson - polynomial" =
        
        glmmTMB::glmmTMB(
          formula = polynomial_formula,
          data = train_data,
          family = poisson(link = "log")
        )
    )
    
    
    ################################################################################
    # 14. Validation predictions
    ################################################################################
    
    validation_predictions <- purrr::map(
      models,
      function(model) {
        
        stats::predict(
          model,
          newdata = validation_data,
          type = "response",
          allow.new.levels = FALSE
        )
      }
    )
    
    
    ## ------------------------------------------------------------------------
    ## Prediction names
    ## ------------------------------------------------------------------------
    
    prediction_names <- c(
      "predicted_poisson",
      "predicted_nb",
      "predicted_polynomial",
      "predicted_region_specific",
      "predicted_gam",
      "predicted_poisson_polynomial"
    )
    
    
    names(validation_predictions) <-
      prediction_names
    
    
    ## ------------------------------------------------------------------------
    ## Add predictions to validation data
    ## ------------------------------------------------------------------------
    
    validation_data <- dplyr::bind_cols(
      validation_data,
      as.data.frame(validation_predictions)
    )
    
    
    ################################################################################
    # 15. Validation statistics
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
    
    
    ## ------------------------------------------------------------------------
    ## Compare models
    ## ------------------------------------------------------------------------
    
    validation_comparison <- purrr::map2_dfr(
      
      prediction_names,
      
      names(models),
      
      ~ calculate_validation_stats(
        
        observed =
          validation_data$insct_r,
        
        predicted =
          validation_data[[.x]]
        
      ) %>%
        
        dplyr::mutate(
          model = .y,
          .before = 1
        )
      
    ) %>%
      
      dplyr::arrange(
        RMSE
      )
    
    
    ################################################################################
    # 16. Select best model
    ################################################################################
    
    best_model_name <- validation_comparison %>%
      
      dplyr::slice_min(
        RMSE,
        n = 1,
        with_ties = FALSE
      ) %>%
      
      dplyr::pull(
        model
      )
    
    
    insect_model <- models[[
      best_model_name
    ]]
    
    
    if (exists("checks") && checks) {
      
      cat("\n")
      cat("============================================================\n")
      cat("Validation comparison\n")
      cat("============================================================\n")
      
      print(validation_comparison)
      
      cat("\n")
      cat("Best model:\n")
      print(best_model_name)
      
      
      library(ggplot2)
      library(dplyr)
      library(tidyr)
      
      # -------------------------------------------------------------------------
      # Observed vs predicted data
      # -------------------------------------------------------------------------
      
      prediction_names <- c(
        "predicted_poisson",
        "predicted_nb",
        "predicted_polynomial",
        "predicted_region_specific",
        "predicted_gam",
        "predicted_poisson_polynomial"
      )
      
      plot_data <- validation_data %>%
        dplyr::select(
          insct_r,
          all_of(prediction_names)
        ) %>%
        tidyr::pivot_longer(
          cols = all_of(prediction_names),
          names_to = "model",
          values_to = "predicted"
        ) %>%
        dplyr::mutate(
          model = dplyr::recode(
            model,
            predicted_poisson = "Poisson - linear",
            predicted_nb = "NB2 - linear",
            predicted_polynomial = "NB2 - polynomial",
            predicted_region_specific = "NB2 - nature-specific CLC",
            predicted_gam = "GAM - negative binomial",
            predicted_poisson_polynomial = "Poisson - polynomial"
          )
        ) %>%
        dplyr::filter(
          is.finite(insct_r),
          is.finite(predicted)
        )
      
      
      # -------------------------------------------------------------------------
      # Observed vs predicted plot
      # -------------------------------------------------------------------------
      
      p <-  ggplot(
        plot_data,
        aes(
          y = insct_r,
          x = predicted
        )
      ) +
        
        geom_point(
          alpha = 0.65,
          size = 2
        ) +
        
        geom_abline(
          slope = 1,
          intercept = 0,
          linetype = "dashed",
          linewidth = 0.8
        ) +
        
        geom_smooth(
          method = "lm",
          se = FALSE,
          linewidth = 0.8
        ) +
        
        facet_wrap(
          ~ model,
          scales = "free"
        ) +
        
        labs(
          y = "Observed insect richness",
          x = "Predicted insect richness",
          title = "Observed versus predicted insect richness",
          subtitle = "Validation predictions from candidate models"
        ) +
        
        theme_bw() +
        
        theme(
          strip.text = element_text(
            face = "bold"
          ),
          plot.title = element_text(
            face = "bold"
          ),
          panel.grid.minor = element_blank()
        )
      
      print(p)
      
      ggsave(filename = file.path(validation_dir,
                                  "validation_plots.pdf"),
             plot = p,
             units = "in",
             width = 10,
             height = 10)
      
      ggsave(filename = file.path(validation_dir,
                                  "validation_plots.png"),
             plot = p,
             units = "in",
             width = 10,
             height = 10)  
      
    }
    
    
    
    
  } else {
    
    
    ################################################################################
    # Insect richness modelling
    #
    # Models:
    #
    # A. Poisson - linear + weight offset
    # B. NB1 - linear
    # C. NB1 - polynomial
    # D. NB1 - nature-specific CLC
    # E. Gaussian - log response
    # F. Gaussian - log response + polynomial
    # G. GAM - negative binomial
    #
    # Validation:
    #   - Stratified 5-fold cross-validation by nature_region
    #   - Out-of-fold predictions
    #   - RMSE
    #   - MAE
    #   - Bias
    #   - Pearson correlation
    #   - R2
    #
    ################################################################################
    
    
    # ==============================================================================
    # 10. Prepare data and stratified cross-validation
    # ==============================================================================
    
    set.seed(seed)
    
    k_folds <- 5
    
    model_data <- model_data %>%
      dplyr::filter(weight >= 0)
    
    ## ------------------------------------------------------------------------
    ## Check weight
    ## ------------------------------------------------------------------------
    
    if (!"weight" %in% names(model_data)) {
      stop(
        "'weight' is not present in model_data. ",
        "The Poisson offset model requires a 'weight' variable."
      )
    }
    
    if (anyNA(model_data$weight)) {
      stop(
        "'weight' contains missing values. ",
        "Remove or appropriately handle missing weights before model fitting."
      )
    }
    
    if (any(!is.finite(model_data$weight))) {
      stop(
        "'weight' contains non-finite values."
      )
    }
    
    if (any(model_data$weight <= 0)) {
      stop(
        "'weight' contains values <= 0. ",
        "The offset uses log(weight), so all weights must be strictly positive."
      )
    }
    
    
    ## ------------------------------------------------------------------------
    ## Check response
    ## ------------------------------------------------------------------------
    
    if (anyNA(model_data$insct_r)) {
      stop("'insct_r' contains missing values.")
    }
    
    if (any(model_data$insct_r < 0)) {
      stop("'insct_r' contains negative values.")
    }
    
    
    ## ------------------------------------------------------------------------
    ## Check nature_region
    ## ------------------------------------------------------------------------
    
    if (!"nature_region" %in% names(model_data)) {
      stop("'nature_region' is not present in model_data.")
    }
    
    model_data$nature_region <- as.factor(model_data$nature_region)
    
    
    ## ------------------------------------------------------------------------
    ## Check that all calibration variables exist
    ## ------------------------------------------------------------------------
    
    missing_calibration_vars <- setdiff(
      calibration_vars,
      names(model_data)
    )
    
    if (length(missing_calibration_vars) > 0) {
      stop(
        "The following calibration variables are missing from model_data: ",
        paste(missing_calibration_vars, collapse = ", ")
      )
    }
    
    
    ## ------------------------------------------------------------------------
    ## Inspect weight
    ## ------------------------------------------------------------------------
    
    if (checks) {
      
      cat("\n")
      cat("============================================================\n")
      cat("Weight diagnostics\n")
      cat("============================================================\n")
      
      print(summary(model_data$weight))
      
      cat("\nNumber of unique weights:",
          length(unique(model_data$weight)), "\n")
      
      cat("Minimum weight:",
          min(model_data$weight), "\n")
      
      cat("Maximum weight:",
          max(model_data$weight), "\n")
      
      cat("Geometric mean weight:",
          exp(mean(log(model_data$weight))), "\n")
    }
    
    
    ## ------------------------------------------------------------------------
    ## Create stratified folds
    ## ------------------------------------------------------------------------
    
    region_counts <- model_data %>%
      dplyr::count(nature_region)
    
    if (any(region_counts$n < k_folds)) {
      
      stop(
        "At least one nature_region has fewer observations than the ",
        k_folds,
        " required folds."
      )
    }
    
    
    model_data <- model_data %>%
      
      dplyr::group_by(nature_region) %>%
      
      dplyr::mutate(
        fold = sample(
          rep(
            seq_len(k_folds),
            length.out = dplyr::n()
          )
        )
      ) %>%
      
      dplyr::ungroup()
    
    
    ## ------------------------------------------------------------------------
    ## Check fold balance
    ## ------------------------------------------------------------------------
    
    fold_check <- model_data %>%
      dplyr::count(
        nature_region,
        fold
      )
    
    if (checks) {
      
      cat("\n")
      cat("============================================================\n")
      cat("Fold allocation\n")
      cat("============================================================\n")
      
      print(fold_check)
    }
    
    
    # ==============================================================================
    # 11. Model formulas
    # ==============================================================================
    
    
    ## ------------------------------------------------------------------------
    ## A. Poisson linear model with weight offset
    ## ------------------------------------------------------------------------
    
    poisson_offset_formula <- stats::as.formula(
      paste(
        "insct_r ~",
        paste(calibration_vars, collapse = " + "),
        "+ offset(log(weight))",
        "+ (1 | nature_region)"
      )
    )
    
    
    ## ------------------------------------------------------------------------
    ## Linear formula
    ## ------------------------------------------------------------------------
    
    linear_formula <- stats::as.formula(
      paste(
        "insct_r ~",
        paste(calibration_vars, collapse = " + "),
        "+ (1 | nature_region)"
      )
    )
    
    
    ## ------------------------------------------------------------------------
    ## Polynomial terms
    ## ------------------------------------------------------------------------
    
    polynomial_terms <- paste0(
      "poly(",
      calibration_vars,
      ", degree = 2, raw = TRUE)"
    )
    
    
    ## ------------------------------------------------------------------------
    ## Polynomial formula
    ## ------------------------------------------------------------------------
    
    polynomial_formula <- stats::as.formula(
      paste(
        "insct_r ~",
        paste(polynomial_terms, collapse = " + "),
        "+ (1 | nature_region)"
      )
    )
    
    
    ## ------------------------------------------------------------------------
    ## Nature-region-specific CLC effects
    ## ------------------------------------------------------------------------
    
    region_specific_formula <- stats::as.formula(
      paste(
        "insct_r ~",
        paste0(
          calibration_vars,
          " * nature_region",
          collapse = " + "
        )
      )
    )
    
    
    ## ------------------------------------------------------------------------
    ## Log-response formula
    ## ------------------------------------------------------------------------
    
    log_formula <- stats::as.formula(
      paste(
        "log(insct_r) ~",
        paste(calibration_vars, collapse = " + "),
        "+ (1 | nature_region)"
      )
    )
    
    
    ## ------------------------------------------------------------------------
    ## Log-response polynomial formula
    ## ------------------------------------------------------------------------
    
    log_polynomial_formula <- stats::as.formula(
      paste(
        "log(insct_r) ~",
        paste(polynomial_terms, collapse = " + "),
        "+ (1 | nature_region)"
      )
    )
    
    
    ## ------------------------------------------------------------------------
    ## GAM formula
    ## ------------------------------------------------------------------------
    
    gam_terms <- paste0(
      "s(",
      calibration_vars,
      ", k = 5)"
    )
    
    
    gam_formula <- stats::as.formula(
      paste(
        "insct_r ~",
        paste(gam_terms, collapse = " + "),
        "+ s(nature_region, bs = 're')"
      )
    )
    
    
    # ==============================================================================
    # 12. Candidate model fitting function
    # ==============================================================================
    
    fit_candidate_models <- function(data) {
      
      list(
        
        # --------------------------------------------------------------------------
        # A. Poisson - linear + weight offset
        # --------------------------------------------------------------------------
        
        "A. Poisson - linear + weight offset" =
          
          glmmTMB::glmmTMB(
            formula = poisson_offset_formula,
            data = data,
            family = poisson(link = "log")
          ),
        
        
        # --------------------------------------------------------------------------
        # B. Negative binomial 1 - linear
        # --------------------------------------------------------------------------
        
        "B. NB1 - linear" =
          
          glmmTMB::glmmTMB(
            formula = linear_formula,
            data = data,
            family = glmmTMB::nbinom1(link = "log")
          ),
        
        
        # --------------------------------------------------------------------------
        # C. Negative binomial 1 - polynomial
        # --------------------------------------------------------------------------
        
        "C. NB1 - polynomial" =
          
          glmmTMB::glmmTMB(
            formula = polynomial_formula,
            data = data,
            family = glmmTMB::nbinom1(link = "log")
          ),
        
        
        # --------------------------------------------------------------------------
        # D. Negative binomial 1 - nature-specific CLC
        # --------------------------------------------------------------------------
        
        "D. NB1 - nature-specific CLC" =
          
          glmmTMB::glmmTMB(
            formula = region_specific_formula,
            data = data,
            family = glmmTMB::nbinom1(link = "log")
          ),
        
        
        # --------------------------------------------------------------------------
        # E. Gaussian - log response
        # --------------------------------------------------------------------------
        
        "E. Gaussian - log response" =
          
          glmmTMB::glmmTMB(
            formula = log_formula,
            data = data,
            family = gaussian()
          ),
        
        
        # --------------------------------------------------------------------------
        # F. Gaussian - log response polynomial
        # --------------------------------------------------------------------------
        
        "F. Gaussian - log response polynomial" =
          
          glmmTMB::glmmTMB(
            formula = log_polynomial_formula,
            data = data,
            family = gaussian()
          ),
        
        
        # --------------------------------------------------------------------------
        # G. GAM - negative binomial
        # --------------------------------------------------------------------------
        
        "G. GAM - negative binomial" =
          
          mgcv::gam(
            formula = gam_formula,
            data = data,
            family = mgcv::nb(),
            method = "REML"
          )
      )
    }
    
    
    # ==============================================================================
    # 13. Prediction function
    # ==============================================================================
    
    predict_model <- function(
    model,
    newdata,
    model_name = NULL
    ) {
      
      ## ------------------------------------------------------------------------
      ## GAM
      ## ------------------------------------------------------------------------
      
      if (inherits(model, "gam")) {
        
        return(
          as.numeric(
            stats::predict(
              model,
              newdata = newdata,
              type = "response"
            )
          )
        )
      }
      
      
      ## ------------------------------------------------------------------------
      ## glmmTMB models
      ## ------------------------------------------------------------------------
      
      if (inherits(model, "glmmTMB")) {
        
        ## --------------------------------------------------------------
        ## Identify log-response models using model name
        ## --------------------------------------------------------------
        
        is_log_model <- !is.null(model_name) &&
          grepl(
            "^E\\.|^F\\.",
            model_name
          )
        
        
        ## --------------------------------------------------------------
        ## Log-response models
        ## --------------------------------------------------------------
        
        if (is_log_model) {
          
          predicted_log <- stats::predict(
            model,
            newdata = newdata,
            type = "response",
            allow.new.levels = FALSE
          )
          
          
          ## ------------------------------------------------------------
          ## Lognormal retransformation
          ##
          ## E[Y] = exp(mu + sigma^2 / 2)
          ## ------------------------------------------------------------
          
          sigma2 <- stats::sigma(model)^2
          
          predicted <- exp(
            predicted_log +
              sigma2 / 2
          )
          
          return(
            as.numeric(predicted)
          )
        }
        
        
        ## --------------------------------------------------------------
        ## Count models
        ##
        ## Includes:
        ##   A. Poisson + weight offset
        ##   B. NB1 linear
        ##   C. NB1 polynomial
        ##   D. NB1 nature-specific CLC
        ## --------------------------------------------------------------
        
        return(
          as.numeric(
            stats::predict(
              model,
              newdata = newdata,
              type = "response",
              allow.new.levels = FALSE
            )
          )
        )
      }
      
      
      ## ------------------------------------------------------------------------
      ## Unsupported model
      ## ------------------------------------------------------------------------
      
      stop(
        "Unsupported model class: ",
        paste(
          class(model),
          collapse = ", "
        )
      )
    }
    
    
    # ==============================================================================
    # 14. Validation statistics
    # ==============================================================================
    
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
      
      
      if (length(observed) < 2) {
        
        return(
          tibble::tibble(
            n = length(observed),
            RMSE = NA_real_,
            MAE = NA_real_,
            bias = NA_real_,
            correlation = NA_real_,
            R2 = NA_real_
          )
        )
      }
      
      
      correlation <- cor(
        observed,
        predicted,
        method = "pearson"
      )
      
      
      tibble::tibble(
        
        n = length(observed),
        
        RMSE =
          sqrt(
            mean(
              (observed - predicted)^2
            )
          ),
        
        MAE =
          mean(
            abs(
              observed - predicted
            )
          ),
        
        bias =
          mean(
            predicted - observed
          ),
        
        correlation =
          correlation,
        
        R2 =
          correlation^2
      )
    }
    
    
    # ==============================================================================
    # 15. Five-fold cross-validation
    # ==============================================================================
    
    if (checks) {
      
      cat("\n")
      cat("============================================================\n")
      cat("Starting 5-fold cross-validation\n")
      cat("============================================================\n")
    }
    
    
    ## ------------------------------------------------------------------------
    ## Storage objects
    ## ------------------------------------------------------------------------
    
    cv_predictions <- list()
    
    cv_statistics <- list()
    
    
    ## ------------------------------------------------------------------------
    ## Cross-validation loop
    ## ------------------------------------------------------------------------
    
    for (fold_i in seq_len(k_folds)) {
      
      if (checks) {
        
        cat("\n")
        cat("------------------------------------------------------------\n")
        cat("Fold:", fold_i, "/", k_folds, "\n")
        cat("------------------------------------------------------------\n")
      }
      
      
      ## ----------------------------------------------------------------------
      ## Training and validation data
      ## ----------------------------------------------------------------------
      
      train_data <- model_data %>%
        dplyr::filter(
          fold != fold_i
        )
      
      
      validation_data <- model_data %>%
        dplyr::filter(
          fold == fold_i
        )
      
      
      ## ----------------------------------------------------------------------
      ## Fit candidate models
      ## ----------------------------------------------------------------------
      
      fitted_models <- fit_candidate_models(
        train_data
      )
      
      
      ## ----------------------------------------------------------------------
      ## Predict validation data
      ## ----------------------------------------------------------------------
      
      fold_predictions <- lapply(
        names(fitted_models),
        function(model_name) {
          
          predict_model(
            model = fitted_models[[model_name]],
            newdata = validation_data,
            model_name = model_name
          )
        }
      )
      
      names(fold_predictions) <- names(fitted_models)
      
      
      ## ----------------------------------------------------------------------
      ## Store predictions
      ## ----------------------------------------------------------------------
      
      for (model_name in names(fold_predictions)) {
        
        cv_predictions[[length(cv_predictions) + 1]] <-
          validation_data %>%
          
          dplyr::select(
            insct_r,
            nature_region,
            fold
          ) %>%
          
          dplyr::mutate(
            
            model =
              model_name,
            
            observed =
              insct_r,
            
            predicted =
              as.numeric(
                fold_predictions[[model_name]]
              )
          )
      }
      
      
      ## ----------------------------------------------------------------------
      ## Calculate fold-specific validation statistics
      ## ----------------------------------------------------------------------
      
      for (model_name in names(fold_predictions)) {
        
        stats_i <- calculate_validation_stats(
          
          observed =
            validation_data$insct_r,
          
          predicted =
            fold_predictions[[model_name]]
        )
        
        
        cv_statistics[[length(cv_statistics) + 1]] <-
          stats_i %>%
          
          dplyr::mutate(
            
            model =
              model_name,
            
            fold =
              fold_i
          )
      }
      
      
      if (checks) {
        
        cat(
          "Completed fold",
          fold_i,
          "\n"
        )
      }
    }
    
    
    # ==============================================================================
    # 16. Combine cross-validation results
    # ==============================================================================
    
    cv_predictions_df <- dplyr::bind_rows(
      cv_predictions
    )
    
    
    cv_statistics_df <- dplyr::bind_rows(
      cv_statistics
    )
    
    
    # ==============================================================================
    # 17. Cross-validation model comparison
    # ==============================================================================
    
    cv_comparison <- cv_statistics_df %>%
      
      dplyr::group_by(model) %>%
      
      dplyr::summarise(
        
        mean_RMSE =
          mean(
            RMSE,
            na.rm = TRUE
          ),
        
        sd_RMSE =
          stats::sd(
            RMSE,
            na.rm = TRUE
          ),
        
        mean_MAE =
          mean(
            MAE,
            na.rm = TRUE
          ),
        
        sd_MAE =
          stats::sd(
            MAE,
            na.rm = TRUE
          ),
        
        mean_bias =
          mean(
            bias,
            na.rm = TRUE
          ),
        
        mean_correlation =
          mean(
            correlation,
            na.rm = TRUE
          ),
        
        mean_R2 =
          mean(
            R2,
            na.rm = TRUE
          ),
        
        .groups = "drop"
      ) %>%
      
      dplyr::arrange(
        mean_RMSE
      )
    
    
    if (checks) {
      
      cat("\n")
      cat("============================================================\n")
      cat("Cross-validation comparison\n")
      cat("============================================================\n")
      
      print(cv_comparison)
    }
    
    
    # ==============================================================================
    # 18. Out-of-fold prediction plot
    # ==============================================================================
    
    oof_plot <- cv_predictions_df %>%
      
      ggplot2::ggplot(
        ggplot2::aes(
          y = observed,
          x = predicted,
          colour = nature_region
        )
      ) +
      
      ggplot2::geom_point(
        alpha = 0.45,
        size = 1.5
      ) +
      
      ggplot2::geom_abline(
        intercept = 0,
        slope = 1,
        linetype = "dashed"
      ) +
      
      ggplot2::facet_wrap(
        ~ model,
        scales = "free"
      ) +
      
      ggplot2::labs(
        
        y = "Observed insect richness",
        
        x = "Predicted insect richness",
        
        colour = "Nature region",
        
        title =
          "Out-of-fold predictions from 5-fold cross-validation"
      ) +
      
      ggplot2::theme_bw()
    
    
    if (checks) {
      print(oof_plot)
    }
    
    
    # ==============================================================================
    # 19. Select best model
    # ==============================================================================
    
    best_model_name <- cv_comparison %>%
      
      dplyr::slice_min(
        mean_RMSE,
        n = 1,
        with_ties = FALSE
      ) %>%
      
      dplyr::pull(
        model
      )
    
    
    if (checks) {
      
      cat("\n")
      cat("============================================================\n")
      cat("Best model\n")
      cat("============================================================\n")
      
      cat(
        "Selected model:",
        best_model_name,
        "\n"
      )
    }
    
    
    # ==============================================================================
    # 20. Fit final model using all calibration data
    # ==============================================================================
    
    final_model_data <- model_data %>%
      
      dplyr::select(
        -fold
      )
    
    
    final_models <- fit_candidate_models(
      final_model_data
    )
    
    
    insect_model <- final_models[[
      best_model_name
    ]]
    
    
    # ==============================================================================
    # 21. Final model summary
    # ==============================================================================
    
    if (checks) {
      
      cat("\n")
      cat("============================================================\n")
      cat("Final selected model summary\n")
      cat("============================================================\n")
      
      print(
        summary(insect_model)
      )
    }
    
    
    # ==============================================================================
    # 22. Save cross-validation results
    # ==============================================================================
    
    if (exists("resultFolder")) {
      
      dir.create(
        resultFolder,
        recursive = TRUE,
        showWarnings = FALSE
      )
      
      
      utils::write.csv(
        cv_comparison,
        file = file.path(
          resultFolder,
          "insect_model_cv_comparison.csv"
        ),
        row.names = FALSE
      )
      
      
      utils::write.csv(
        cv_statistics_df,
        file = file.path(
          resultFolder,
          "insect_model_cv_statistics.csv"
        ),
        row.names = FALSE
      )
      
      
      utils::write.csv(
        cv_predictions_df,
        file = file.path(
          resultFolder,
          "insect_model_oof_predictions.csv"
        ),
        row.names = FALSE
      )
      
      
      ggplot2::ggsave(
        filename = file.path(
          resultFolder,
          "insect_model_oof_predictions.png"
        ),
        plot = oof_plot,
        width = 14,
        height = 10,
        dpi = 300
      )
      
      
      saveRDS(
        insect_model,
        file = file.path(
          resultFolder,
          "final_insect_model.rds"
        )
      )
    }
    
    
    # ==============================================================================
    # 23. Final diagnostic information
    # ==============================================================================
    
    if (checks) {
      
      cat("\n")
      cat("============================================================\n")
      cat("Model comparison complete\n")
      cat("============================================================\n")
      
      cat("\nBest model:\n")
      print(best_model_name)
      
      cat("\nCross-validation results:\n")
      print(cv_comparison)
      
      cat("\nWeight summary:\n")
      print(summary(model_data$weight))
    }
    
    
    ################################################################################
    # End of script
    ################################################################################
    
  }
  
  ################################################################################
  # 16. Create nature-region raster
  ################################################################################
  
  region_levels <- levels(
    model_data$nature_region
  )
  
  
  nature_regions_raster_crs <- nature_regions %>%
    
    sf::st_transform(
      terra::crs(
        allCalibrationRaster
      )
    ) %>%
    
    dplyr::mutate(
      
      nature_region = factor(
        nature_region,
        levels = region_levels
      ),
      
      nature_region_id =
        as.integer(
          nature_region
        )
    )
  
  
  if (
    any(
      is.na(
        nature_regions_raster_crs$
        nature_region_id
      )
    )
  ) {
    
    stop(
      "Nature-region values could not be ",
      "matched to model_data levels."
    )
  }
  
  
  nature_region_raster <- terra::rasterize(
    
    terra::vect(
      nature_regions_raster_crs
    ),
    
    allCalibrationRaster[[1]],
    
    field =
      "nature_region_id"
  )
  
  
  names(
    nature_region_raster
  ) <- "nature_region"
  
  
  region_levels_df <- data.frame(
    ID = 1:5,
    nature_region = region_levels
  )
  
  nature_region_raster <- terra::as.factor(
    nature_region_raster
  )
  
  levels(nature_region_raster) <- region_levels_df
  
  ################################################################################
  # 17. Combine predictors and region raster
  ################################################################################
  
  
  ################################################################################
  # Scale prediction raster using training-data scaling parameters
  ################################################################################
  
  scale_raster <- function(
    raster,
    scaling_parameters
  ) {
    
    # --------------------------------------------------------------------------
    # Variables that need scaling
    # --------------------------------------------------------------------------
    
    scaling_vars <- scaling_parameters$variable
    
    
    # --------------------------------------------------------------------------
    # Check that all required variables are present
    # --------------------------------------------------------------------------
    
    missing_vars <- setdiff(
      scaling_vars,
      names(raster)
    )
    
    if (length(missing_vars) > 0) {
      
      stop(
        "The following variables are missing from prediction_raster: ",
        paste(
          missing_vars,
          collapse = ", "
        )
      )
    }
    
    
    # --------------------------------------------------------------------------
    # Check that SDs are valid
    # --------------------------------------------------------------------------
    
    if (any(
      is.na(scaling_parameters$sd) |
      scaling_parameters$sd <= 0
    )) {
      
      stop(
        "Scaling parameters contain missing or non-positive SD values."
      )
    }
    
    
    # --------------------------------------------------------------------------
    # Apply training-data scaling
    #
    # z = (x - training_mean) / training_sd
    # --------------------------------------------------------------------------
    
    for (v in scaling_vars) {
      
      mean_v <- scaling_parameters$mean[
        scaling_parameters$variable == v
      ]
      
      sd_v <- scaling_parameters$sd[
        scaling_parameters$variable == v
      ]
      
      raster[[v]] <- (
        raster[[v]] - mean_v
      ) / sd_v
    }
    
    
    return(raster)
  }
  
  prediction_raster_scaled <- scale_raster(
    raster = calibration_prediction_raster,
    scaling_parameters = scaling_parameters
  )
  
  
  prediction_raster <- c(
    
    prediction_raster_scaled,
    
    nature_region_raster
    
  )
  
  
  
  
  ################################################################################
  # 18. Prediction function
  ################################################################################
  ################################################################################
  # Prediction function
  #
  # Uses:
  #   - fixed effects
  #   - nature_region random effect
  #
  # Excludes:
  #   - year random effect
  ################################################################################
  
  predict_insect_richness <- function(
    model,
    data
  ) {
    
    # --------------------------------------------------------------------------
    # glmmTMB
    # --------------------------------------------------------------------------
    
    if (inherits(model, "glmmTMB")) {
      
      # Ensure region is a factor with the same levels as during fitting
      data$nature_region <- factor(
        data$nature_region,
        levels = region_levels
      )
      
      
      # ------------------------------------------------------------------------
      # Fixed effects
      # ------------------------------------------------------------------------
      
      beta <- glmmTMB::fixef(model)$cond
      
      
      # Get fixed-effects-only formula
      fixed_formula <- stats::formula(
        model,
        fixed.only = TRUE
      )
      
      
      # Extract RHS
      rhs <- fixed_formula[[3]]
      
      
      # Create RHS-only formula
      rhs_formula <- stats::as.formula(
        paste(
          "~",
          paste(
            deparse(rhs),
            collapse = " "
          )
        )
      )
      
      
      # ------------------------------------------------------------------------
      # Model matrix
      # ------------------------------------------------------------------------
      
      X <- stats::model.matrix(
        rhs_formula,
        data = data
      )
      
      
      # Check that all fitted coefficients are present
      missing_columns <- setdiff(
        names(beta),
        colnames(X)
      )
      
      if (length(missing_columns) > 0) {
        
        stop(
          "Missing model-matrix columns: ",
          paste(
            missing_columns,
            collapse = ", "
          )
        )
      }
      
      
      # Keep exactly the columns used by the fitted model
      X <- X[
        ,
        names(beta),
        drop = FALSE
      ]
      
      
      # ------------------------------------------------------------------------
      # Fixed-effect linear predictor
      # ------------------------------------------------------------------------
      
      eta <- as.vector(
        X %*% beta
      )
      
      
      # ------------------------------------------------------------------------
      # Add nature-region random effect
      # ------------------------------------------------------------------------
      
      region_re <- glmmTMB::ranef(
        model
      )$cond$nature_region
      
      
      region_effect <- region_re[
        as.character(data$nature_region),
        "(Intercept)"
      ]
      
      
      # Check for unmatched regions
      if (anyNA(region_effect)) {
        
        stop(
          "Some nature_region values do not match the fitted model levels."
        )
      }
      
      
      eta <- eta + region_effect
      
      
      # ------------------------------------------------------------------------
      # Convert from log scale to response scale
      # ------------------------------------------------------------------------
      
      return(
        exp(eta)
      )
    }
    
    
    # --------------------------------------------------------------------------
    # GAM
    # --------------------------------------------------------------------------
    
    if (inherits(model, "gam")) {
      
      data$nature_region <- factor(
        data$nature_region,
        levels = region_levels
      )
      
      return(
        stats::predict(
          model,
          newdata = data,
          type = "response",
          exclude = "s(year)"
        )
      )
    }
    
    
    # --------------------------------------------------------------------------
    # Unsupported model
    # --------------------------------------------------------------------------
    
    stop(
      "Unsupported model class: ",
      paste(
        class(model),
        collapse = ", "
      )
    )
  }
  
  ################################################################################
  # 19. Predict calibrated insect richness
  ################################################################################
  
  if (
    is.null(output_folder)
  ) {
    
    predicted_insect_richness <-
      terra::predict(
        
        prediction_raster,
        
        insect_model,
        
        fun = predict_insect_richness,
        
        na.rm = TRUE,
        
        cores = cores,
        
        filename =
          file.path(
            output_folder,
            "predicted_insect_richness.tif"
          ),
        
        overwrite = overwrite
        
      )
    
  } else {
    
    if (
      !dir.exists(output_folder)
    ) {
      
      dir.create(
        output_folder,
        recursive = TRUE
      )
    }
    
    
    predicted_insect_richness <-
      terra::predict(
        
        prediction_raster,
        
        insect_model,
        
        fun = predict_insect_richness,
        
        na.rm = TRUE,
        
        cores = cores,
        
        filename =
          file.path(
            output_folder,
            "predicted_insect_richness.tif"
          ),
        
        overwrite = overwrite
        
      )
  }
  
  
  names(
    predicted_insect_richness
  ) <-
    "predicted_insect_richness"
  
  ################################################################################
  # 20. Return
  ################################################################################
  
  return(
    list(
      
      predicted_insect_richness =
        predicted_insect_richness,
      
      insect_model =
        insect_model,
      
      best_model_name =
        best_model_name,
      
      models =
        models,
      
      validation_comparison =
        validation_comparison,
      
      validation_data =
        validation_data,
      
      train_data =
        train_data,
      
      calibration_prediction_raster =
        calibration_prediction_raster,
      
      nature_region_raster =
        nature_region_raster
      
    )
  )
}

# 
# insect_calibration <- fit_predict_insect_calibration(
#   
#   observed_richness =
#     observed_richness,
#   
#   allCalibrationRaster =
#     allCalibrationRaster,
#   
#   plant_expected_richness =
#     plant_expected_richness,
#   
#   richness =
#     indicator_results$allRasts %>%
#     tidyterra::select(richness),
#   
#   nature_regions =
#     nature_regions,
#   
#   calibration_vars =
#     calibration_vars,
#   
#   clc_vars =
#     clc_vars,
#   
#   output_folder =
#     plotFolder,
#   
#   cores = 1
# )
