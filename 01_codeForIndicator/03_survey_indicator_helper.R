################################################################################
# Shared helper for survey-location indicator estimates (ASO, ANO, ASO + ANO)
#
# Uses the SAME raster and layer as 03b_get_indices_plot.R:
#   <plotFolder>/calibrated_insect_richness_normalised.tif
#
# For a set of survey points the function:
#   1. assigns each point to a Nature region,
#   2. extracts the indicator value at each point,
#   3. calculates descriptive statistics from those point values,
#   4. bootstraps the beta-regression mean (rows resampled with replacement,
#      bootstrap size = number of valid locations). The national model has a
#      random intercept for Nature region, as in 03b; regional models are
#      intercept-only.
#
# Output columns match 03b so 03e_processIndicesResults.R can bind them.
################################################################################

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(dplyr)
  library(tibble)
  library(betareg)
  library(glmmTMB)
})

survey_indicator_settings <- list(
  raster_file = file.path(plotFolder, "calibrated_insect_richness_normalised.tif"),
  indicator_name = "Calibrated regional indicator",
  B = 500L,
  threshold = 0.6,
  seed = 1000L
)

.adjust_beta <- function(y) {
  n <- length(y)
  if (any(y <= 0 | y >= 1)) y <- (y * (n - 1) + 0.5) / n
  y
}

.estimate_beta_mean <- function(y, region = NULL, national = FALSE) {
  if (length(y) < 2L) return(NA_real_)
  if (all(y == 1)) return(1)
  if (all(y == 0)) return(0)
  y <- .adjust_beta(y)

  if (national) {
    dat <- data.frame(y = y, nature_region = droplevels(factor(region)))
    fit <- tryCatch(
      glmmTMB::glmmTMB(
        y ~ 1 + (1 | nature_region), data = dat,
        family = glmmTMB::beta_family(link = "logit")
      ),
      error = function(e) NULL
    )
    if (is.null(fit)) return(mean(y))
    return(stats::plogis(glmmTMB::fixef(fit)$cond[["(Intercept)"]]))
  }

  fit <- tryCatch(
    betareg::betareg(y ~ 1, data = data.frame(y = y), link = "logit"),
    error = function(e) NULL
  )
  if (is.null(fit)) return(mean(y))
  unname(stats::plogis(stats::coef(fit, model = "mean")[["(Intercept)"]]))
}

.estimate_group <- function(dat, label, national, B, threshold, seed, indicator) {
  n <- nrow(dat)
  message("\nProcessing: ", label, " (", n, " locations)")

  est <- vapply(seq_len(B), function(b) {
    if (b == 1L || b %% 50L == 0L || b == B) message("  Bootstrap ", b, "/", B)
    set.seed(seed + b)
    idx <- sample.int(n, size = n, replace = TRUE)
    .estimate_beta_mean(dat$y[idx], dat$nature_region[idx], national)
  }, numeric(1))

  tibble::tibble(
    nature_region = label,
    indicator = indicator,
    global_mean = mean(dat$y),
    global_sd = stats::sd(dat$y),
    proportion_ge_0.6 = mean(dat$y >= threshold),
    percentage_ge_0.6 = mean(dat$y >= threshold) * 100,
    n_valid_cells = NA_real_,
    n_sampled_cells = n,
    bootstrap_mean = mean(est, na.rm = TRUE),
    bootstrap_sd = stats::sd(est, na.rm = TRUE),
    lower_quantile = unname(stats::quantile(est, 0.025, na.rm = TRUE)),
    upper_quantile = unname(stats::quantile(est, 0.975, na.rm = TRUE)),
    n_successful_bootstraps = sum(is.finite(est))
  )
}

run_survey_indicator <- function(locations_sf, out_prefix, force = FALSE,
                                 settings = survey_indicator_settings) {
  output_rds <- file.path(plotFolder, paste0(out_prefix, "_region_results.rds"))
  output_gpkg <- file.path(plotFolder, paste0(out_prefix, "_region_map.gpkg"))

  if (file.exists(output_rds) && file.exists(output_gpkg) && !force) {
    message("Existing results found; skipping: ", output_rds)
    return(invisible(readRDS(output_rds)))
  }

  stopifnot(file.exists(settings$raster_file))
  stopifnot(inherits(locations_sf, "sf"), "nature_region" %in% names(nature_regions))

  r <- terra::rast(settings$raster_file)
  if (terra::nlyr(r) != 1L) stop("Expected a single-layer indicator raster.")
  names(r) <- settings$indicator_name

  regions_sf <- sf::st_as_sf(nature_regions) |>
    sf::st_transform(terra::crs(r)) |>
    dplyr::select(nature_region)

  pts <- sf::st_transform(sf::st_geometry(locations_sf), terra::crs(r)) |>
    sf::st_as_sf()

  # Points on a border are assigned to the first matching region
  hits <- sf::st_intersects(pts, regions_sf)
  region_index <- vapply(
    hits, function(x) if (length(x)) x[[1L]] else NA_integer_, integer(1)
  )

  vals <- terra::extract(r, terra::vect(pts))[[2]]

  dat <- tibble::tibble(
    nature_region = as.character(regions_sf$nature_region[region_index]),
    y = vals
  )

  message(
    out_prefix, ": ", nrow(dat), " locations; ",
    sum(is.na(dat$nature_region)), " outside all regions; ",
    sum(is.na(dat$y)), " with NA raster value."
  )
  dat <- dat[!is.na(dat$nature_region) & is.finite(dat$y), ]
  if (nrow(dat) < 2L) stop("Fewer than two valid locations for ", out_prefix)
  if (any(dat$y < 0 | dat$y > 1)) {
    stop("Indicator values outside [0, 1] at survey locations.")
  }

  args <- list(
    B = settings$B, threshold = settings$threshold,
    indicator = settings$indicator_name
  )

  national <- do.call(.estimate_group, c(
    list(dat = dat, label = "National", national = TRUE,
         seed = settings$seed), args
  ))

  region_results <- dplyr::bind_rows(lapply(
    sort(unique(dat$nature_region)), function(reg) {
      do.call(.estimate_group, c(
        list(dat = dat[dat$nature_region == reg, ], label = reg,
             national = FALSE,
             seed = settings$seed + match(reg, regions_sf$nature_region) * 10000L),
        args
      ))
    }
  ))

  results <- dplyr::bind_rows(national, region_results)
  saveRDS(results, output_rds)

  region_map <- regions_sf |>
    dplyr::left_join(region_results, by = "nature_region")
  sf::st_write(region_map, output_gpkg, delete_dsn = TRUE, quiet = TRUE)

  invisible(results)
}
