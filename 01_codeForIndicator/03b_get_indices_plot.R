################################################################################
# Indicator estimation using fixed samples and bootstrap resampling
#
# National and regional estimates use the ORIGINAL raster resolution.
# For each area, sample up to 100,000 complete cells WITHOUT replacement.
# Then bootstrap 10,000 observations WITH replacement from that fixed sample.
# The national beta mixed model includes a random intercept for Nature region;
# regional models remain intercept-only beta regressions.
#
# Required objects in the calling environment:
#   plotFolder    : folder containing all_indicator_raster.tif
#   nature_regions: sf object with a nature_region column and valid CRS
#
# Required raster layers:
#   regional_clamped, calibrated_regional_clamped,
#   regional_clamped_without_interaction, global_clamped,
#   calibrated_global_clamped, global_clamped_without_interaction
################################################################################

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(dplyr)
  library(purrr)
  library(tibble)
  library(betareg)
  library(glmmTMB)
})

# --------------------------- Settings ---------------------------------------
B <- 500L
initial_sample_size <- 500000L   # initial sample: without replacement
bootstrap_sample_size <- 10000L # bootstrap sample: with replacement
bootstrap_seed <- 1000L
threshold <- 0.6
force_rebuild_mask <- TRUE      # TRUE rebuilds the mountain mask

stopifnot(exists("plotFolder"), dir.exists(plotFolder))
stopifnot(exists("nature_regions"))

#raster_file <- file.path(plotFolder, "all_indicator_raster.tif")
raster_file <- file.path(plotFolder, "calibrated_insect_richness_normalised.tif")
masked_file <- file.path(plotFolder, "all_indicator_raster_no_mountains.tif")
results_file <- file.path(plotFolder, "region_results.rds")
samples_file <- file.path(plotFolder, "indicator_fixed_samples.rds")
bootstrap_file <- file.path(plotFolder, "indicator_bootstrap_draws.rds")
map_file <- file.path(plotFolder, "region_map.gpkg")

indicator_names <- c(
  #"Regional indicator",
  "Calibrated regional indicator"#,
  #"Regional indicator without interaction",
  #"Global indicator",
 # "Calibrated global indicator",
  #"Global indicator without interaction"
)
indicator_layer_names <- c(
 # "regional_clamped",
  "calibrated_insect_richness_normalised"#,
 # "regional_clamped_without_interaction",
 # "global_clamped",
 # "calibrated_global_clamped",
 # "global_clamped_without_interaction"
)

if (!file.exists(raster_file)) stop("Raster not found: ", raster_file)
allRasts <- terra::rast(raster_file)
names(allRasts) <- indicator_layer_names
if (!all(indicator_layer_names %in% names(allRasts))) {
  stop("Missing indicator layer(s): ",
       paste(setdiff(indicator_layer_names, names(allRasts)), collapse = ", "))
}

# --------------------- Exclude mountain areas -------------------------------
if (file.exists(masked_file) && !force_rebuild_mask) {
  message("Using existing masked raster: ", masked_file)
  allRasts_masked <- terra::rast(masked_file)
} else {
  mountain_file <- file.path(
    "R:/GeoSpatialData", "LandCover", "Norway_MountainAreas",
    "Original", "Fjell_norge_publisering.gdb"
  )
  if (!file.exists(mountain_file)) {
    stop("Mountain geodatabase not found: ", mountain_file,
         ". If the masked raster already exists, set force_rebuild_mask <- FALSE.")
  }
  fjell_grid50 <- sf::st_read(
    dsn = mountain_file, layer = "fjell_grid50", quiet = TRUE
  )
  mountain_mask <- terra::rasterize(
    terra::vect(fjell_grid50), allRasts[[1]], field = 1, background = 0
  )
  allRasts_masked <- terra::mask(
    allRasts, mountain_mask, maskvalues = 1, updatevalue = NA
  )
  terra::writeRaster(
    allRasts_masked, masked_file, overwrite = TRUE,
    gdal = c("COMPRESS=LZW", "TILED=YES")
  )
}

if (!all(indicator_layer_names %in% names(allRasts_masked))) {
  stop("Masked raster is missing one or more required indicator layers.")
}
indicator_rasters <- allRasts[[indicator_layer_names]]
names(indicator_rasters) <- indicator_names

# ----------------------- Prepare regions ------------------------------------
nature_regions <- sf::st_as_sf(nature_regions)
if (!"nature_region" %in% names(nature_regions)) {
  stop("nature_regions must contain a 'nature_region' column.")
}
if (is.na(sf::st_crs(nature_regions))) {
  stop("nature_regions has no CRS. Assign its correct CRS first.")
}
nature_regions <- sf::st_transform(nature_regions, crs = terra::crs(indicator_rasters))
regions <- unique(as.character(nature_regions$nature_region))
regions <- regions[!is.na(regions)]
if (!length(regions)) stop("No valid Nature-region names found.")

# ----------------------- Sampling helpers -----------------------------------
area_raster <- function(r, polygon = NULL) {
  if (is.null(polygon)) return(r)
  p <- terra::vect(polygon)
  terra::mask(terra::crop(r, p), p)
}

# Draw up to n valid cells without replacement.
# Avoid counting all valid raster cells first, because that requires a full scan.
# exhaustive=TRUE asks terra to keep searching when missing values are common.
sample_fixed_values <- function(
    r, n = initial_sample_size, seed = 1000L, label = "area"
) {
  message(
    label, ": sampling up to ", format(n, big.mark = ","),
    " valid cells without replacement."
  )

  set.seed(seed)
  #x <- terra::spatSample(
  #  r,
   # size = as.integer(n),
   # method = "regular",
   # replace = FALSE,
   # na.rm = TRUE,
   # as.df = TRUE,
   # xy = FALSE,
   # exhaustive = TRUE
  #)

  x <- terra::values(r, na.rm = TRUE)[, 1]

  x <- as.data.frame(x)
  x <- x[stats::complete.cases(x), , drop = FALSE]

  if (!nrow(x)) stop("Sampling returned no complete rows for ", label, ".")
  if (ncol(x) != length(indicator_names)) {
    stop(
      "Unexpected number of sampled columns for ", label,
      ": expected ", length(indicator_names), ", got ", ncol(x), "."
    )
  }

  names(x) <- indicator_names
  message(label, ": obtained ", format(nrow(x), big.mark = ","), " rows.")
  x
}

# Beta regression requires responses strictly between 0 and 1. If boundary
# values occur, apply the Smithson-Verkuilen adjustment to that fixed sample.
adjust_beta_response <- function(y) {
  if (any(!is.finite(y))) stop("Non-finite indicator values encountered.")
  if (any(y < 0 | y > 1)) stop("Indicator values outside [0, 1] encountered.")
  n <- length(y)
  if (n < 2L) stop("At least two observations are required for beta regression.")
  if (any(y <= 0 | y >= 1)) y <- (y * (n - 1) + 0.5) / n
  y
}

# Fit beta regression for a fixed sample. The national model includes a
# random intercept for Nature region; regional models remain intercept-only.
estimate_beta_mean <- function(y, region = NULL, national = FALSE) {
  if (any(!is.finite(y))) stop("Non-finite indicator values encountered.")
  if (any(y < 0 | y > 1)) stop("Indicator values outside [0, 1] encountered.")
  if (length(y) < 2L) stop("At least two observations are required.")

  if (all(y == 1)) return(1)
  if (all(y == 0)) return(0)

  n <- length(y)
  if (any(y <= 0 | y >= 1)) {
    y <- (y * (n - 1) + 0.5) / n
  }

  if (national) {
    if (is.null(region)) {
      stop("Region labels are required for the national beta mixed model.")
    }
    dat <- data.frame(y = y, nature_region = droplevels(factor(region)))
    fit <- tryCatch(
      glmmTMB::glmmTMB(
        y ~ 1 + (1 | nature_region),
        data = dat,
        family = glmmTMB::beta_family(link = "logit")
      ),
      error = function(e) NULL
    )
    if (is.null(fit)) return(mean(y))

    # This is the inverse-logit fixed intercept: the fitted mean at a
    # zero random effect, not automatically the marginal national mean.
    return(stats::plogis(glmmTMB::fixef(fit)$cond[["(Intercept)"]]))
  }

  dat <- data.frame(y = y)
  fit <- tryCatch(
    betareg::betareg(y ~ 1, data = dat, link = "logit"),
    error = function(e) NULL
  )
  if (is.null(fit)) return(mean(y))
  unname(stats::plogis(stats::coef(fit, model = "mean")[["(Intercept)"]]))
}

# Bootstrap from a fixed sample: resample rows with replacement. For the
# national model, keep the region label attached to each resampled observation.
bootstrap_fixed_sample <- function(x, B = 100L, bootstrap_n = 10000L,
                                   seed = 1000L, label = "area",
                                   national = FALSE) {
  if (nrow(x) < 2L) stop("Too few sampled cells for ", label)
  if (national && !"nature_region" %in% names(x)) {
    stop("National sample must contain a 'nature_region' column.")
  }

  ylist <- lapply(indicator_names, function(nm) adjust_beta_response(x[[nm]]))
  names(ylist) <- indicator_names
  out <- vector("list", B)

  for (b in seq_len(B)) {
    if (b == 1L || b %% 10L == 0L || b == B) {
      message("  ", label, " bootstrap ", b, "/", B)
    }
    set.seed(seed + b)
    idx <- sample.int(nrow(x), size = bootstrap_n, replace = TRUE)
    sampled_regions <- if (national) x$nature_region[idx] else NULL

    est <- purrr::map_dbl(indicator_names, function(nm) {
      estimate_beta_mean(
        ylist[[nm]][idx],
        region = sampled_regions,
        national = national
      )
    })

    out[[b]] <- tibble::tibble(
      replicate = b, indicator = indicator_names, estimate = est
    )
  }
  dplyr::bind_rows(out)
}

summarise_bootstrap <- function(boot) {
  boot %>%
    group_by(indicator) %>%
    summarise(
      bootstrap_mean = mean(estimate, na.rm = TRUE),
      bootstrap_sd = sd(estimate, na.rm = TRUE),
      lower_quantile = unname(quantile(estimate, 0.025, na.rm = TRUE)),
      upper_quantile = unname(quantile(estimate, 0.975, na.rm = TRUE)),
      n_successful_bootstraps = sum(is.finite(estimate)),
      .groups = "drop"
    )
}


# Calculate regional estimates from a fixed sample.
# To reduce runtime, descriptive statistics are calculated from the sample
# instead of scanning every raster cell for exact full-raster summaries.
calculate_area <- function(r, polygon = NULL, label = "region", seed = 1000L) {
  message("\nProcessing region: ", label)
  rr <- area_raster(r, polygon)

  fixed_sample <- sample_fixed_values(
    rr, n = initial_sample_size, seed = seed, label = label
  )
  fixed_sample$nature_region <- label

  boot <- bootstrap_fixed_sample(
    fixed_sample, B = B, bootstrap_n = bootstrap_sample_size,
    seed = seed + 500000L, label = label, national = FALSE
  )
  boot_summary <- summarise_bootstrap(boot)

  descriptive <- purrr::map_dfr(indicator_names, function(nm) {
    x <- fixed_sample[[nm]]
    tibble::tibble(
      indicator = nm,
      global_mean = mean(x),
      global_sd = stats::sd(x),
      proportion_ge_0.6 = mean(x >= threshold),
      percentage_ge_0.6 = mean(x >= threshold) * 100
    )
  })

  results <- descriptive %>%
    dplyr::left_join(boot_summary, by = "indicator") %>%
    dplyr::mutate(
      nature_region = label,
      n_valid_cells = NA_real_,
      n_sampled_cells = nrow(fixed_sample)
    ) %>%
    dplyr::select(
      nature_region, indicator,
      global_mean, global_sd,
      proportion_ge_0.6, percentage_ge_0.6,
      n_valid_cells, n_sampled_cells,
      bootstrap_mean, bootstrap_sd,
      lower_quantile, upper_quantile,
      n_successful_bootstraps
    )

  list(results = results, sample = fixed_sample, bootstrap = boot)
}

if (!file.exists(map_file
)) {

# ------------------------ Regional estimates FIRST ----------------------------
# Run regional sampling first. The national sample will be pooled from these
# already-created samples, avoiding a separate national raster scan/sample.
region_outputs <- setNames(vector("list", length(regions)), regions)

for (i in seq_along(regions)) {
  reg <- regions[[i]]
  poly <- nature_regions %>% dplyr::filter(as.character(nature_region) == reg)

  if (!nrow(poly)) {
    warning("No polygon found for region: ", reg)
    next
  }

  region_outputs[[reg]] <- calculate_area(
    indicator_rasters, polygon = poly, label = reg,
    seed = bootstrap_seed + i * 10000L
  )
}

region_outputs <- region_outputs[!vapply(region_outputs, is.null, logical(1))]
if (!length(region_outputs)) {
  stop("No regional samples were created; cannot calculate the national estimate.")
}

region_results <- purrr::map_dfr(region_outputs, "results")

# ------------------------ National from regional samples ----------------------
# Pool the existing samples. No new national raster sampling or full-raster
# national statistics are calculated.
regional_samples <- lapply(region_outputs, `[[`, "sample")
national_sample <- dplyr::bind_rows(regional_samples)
national_sample <- national_sample[
  stats::complete.cases(national_sample[, indicator_names, drop = FALSE]) &
    !is.na(national_sample$nature_region), , drop = FALSE
]

if (nrow(national_sample) < 2L) {
  stop("The pooled regional sample contains fewer than two complete rows.")
}

message(
  "\nNational estimate: pooled ",
  format(nrow(national_sample), big.mark = ","),
  " rows from regional samples."
)

national_bootstrap <- bootstrap_fixed_sample(
  national_sample, B = B, bootstrap_n = bootstrap_sample_size,
  seed = bootstrap_seed + 500000L, label = "National", national = TRUE
)
national_bootstrap_summary <- summarise_bootstrap(national_bootstrap)

national_descriptive <- purrr::map_dfr(indicator_names, function(nm) {
  x <- national_sample[[nm]]
  tibble::tibble(
    indicator = nm,
    global_mean = mean(x),
    global_sd = stats::sd(x),
    proportion_ge_0.6 = mean(x >= threshold),
    percentage_ge_0.6 = mean(x >= threshold) * 100
  )
})

national_results <- national_descriptive %>%
  dplyr::left_join(national_bootstrap_summary, by = "indicator") %>%
  dplyr::mutate(
    nature_region = "National",
    n_valid_cells = NA_real_,
    n_sampled_cells = nrow(national_sample)
  ) %>%
  dplyr::select(
    nature_region, indicator,
    global_mean, global_sd,
    proportion_ge_0.6, percentage_ge_0.6,
    n_valid_cells, n_sampled_cells,
    bootstrap_mean, bootstrap_sd,
    lower_quantile, upper_quantile,
    n_successful_bootstraps
  )

allIndicatorsEstimated <- bind_rows(national_results, region_results)
saveRDS(allIndicatorsEstimated, results_file)

# Save fixed samples and bootstrap draws to allow auditing/re-summarising.
saveRDS(
  list(
    national = national_sample,
    regions = lapply(region_outputs, `[[`, "sample"),
    settings = list(
      B = B, initial_sample_size = initial_sample_size,
      bootstrap_sample_size = bootstrap_sample_size,
      seed = bootstrap_seed, threshold = threshold,
      raster_resolution = terra::res(indicator_rasters),
      source_raster = raster_file, masked_raster = masked_file,
      national_method = "Pooled regional samples; contribution is proportional to regional sample size"
    )
  ),
  samples_file
)
saveRDS(
  list(
    national = national_bootstrap,
    regions = lapply(region_outputs, `[[`, "bootstrap")
  ),
  bootstrap_file
)

# Save regional results as GeoPackage (preserves long field names better than
# ESRI Shapefile). Results remain long-form: one row per region and indicator.
region_map <- nature_regions %>% left_join(region_results, by = "nature_region")
if (file.exists(map_file)) unlink(map_file)
sf::st_write(region_map, map_file, quiet = TRUE)

message("\nCompleted.")
message("Results: ", results_file)
message("Fixed samples: ", samples_file)
message("Bootstrap draws: ", bootstrap_file)
message("Regional map: ", map_file)
message("National bootstrap uses pooled regional samples; descriptive statistics are sample-based.")
}
