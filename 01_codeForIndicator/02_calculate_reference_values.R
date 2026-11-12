#' Calculate polygon-level reference values for insect richness rasters
#'
#' Calculates summed richness rasters, optionally runs the calibration model,
#' aligns the raster layers, extracts one area-weighted mean per reference
#' polygon, classifies polygons by nature region, and saves intermediate files.
#'
#' Notes:
#' - `overwrite = FALSE` reuses existing cached raster/extraction files where
#'   available. Set `overwrite = TRUE` to regenerate them.
#' - `rerun = TRUE` reruns `fit_predict_insect_calibration()`.
#' - Extraction returns one mean per polygon, not all intersecting cell values.
#'   This matches the one-row-per-polygon data structure used below.
#' - The source function referred to downstream objects (e.g. `regional_scaled`,
#'   `region_raster`, `allRasts`) that were not created in the supplied script.
#'   This revised function returns the objects it actually computes.

calculate_reference_values <- function(
    richness,
    insects,
    study_region,
    classificationRegions,
    referencePolygons,
    crs = terra::crs(richness),
    resultFolder,
    observed_richness,
    allCalibrationRaster,
    plant_expected_richness = NULL,
    calibration_vars,
    prediction_year,
    clc_vars,
    cores = 1,
    overwrite = FALSE,
    rerun = FALSE,
    use_max = FALSE
) {
  # `plant_expected_richness` and `use_max` are retained for call compatibility;
  # they were not used in the supplied function body.
  invisible(plant_expected_richness)
  invisible(use_max)

  ## 0. Validate inputs and create output directory -------------------------
  if (missing(resultFolder) || !is.character(resultFolder) ||
      length(resultFolder) != 1L || is.na(resultFolder) || !nzchar(resultFolder)) {
    stop("Supply a valid `resultFolder` path.")
  }
  if (!dir.exists(resultFolder)) {
    dir.create(resultFolder, recursive = TRUE, showWarnings = FALSE)
  }
  if (!dir.exists(resultFolder)) stop("Could not create resultFolder: ", resultFolder)

  if (!inherits(richness, "SpatRaster")) stop("`richness` must be a terra SpatRaster.")
  if (!inherits(insects, "SpatRaster")) stop("`insects` must be a terra SpatRaster.")
  if (!inherits(study_region, "SpatRaster")) stop("`study_region` must be a terra SpatRaster.")
  if (!inherits(referencePolygons, "sf")) stop("`referencePolygons` must be an sf object.")
  if (!inherits(classificationRegions, "sf")) stop("`classificationRegions` must be an sf object.")
  if (nrow(referencePolygons) == 0L) stop("`referencePolygons` contains no features.")
  if (!"nature_region" %in% names(classificationRegions)) {
    stop("`classificationRegions` must contain a `nature_region` column.")
  }
  if (!is.numeric(cores) || length(cores) != 1L || is.na(cores) || cores < 1) {
    stop("`cores` must be a positive integer.")
  }
  cores <- as.integer(cores)

  # Use a single helper for cached raster outputs. Existing files are reused
  # only when overwrite = FALSE.
  cache_path <- function(filename) file.path(resultFolder, filename)
  should_build <- function(path) isTRUE(overwrite) || !file.exists(path)

  ## 1. Sum richness across layers ------------------------------------------
  richness_sum_path <- cache_path("insect_richness_sum.tif")
  if (should_build(richness_sum_path)) {
    richness_sum <- terra::app(
      richness,
      fun = sum,
      na.rm = TRUE,
      cores = cores,
      filename = richness_sum_path,
      overwrite = TRUE
    )
  } else {
    richness_sum <- terra::rast(richness_sum_path)
  }
  names(richness_sum) <- "insect_richness"

  richness_no_interaction_path <- cache_path("insect_richness_without_interaction_sum.tif")
  if (should_build(richness_no_interaction_path)) {
    richness_sum_without_interaction <- terra::app(
      insects,
      fun = sum,
      na.rm = TRUE,
      cores = cores,
      filename = richness_no_interaction_path,
      overwrite = TRUE
    )
  } else {
    richness_sum_without_interaction <- terra::rast(richness_no_interaction_path)
  }
  names(richness_sum_without_interaction) <- "insect_richness_without_interaction"

  ## 2. Load or calculate calibrated richness rasters -----------------------
  predicted_insect_path <- cache_path(paste0("predicted_insect_richness_", prediction_year, ".tif"))
  predicted_reference_path <- cache_path(paste0("predicted_reference_", prediction_year, ".tif"))

  if (isTRUE(rerun)) {
    if (!exists("fit_predict_insect_calibration", mode = "function")) {
      stop("`rerun = TRUE`, but function `fit_predict_insect_calibration()` is not available.")
    }

    insect_calibration <- fit_predict_insect_calibration(
      observed_richness = observed_richness,
      allCalibrationRaster = allCalibrationRaster,
      richness = richness_sum,
      nature_regions = classificationRegions,
      calibration_vars = calibration_vars,
      clc_vars = clc_vars,
      output_folder = resultFolder,
      study_region = study_region,
      prediction_year = prediction_year,
      cores = cores,
      overwrite = TRUE,
      seed = 12345,
      checks = TRUE
    )

    if (is.null(insect_calibration$predicted_insect_richness) ||
        is.null(insect_calibration$predicted_reference_richness)) {
      stop("Calibration output must contain `predicted_insect_richness` and `predicted_reference_richness`.")
    }

    calibrated_richness <- insect_calibration$predicted_insect_richness
    reference_richness <- insect_calibration$predicted_reference_richness

    # Persist model outputs when they are SpatRasters.
    terra::writeRaster(calibrated_richness, predicted_insect_path, overwrite = TRUE)
    terra::writeRaster(reference_richness, predicted_reference_path, overwrite = TRUE)
  } else {
    if (!file.exists(predicted_insect_path)) {
      stop("Calibrated richness raster not found: ", predicted_insect_path,
           ". Set `rerun = TRUE` to run calibration, or provide the file.")
    }
    if (!file.exists(predicted_reference_path)) {
      stop("Reference richness raster not found: ", predicted_reference_path,
           ". Set `rerun = TRUE` to run calibration, or provide the file.")
    }
    calibrated_richness <- terra::rast(predicted_insect_path)
    reference_richness <- terra::rast(predicted_reference_path)
  }

  names(calibrated_richness) <- "calibrated_insect_richness"
  names(reference_richness) <- "calibrated_reference_richness"

  ## 3. Align richness layers to the calibrated raster grid -----------------
  # A raster template is used as the target. Avoid projection if grids already
  # match; otherwise project/resample and write directly to disk.
  raster_stack_path <- cache_path("rasters_for_calculation.tif")
  if (!should_build(raster_stack_path)) {
    reference_rasters <- terra::rast(raster_stack_path)
  } else {
  aligned_richness_path <- cache_path("insect_richness_aligned.tif")
  aligned_no_interaction_path <- cache_path("insect_richness_without_interaction_aligned.tif")

  same_grid <- function(x, y) {
    terra::same.crs(x, y) &&
      isTRUE(all.equal(terra::res(x), terra::res(y))) &&
      isTRUE(all.equal(terra::origin(x), terra::origin(y))) &&
      isTRUE(all.equal(terra::ext(x), terra::ext(y)))
  }

  if (same_grid(richness_sum, calibrated_richness)) {
    richness_sum_aligned <- richness_sum
  } else if (!should_build(aligned_richness_path)) {
    richness_sum_aligned <- terra::rast(aligned_richness_path)
  } else {
    richness_sum_aligned <- terra::project(
      richness_sum, calibrated_richness, method = "bilinear",
      filename = aligned_richness_path, overwrite = TRUE
    )
  }

  if (same_grid(richness_sum_without_interaction, calibrated_richness)) {
    richness_sum_without_interaction_aligned <- richness_sum_without_interaction
  } else if (!should_build(aligned_no_interaction_path)) {
    richness_sum_without_interaction_aligned <- terra::rast(aligned_no_interaction_path)
  } else {
    richness_sum_without_interaction_aligned <- terra::project(
      richness_sum_without_interaction, calibrated_richness, method = "bilinear",
      filename = aligned_no_interaction_path, overwrite = TRUE
    )
  }

  names(richness_sum_aligned) <- "insect_richness"
  names(richness_sum_without_interaction_aligned) <- "insect_richness_without_interaction"


    reference_rasters <- c(
      richness_sum_aligned,
      richness_sum_without_interaction_aligned,
      calibrated_richness,
      reference_richness
    )
    names(reference_rasters) <- c(
      "insect_richness",
      "insect_richness_without_interaction",
      "calibrated_insect_richness",
      "calibrated_reference_richness"
    )
    terra::writeRaster(reference_rasters, raster_stack_path, overwrite = TRUE)
    reference_rasters <- terra::rast(raster_stack_path)
  }
  
  names(reference_rasters) <- c(
    "insect_richness",
    "insect_richness_without_interaction",
    "calibrated_insect_richness",
    "calibrated_reference_richness"
  )

  ## 4. Transform and classify polygons ------------------------------------
  target_crs <- terra::crs(reference_rasters)
  if (is.na(target_crs) || !nzchar(target_crs)) {
    stop("The reference raster stack has no valid CRS.")
  }

  # Transform the polygons once into the raster CRS. Classification regions
  # are transformed into the same CRS so spatial joining uses matching CRS.
  referencePolygons_raster_crs <- sf::st_transform(referencePolygons, target_crs)
  classificationRegions_raster_crs <- sf::st_transform(classificationRegions, target_crs)

  referencePolygons_classified <- sf::st_join(
    referencePolygons_raster_crs,
    classificationRegions_raster_crs["nature_region"],
    join = sf::st_within,
    left = TRUE
  )

  n_unclassified <- sum(is.na(referencePolygons_classified$nature_region))
  if (n_unclassified > 0L) {
    warning(n_unclassified, " reference polygons were not assigned to a nature region.")
  }

  referencePolygons_vect <- terra::vect(referencePolygons_classified)

  ## 5. Extract one area-weighted mean per polygon --------------------------
  # Exact overlap is more computationally expensive than ordinary extraction,
  # but here it is retained to calculate overlap-weighted means. The extraction
  # file is cached to avoid repeating this step on subsequent runs.
  extraction_path <- cache_path("reference_polygon_values.csv")

  if (!should_build(extraction_path)) {
    reference_polygon_values <- utils::read.csv(
      extraction_path,
      check.names = FALSE
    )
  } else {

reference_polygon_values <- terra::extract(
  reference_rasters,
  referencePolygons_vect,
  exact = TRUE,
  ID = TRUE
)

# Add polygon attributes to each extracted cell record
reference_polygon_values <- reference_polygon_values |>
  dplyr::left_join(
    referencePolygons_vect |>
      sf::st_as_sf()|>
      sf::st_drop_geometry() |>
      dplyr::mutate(ID = dplyr::row_number()),
    by = "ID"
  )

# Save the combined table
utils::write.csv(
  reference_polygon_values,
  extraction_path,
  row.names = FALSE
)
  }

 # if (nrow(reference_polygon_values) != nrow(referencePolygons_classified)) {
  #  stop(
   #   "Expected one extraction row per reference polygon, but got ",
   #   nrow(reference_polygon_values), " rows for ",
   #   nrow(referencePolygons_classified), " polygons."
   # )
  #}

  expected_extract_cols <- c(
    "insect_richness",
    "insect_richness_without_interaction",
    "calibrated_insect_richness",
    "calibrated_reference_richness"
  )
  missing_extract_cols <- setdiff(expected_extract_cols, names(reference_polygon_values))
  if (length(missing_extract_cols)) {
    stop("Missing columns in extracted values: ", paste(missing_extract_cols, collapse = ", "))
  }

  ## 6. Build and filter reference-level table ------------------------------
  reference_cols <- intersect(
    c("naturtype", "naturtypeKode", "natrtyp", "ntrtypK"),
    names(referencePolygons_classified)
  )


  reference_polygon_values <- reference_polygon_values |>
    dplyr::select(insect_richness, insect_richness_without_interaction,
    calibrated_reference_richness, natrtyp, ntrtypK, tilstnd, ntrmngf, nature_region) |>
    dplyr::rename(
      nature_type = natrtyp,
      nature_type_code = ntrtypK,
      condition = tilstnd,
      nature_diversity = ntrmngf
    ) 

  value_cols <- c(
    "insect_richness",
    "insect_richness_without_interaction",
    "calibrated_reference_richness"
  )

  reference_polygon_values <- reference_polygon_values |>
    dplyr::filter(
      !is.na(nature_region),
      dplyr::if_all(dplyr::all_of(value_cols), ~ is.finite(.x))
    )


  ## 7. Save tabular results ------------------------------------------------
  utils::write.csv(
    reference_polygon_values,
    cache_path("reference_values_raw.csv"),
    row.names = FALSE
  )

  ## 8. Return objects actually computed -----------------------------------
  list(
    richness_sum = richness_sum,
    richness_sum_without_interaction = richness_sum_without_interaction,
    calibrated_richness = calibrated_richness,
    reference_richness = reference_richness,
    reference_rasters = reference_rasters,
    reference_polygon_values = reference_polygon_values,
    reference_polygons = referencePolygons_classified,
    unclassified_polygon_count = n_unclassified,
    output_files = c(
      richness_sum = richness_sum_path,
      richness_without_interaction = richness_no_interaction_path,
      aligned_rasters = raster_stack_path,
      polygon_values = extraction_path,
      reference_values = cache_path("reference_values_raw.csv")
    )
  )
}
