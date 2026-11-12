calculate_scaled_indicator <- function(
    classificationRegions,
    overwrite = TRUE,
    rerun = FALSE,
    use_max = FALSE,
    allRasts = NULL,
    ref_vals,
    output_dir,
    value_layer = NULL,
    truncation = TRUE,
    cores = 1
) {
  # Dependencies used: terra, sf, ecTools

  ## 0. Validate inputs ------------------------------------------------------

  if (!inherits(classificationRegions, "sf")) {
    stop("classificationRegions must be an sf object.")
  }
  if (is.null(allRasts) || !inherits(allRasts, "SpatRaster")) {
    stop("allRasts must be a terra SpatRaster containing the raster to normalise.")
  }
  if (!is.data.frame(ref_vals)) {
    stop("ref_vals must be a data.frame.")
  }

  if(!is.null(ref_vals$X_60)){
  required_ref_cols <- c("nature_region", "X_0", "X_60", "X_100", "region_id")
  missing_ref_cols <- setdiff(required_ref_cols, names(ref_vals))
  if (length(missing_ref_cols)) {
    stop("ref_vals is missing required columns: ", paste(missing_ref_cols, collapse = ", "))
  }
  if (anyDuplicated(ref_vals$region_id) || anyDuplicated(as.character(ref_vals$nature_region))) {
    stop("ref_vals must contain exactly one row per region_id and nature_region.")
  }
  if (!is.numeric(ref_vals$X_0) || !is.numeric(ref_vals$X_60) || !is.numeric(ref_vals$X_100)) {
    stop("X_0, X_60, and X_100 in ref_vals must be numeric.")
  }
    
  } else {

  required_ref_cols <- c("nature_region", "X_0", "X_100", "region_id")
  missing_ref_cols <- setdiff(required_ref_cols, names(ref_vals))
  if (length(missing_ref_cols)) {
    stop("ref_vals is missing required columns: ", paste(missing_ref_cols, collapse = ", "))
  }
  if (anyDuplicated(ref_vals$region_id) || anyDuplicated(as.character(ref_vals$nature_region))) {
    stop("ref_vals must contain exactly one row per region_id and nature_region.")
  }
  if (!is.numeric(ref_vals$X_0) || !is.numeric(ref_vals$X_100)) {
    stop("X_0, and X_100 in ref_vals must be numeric.")
  }
  }
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  # Generic layer selection: if multiple layers exist, require the caller to
  # specify which one to normalise instead of silently choosing the wrong layer.
  raster_names <- names(allRasts)
  if (is.null(value_layer)) {
    if (terra::nlyr(allRasts) != 1L) {
      stop(
        "allRasts contains multiple layers. Specify value_layer explicitly. Available layers: ",
        paste(raster_names, collapse = ", ")
      )
    }
    value_layer <- raster_names[1]
  }
  if (length(value_layer) != 1L || is.na(value_layer) || !value_layer %in% raster_names) {
    stop(
      "value_layer must be one layer name in allRasts. Available layers: ",
      paste(raster_names, collapse = ", ")
    )
  }

  # Validate reference anchors. ec_normalise supports vector-valued anchors,
  # but each region's three anchor values must be ordered and finite.
  has_x60 <- !is.null(ref_vals$X_60)
  anchor_cols <- if (has_x60) c("X_0", "X_60", "X_100") else c("X_0", "X_100")
  anchors <- as.matrix(ref_vals[, anchor_cols])
  storage.mode(anchors) <- "double"
  valid_anchors <- apply(anchors, 1L, function(z) {
    all(is.finite(z)) && all(diff(z) > 0)
  })
  if (!all(valid_anchors)) {
    stop("Invalid anchors in ref_vals for: ",
         paste(as.character(ref_vals$nature_region[!valid_anchors]), collapse = ", "))
  }

  # Build safe, layer-specific output names.
  safe_layer <- gsub("[^A-Za-z0-9_-]+", "_", value_layer)
  output_file <- file.path(output_dir, paste0(safe_layer, "_normalised.tif"))
  region_file <- file.path(output_dir, "nature_region_id.tif")

  # Reuse output only when explicitly requested via overwrite = FALSE and
  # rerun = FALSE. rerun = TRUE forces regeneration.
  if (file.exists(output_file) && !overwrite && !rerun) {
    normalised_rast <- terra::rast(output_file)
    region_raster <- if (file.exists(region_file)) terra::rast(region_file) else NULL
    return(list(
      allRasts = allRasts,
      input_layer = value_layer,
      region_raster = region_raster,
      scaled_richness = normalised_rast,
      scaled_richness_masked = normalised_rast,
      scaled_richness_clamped = normalised_rast,
      output_file = output_file,
      reference_values = ref_vals
    ))
  }

  ## 1. Prepare region polygons and rasterise IDs ----------------------------
  if (!"nature_region" %in% names(classificationRegions)) {
    stop("classificationRegions must contain a 'nature_region' column.")
  }
  if (is.na(sf::st_crs(classificationRegions))) {
    stop("classificationRegions has no defined CRS.")
  }

  classificationRegions <- sf::st_make_valid(classificationRegions)
  value_rast <- allRasts[[value_layer]]
  target_crs <- terra::crs(value_rast)
  if (is.na(target_crs) || !nzchar(target_crs)) {
    stop("The selected input raster has no defined CRS.")
  }
  classificationRegions <- sf::st_transform(classificationRegions, target_crs)

  region_names <- unique(as.character(classificationRegions$nature_region))
  unknown_regions <- setdiff(region_names, as.character(ref_vals$nature_region))
  if (length(unknown_regions)) {
    stop("classificationRegions contains regions missing from ref_vals: ",
         paste(unknown_regions, collapse = ", "))
  }

  classificationRegions$region_id <- ref_vals$region_id[
    match(as.character(classificationRegions$nature_region),
          as.character(ref_vals$nature_region))
  ]
  regions_vect <- terra::vect(classificationRegions)

  region_raster <- terra::rasterize(
    regions_vect,
    value_rast,
    field = "region_id",
    background = NA,
    touches = TRUE,
    filename = region_file,
    overwrite = TRUE,
    wopt = list(gdal = c("COMPRESS=LZW", "TILED=YES", "BIGTIFF=IF_SAFER"))
  )
  names(region_raster) <- "region_id"


  ## 3. Normalise in a single raster pass using ecTools::ec_normalise --------
  # ec_normalise accepts vector-valued anchors (as in the package examples).
  # We only pass cells with valid value, region ID, mask, and matched anchors,
  # avoiding NA anchors that can trigger ec_normalise's ordering check.
  input_stack <- c(value_rast, region_raster)
  names(input_stack) <- c("value", "region_id")

  normalised_rast <- terra::lapp(
    input_stack,
    fun = function(value, region_id) {
      out <- rep(NA_real_, length(value))
      idx <- match(region_id, ref_vals$region_id)
      ok <- is.finite(value) & !is.na(idx) 

      if (any(ok)) {
        x <- value[ok]
        j <- idx[ok]
        result <- ecTools::ec_normalise(
          variable = x,
          x0 = ref_vals$X_0[j],
          x60 = if (has_x60) ref_vals$X_60[j] else NULL,
          x100 = ref_vals$X_100[j],
          truncation = truncation
        )
        if (length(result) != length(x)) {
          stop(
            "ecTools::ec_normalise returned ", length(result),
            " values for ", length(x), " input cells. Check ec_normalise vectorisation."
          )
        }
        out[ok] <- as.numeric(result)
      }
      out
    },
    cores = cores,
    filename = output_file,
    overwrite = TRUE,
    wopt = list(gdal = c("COMPRESS=LZW", "TILED=YES", "BIGTIFF=IF_SAFER"))
  )
  names(normalised_rast) <- paste0(value_layer, "_normalised")

  list(
    allRasts = allRasts,
    input_layer = value_layer,
    region_raster = region_raster,
    scaled_richness = normalised_rast,
    scaled_richness_masked = normalised_rast,
    scaled_richness_clamped = normalised_rast,
    output_file = output_file,
    reference_values = ref_vals
  )
  }
