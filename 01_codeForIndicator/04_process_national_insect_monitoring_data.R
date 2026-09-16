################################################################################
# Process national insect monitoring data
#
# This script:
#
#   1. Checks whether the processed national insect monitoring dataset already
#      exists.
#
#   2. If necessary, downloads the national insect monitoring dataset from GBIF.
#
#   3. Reads the hierarchical monitoring data:
#
#          Year / locality
#              │
#              └── Locality sampling event
#                    │
#                    └── Sampling trap
#                          │
#                          └── Identification
#                                │
#                                └── Occurrence
#
#   4. Extracts and joins the information from these different levels.
#
#   5. Converts insect occurrence coordinates into spatial points.
#
#   6. Retains insect species represented in the insect SDM.
#
#   7. Restricts the data to semi-natural habitats.
#
#   8. Calculates observed insect species richness for each locality and year.
#
#   9. Calculates associated sampling/environmental variables:
#
#          - minimum sampling temperature
#          - maximum sampling temperature
#          - mean trap net wet weight
#          - sampling duration
#
#  10. Joins these variables back to the spatial monitoring locations.
#
#  11. Saves the resulting dataset as a Shapefile.
#
# Final object:
#
#     national_monitoring_insect_richness
#
################################################################################


# ==============================================================================
# 0. Check whether the processed output already exists
# ==============================================================================
#
# The filename contains `var_of_interest`, allowing different versions of the
# processed dataset to be saved separately if required.
# ==============================================================================

national_monitoring_file <- file.path(
  resultFolder,
  "polygons",
  paste0(
    "ninaInsect_monitoring_with_richness_for_",
    var_of_interest,
    ".shp"
  )
)


if (!file.exists(national_monitoring_file)) {
  
  
  # ============================================================================
  # 1. Download and process the national insect monitoring data
  # ============================================================================
  #
  # The raw data are downloaded from the GBIF endpoint associated with the
  # National Insect Monitoring dataset.
  #
  # The original workflow/reference for formatting the data is provided in
  # the NINA national insect monitoring repository.
  #
  # The raw GBIF export is only downloaded and processed if the intermediate
  # file `out/nationalInsectMonitoring.csv` does not already exist.
  # ============================================================================
  
  if (!file.exists("out/nationalInsectMonitoring.csv")) {
    
    
    # --------------------------------------------------------------------------
    # 1.1 Load packages required for processing the GBIF metadata/export
    # --------------------------------------------------------------------------
    
    suppressPackageStartupMessages({
      
      require(tidyverse)
      require(tidyjson)
      require(xml2)
      
    })
    
    
    # ==========================================================================
    # 2. Identify the GBIF dataset
    # ==========================================================================
    #
    # `dataset_id` identifies the National Insect Monitoring dataset in GBIF.
    # ==========================================================================
    
    dataset_id <- "19fe96b0-0cf3-4a2e-90a5-7c1c19ac94ee"
    
    
    # ==========================================================================
    # 3. Retrieve the GBIF citation metadata
    # ==========================================================================
    #
    # The GBIF API provides the citation associated with the dataset.
    #
    # This is useful for preserving the appropriate dataset citation when the
    # data are subsequently used in analyses or publications.
    # ==========================================================================
    
    tmp <- tempfile()
    
    
    # --------------------------------------------------------------------------
    # Retrieve the dataset metadata document
    # --------------------------------------------------------------------------
    
    download.file(
      paste0(
        "http://api.gbif.org/v1/dataset/",
        dataset_id,
        "/document"
      ),
      tmp
    )
    
    
    # --------------------------------------------------------------------------
    # Read the XML metadata
    # --------------------------------------------------------------------------
    
    meta <- read_xml(tmp) %>%
      as_list()
    
    
    # --------------------------------------------------------------------------
    # Extract the GBIF citation
    # --------------------------------------------------------------------------
    
    gbif_citation <-
      meta$eml$additionalMetadata$metadata$gbif$citation[[1]]
    
    
    # ==========================================================================
    # 4. Retrieve the GBIF data endpoint
    # ==========================================================================
    
    dataset_url <- paste0(
      "http://api.gbif.org/v1/dataset/",
      dataset_id,
      "/endpoint"
    )
    
    
    dataset <- RJSONIO::fromJSON(
      dataset_url
    )
    
    
    # Extract the actual download endpoint
    endpoint_url <- dataset[[1]]$url
    
    
    # ==========================================================================
    # 5. Create a directory for the GBIF download
    # ==========================================================================
    
    dir.create(
      "GBIF_data",
      showWarnings = FALSE,
      recursive = TRUE
    )
    
    
    # ==========================================================================
    # 6. Set a longer download timeout
    # ==========================================================================
    #
    # The GBIF export may be relatively large, so allow up to 10 minutes for
    # the download.
    # ==========================================================================
    
    options(
      timeout = 600
    )
    
    
    # ==========================================================================
    # 7. Download the GBIF export
    # ==========================================================================
    
    download.file(
      url = endpoint_url,
      destfile = "GBIF_data/gbif_download.zip",
      mode = "wb",
      method = "libcurl"
    )
    
    
    # ==========================================================================
    # 8. Unzip the GBIF export
    # ==========================================================================
    
    unzip(
      "GBIF_data/gbif_download.zip",
      exdir = "GBIF_data"
    )
    
    
    # ==========================================================================
    # 9. Read the event-level data
    # ==========================================================================
    #
    # `event.txt` contains information describing the different levels of the
    # monitoring sampling hierarchy.
    # ==========================================================================
    
    event_raw <- readr::read_delim(
      "GBIF_data/event.txt",
      delim = "\t",
      locale = locale(
        encoding = "UTF-8"
      ),
      progress = FALSE,
      show_col_types = FALSE,
      guess_max = 10000
    )
    
    
    # ==========================================================================
    # 10. Read the occurrence-level data
    # ==========================================================================
    #
    # `occurrence.txt` contains the taxonomic observations associated with
    # the monitoring events.
    # ==========================================================================
    
    occurrence_raw <- readr::read_delim(
      "GBIF_data/occurrence.txt",
      delim = "\t",
      locale = locale(
        encoding = "UTF-8"
      ),
      progress = FALSE,
      show_col_types = FALSE,
      guess_max = 10000
    )
    
    
    # ==========================================================================
    # 11. Separate the monitoring hierarchy
    # ==========================================================================
    #
    # The event data use `samplingProtocol` to distinguish four hierarchical
    # levels of the monitoring design.
    #
    # Level 1 = identification
    # Level 2 = sampling trap
    # Level 3 = locality sampling event
    # Level 4 = year/locality
    # ==========================================================================
    
    identifications_raw <- event_raw %>%
      filter(
        samplingProtocol == "Level_1"
      )
    
    
    sampling_trap_raw <- event_raw %>%
      filter(
        samplingProtocol == "Level_2"
      )
    
    
    locality_sampling_raw <- event_raw %>%
      filter(
        samplingProtocol == "Level_3"
      )
    
    
    year_locality_raw <- event_raw %>%
      filter(
        samplingProtocol == "Level_4"
      )
    
    
    # ==========================================================================
    # 12. Process identification-level information
    # ==========================================================================
    #
    # Dynamic properties are stored separately within the GBIF export.
    # These properties are unpacked and then joined back to the identification
    # records.
    # ==========================================================================
    
    identifications <- identifications_raw %>%
      select(
        -dynamicProperties
      )
    
    
    tempDynamic <- identifications_raw %>%
      select(
        dyn = dynamicProperties
      ) %>%
      unlist() %>%
      tidyjson::spread_all() %>%
      as_tibble() %>%
      select(
        -document.id
      )
    
    
    # --------------------------------------------------------------------------
    # Join dynamic identification properties back to the identification table
    # --------------------------------------------------------------------------
    
    identifications <- identifications %>%
      left_join(
        tempDynamic,
        by = c(
          "id" = "id"
        )
      )
    
    
    # --------------------------------------------------------------------------
    # Check that joining did not change the number of records
    #
    # If the number of rows changes, the join structure is not behaving as
    # expected and the script stops.
    # --------------------------------------------------------------------------
    
    if (
      !nrow(identifications) == nrow(identifications_raw) &
      nrow(identifications) == nrow(tempDynamic)
    ) {
      stop(
        "Identification rows don't match!"
      )
    }
    
    
    # --------------------------------------------------------------------------
    # Retain relevant identification variables
    # --------------------------------------------------------------------------
    
    identifications <- identifications %>%
      select(
        identification_id = eventID,
        sampling_trap_id = parentEventID,
        identification_name,
        identification_comment,
        read_abundance,
        ekstraksjonsdato,
        ekstraksjonskit,
        ekstraksjonskommentar
      )
    
    
    # ==========================================================================
    # 13. Process sampling-trap information
    # ==========================================================================
    
    sampling_trap <- sampling_trap_raw %>%
      
      # Extract the trap name from location remarks
      mutate(
        trap_name = stringr::str_split_i(
          locationRemarks,
          ": ",
          2
        )
      ) %>%
      
      select(
        -dynamicProperties
      )
    
    
    # --------------------------------------------------------------------------
    # Unpack dynamic trap properties
    # --------------------------------------------------------------------------
    
    tempDynamic <- sampling_trap_raw %>%
      select(
        dyn = dynamicProperties
      ) %>%
      unlist() %>%
      tidyjson::spread_all() %>%
      as_tibble() %>%
      select(
        -document.id
      )
    
    
    # --------------------------------------------------------------------------
    # Join dynamic properties to trap information
    # --------------------------------------------------------------------------
    
    sampling_trap <- sampling_trap %>%
      left_join(
        tempDynamic,
        by = c(
          "id" = "id"
        )
      )
    
    
    # --------------------------------------------------------------------------
    # Check that the join preserved the expected number of rows
    # --------------------------------------------------------------------------
    
    if (
      !nrow(sampling_trap) == nrow(sampling_trap_raw) &
      nrow(sampling_trap) == nrow(tempDynamic)
    ) {
      stop(
        "Sampling trap rows don't match!"
      )
    }
    
    
    # --------------------------------------------------------------------------
    # Retain relevant trap variables
    # --------------------------------------------------------------------------
    
    sampling_trap <- sampling_trap %>%
      select(
        sampling_trap_id = eventID,
        locality_sampling_id = parentEventID,
        trap_name,
        sample_size_value = sampleSizeValue,
        sample_size_unit = sampleSizeUnit,
        trap_type,
        trap_model,
        liquid_name,
        ethanol_prc_at_lab,
        net_wet_weight,
        
        # Trap coordinates
        point_latitude = decimalLatitude,
        point_longitude = decimalLongitude,
        
        # Spatial footprint
        point_wkt = footprintWKT,
        
        # Coordinate uncertainty
        coord_unc_in_meters = coordinateUncertaintyInMeters,
        
        # Coordinate reference system
        coordinate_system = geodeticDatum
      )
    
    
    # ==========================================================================
    # 14. Process locality-sampling information
    # ==========================================================================
    
    locality_sampling <- locality_sampling_raw %>%
      
      mutate(
        sampling_event = stringr::str_split_i(
          locationRemarks,
          ": ",
          2
        )
      ) %>%
      
      select(
        -dynamicProperties
      )
    
    
    # --------------------------------------------------------------------------
    # Unpack dynamic locality-sampling properties
    # --------------------------------------------------------------------------
    
    tempDynamic <- locality_sampling_raw %>%
      select(
        dyn = dynamicProperties
      ) %>%
      unlist() %>%
      tidyjson::spread_all() %>%
      as_tibble() %>%
      select(
        -document.id
      )
    
    
    # --------------------------------------------------------------------------
    # Join dynamic properties back to locality-sampling records
    # --------------------------------------------------------------------------
    
    locality_sampling <- locality_sampling %>%
      left_join(
        tempDynamic,
        by = c(
          "id" = "id"
        )
      )
    
    
    # --------------------------------------------------------------------------
    # Validate the join
    # --------------------------------------------------------------------------
    
    if (
      !nrow(locality_sampling) == nrow(locality_sampling_raw) &
      nrow(locality_sampling) == nrow(tempDynamic)
    ) {
      stop(
        "Locality sampling rows don't match!"
      )
    }
    
    
    # ==========================================================================
    # 15. Convert locality sampling times
    # ==========================================================================
    #
    # `eventTime` contains a start/end interval.
    #
    # These are separated into:
    #
    #     start_time
    #     end_time
    #
    # and converted to POSIXct date-time values.
    # ==========================================================================
    
    locality_sampling <- locality_sampling %>%
      
      separate(
        eventTime,
        c(
          "start_time",
          "end_time"
        ),
        sep = "/"
      ) %>%
      
      mutate(
        start_time = as.POSIXct(start_time),
        end_time = as.POSIXct(end_time)
      )
    
    
    # Optional inspection of the sampling period
    locality_sampling %>%
      select(
        start_time,
        end_time
      )
    
    
    # ==========================================================================
    # 16. Retain relevant locality-sampling variables
    # ==========================================================================
    
    locality_sampling <- locality_sampling %>%
      select(
        locality_sampling_id = eventID,
        year_locality_id = parentEventID,
        sampling_event,
        sample_size_value = sampleSizeValue,
        sample_size_unit = sampleSizeUnit,
        start_time,
        end_time,
        sampling_min_temp,
        sampling_max_temp,
        sampling_avg_temp
      )
    
    
    # ==========================================================================
    # 17. Process year/locality-level information
    # ==========================================================================
    
    year_locality <- year_locality_raw %>%
      
      # Keep the original locality field but give it a clearer name
      select(
        locality_vernacular_name = locality,
        everything()
      ) %>%
      
      select(
        -dynamicProperties
      )
    
    
    # --------------------------------------------------------------------------
    # Unpack dynamic year/locality properties
    # --------------------------------------------------------------------------
    
    tempDynamic <- year_locality_raw %>%
      select(
        dyn = dynamicProperties
      ) %>%
      unlist() %>%
      tidyjson::spread_all() %>%
      as_tibble() %>%
      select(
        -document.id
      )
    
    
    # --------------------------------------------------------------------------
    # Join dynamic properties
    # --------------------------------------------------------------------------
    
    year_locality <- year_locality %>%
      left_join(
        tempDynamic,
        by = c(
          "id" = "id"
        )
      )
    
    
    # --------------------------------------------------------------------------
    # Validate the join
    # --------------------------------------------------------------------------
    
    if (
      !nrow(year_locality) == nrow(year_locality_raw) &
      nrow(year_locality) == nrow(tempDynamic)
    ) {
      stop(
        "Year locality rows don't match!"
      )
    }
    
    
    # ==========================================================================
    # 18. Convert the year/locality sampling period
    # ==========================================================================
    
    year_locality <- year_locality %>%
      
      separate(
        eventTime,
        c(
          "start_time",
          "end_time"
        ),
        sep = "/"
      ) %>%
      
      mutate(
        season_start_time = as.POSIXct(start_time),
        season_end_time = as.POSIXct(end_time)
      ) %>%
      
      mutate(
        ssbid = as.character(
          ssbid
        )
      )
    
    
    # Optional inspection of the season sampling period
    year_locality %>%
      select(
        start_time,
        end_time
      )
    
    
    # ==========================================================================
    # 19. Retain relevant year/locality variables
    # ==========================================================================
    
    year_locality <- year_locality %>%
      select(
        year_locality_id = eventID,
        season_sample_size_value = sampleSizeValue,
        season_sample_size_unit = sampleSizeUnit,
        season_start_time,
        season_end_time,
        locality,
        year,
        habitat_type,
        region_name,
        ssbid,
        ano_flate_id,
        no_herb_spec,
        avg_prc_cov_herb_species,
        no_tree_spec,
        dom_tree_spec,
        avg_dom_tree_age,
        centroid_latitude = decimalLatitude,
        centroid_longitude = decimalLongitude,
        polygon_wkt = footprintWKT,
        coordinate_system = geodeticDatum
      )
    
    
    # ==========================================================================
    # 20. Process insect occurrence information
    # ==========================================================================
    #
    # The occurrence table contains the taxonomic identity of each recorded
    # insect occurrence.
    # ==========================================================================
    
    occurrence <- occurrence_raw %>%
      select(
        occurrence_id = occurrenceID,
        identification_id = eventID,
        quantity = organismQuantity,
        quantity_type = organismQuantityType,
        scientific_name = scientificName,
        vernacular_name = vernacularName,
        scientific_name_id = taxonID,
        class,
        order,
        family,
        genus,
        specific_epithet = specificEpithet
      )
    
    
    # ==========================================================================
    # 21. Save intermediate tables
    # ==========================================================================
    #
    # These files provide useful checkpoints and make it possible to inspect
    # each level of the processed monitoring hierarchy independently.
    # ==========================================================================
    
    system(
      "mkdir -p 'out'"
    )
    
    
    write_csv(
      occurrence,
      file = "out/occurrence.csv"
    )
    
    write_csv(
      identifications,
      file = "out/identifications.csv"
    )
    
    write_csv(
      sampling_trap,
      file = "out/sampling_trap.csv"
    )
    
    write_csv(
      locality_sampling,
      file = "out/locality_sampling.csv"
    )
    
    write_csv(
      year_locality,
      file = "out/year_locality.csv"
    )
    
    
    # ==========================================================================
    # 22. Join the complete monitoring hierarchy
    # ==========================================================================
    #
    # The tables are joined from the lowest level (occurrence) upward:
    #
    #     occurrence
    #         ↓
    #     identification
    #         ↓
    #     sampling trap
    #         ↓
    #     locality sampling
    #         ↓
    #     year/locality
    #
    # This produces one integrated table containing both taxonomic observations
    # and information about where, when, and under what conditions the insects
    # were sampled.
    # ==========================================================================
    
    spec_join <- occurrence %>%
      
      left_join(
        identifications,
        by = c(
          "identification_id" =
            "identification_id"
        )
      ) %>%
      
      left_join(
        sampling_trap,
        by = c(
          "sampling_trap_id" =
            "sampling_trap_id"
        )
      ) %>%
      
      left_join(
        locality_sampling,
        by = c(
          "locality_sampling_id" =
            "locality_sampling_id"
        )
      ) %>%
      
      left_join(
        year_locality,
        by = c(
          "year_locality_id" =
            "year_locality_id"
        )
      ) %>%
      
      # ------------------------------------------------------------------------
    # Calculate Julian day from the end of the sampling event
    # ------------------------------------------------------------------------
    
    mutate(
      julian_day = lubridate::yday(
        end_time
      ),
      
      # Convert year to a factor
      year = forcats::as_factor(
        year
      )
    )
    
    
    # ==========================================================================
    # 23. Save the complete processed monitoring dataset
    # ==========================================================================
    
    write_csv(
      spec_join,
      file = "out/nationalInsectMonitoring.csv"
    )
    
    
    # Store the processed data in the main object
    nationalInsectMonitoring <- spec_join
    
    
  } else {
    
    
    # ==========================================================================
    # 24. If the processed CSV already exists, load it instead
    # ==========================================================================
    
    nationalInsectMonitoring <- readr::read_csv(
      "out/nationalInsectMonitoring.csv",
      show_col_types = FALSE
    )
  }
  
  
  # ============================================================================
  # 25. Convert insect monitoring records into spatial points
  # ============================================================================
  #
  # Trap coordinates are originally stored as longitude/latitude in WGS84.
  #
  # They are converted to sf points and then transformed to the analysis CRS.
  # ============================================================================
  
  insectData <- nationalInsectMonitoring %>%
    
    sf::st_as_sf(
      coords = c(
        "point_longitude",
        "point_latitude"
      ),
      crs = 4326
    ) %>%
    
    sf::st_transform(
      crs
    ) %>%
    
    # --------------------------------------------------------------------------
  # Standardise the scientific name
  #
  # Only the first two words are retained:
  #
  #     Genus species
  # --------------------------------------------------------------------------
  
  mutate(
    simpleScientificName = stringr::str_extract(
      scientific_name,
      "^[A-Za-z]+\\s+[a-z]+"
    )
  )
  
  
  # ============================================================================
  # 26. Identify insect species represented in the insect SDM
  # ============================================================================
  #
  # `species_names` contains the insect species for which prediction rasters
  # have been generated.
  #
  # Only observations belonging to these species are retained.
  # ============================================================================
  
  nina_insects <- insectData %>%
    
    mutate(
      insect_name = gsub(
        " ",
        "_",
        simpleScientificName
      )
    ) %>%
    
    filter(
      insect_name %in% species_names
    )
  
  
  # ============================================================================
  # 27. Restrict the observations to semi-natural habitats
  # ============================================================================
  #
  # Only monitoring records with:
  #
  #     habitat_type == "Semi-nat"
  #
  # are retained for the richness calculation.
  #
  # This makes the observed richness measure specific to the semi-natural
  # habitat of interest.
  # ============================================================================
  
  nina_seminat <- nina_insects %>%
    
    dplyr::filter(
      habitat_type == "Semi-nat"
    )
  
  
  # ============================================================================
  # 28. Calculate observed insect species richness
  # ============================================================================
  #
  # Richness is calculated at:
  #
  #     locality × year
  #
  # level.
  #
  # For each locality and year:
  #
  #     insect_richness =
  #         number of distinct insect species recorded
  #
  # Additional sampling variables are summarised at the same level.
  # ============================================================================
  
  richness <- nina_seminat %>%
    
    # Remove geometry because the richness calculation is tabular
    sf::st_drop_geometry() %>%
    
    dplyr::group_by(
      locality,
      year
    ) %>%
    
    dplyr::summarise(
      
      # ------------------------------------------------------------------------
      # Number of distinct insect species recorded
      # ------------------------------------------------------------------------
      
      insect_richness =
        n_distinct(
          scientific_name
        ),
      
      
      # ------------------------------------------------------------------------
      # Mean minimum temperature during sampling
      # ------------------------------------------------------------------------
      
      mean_sampling_min_temp =
        mean(
          sampling_min_temp,
          na.rm = TRUE
        ),
      
      
      # ------------------------------------------------------------------------
      # Mean maximum temperature during sampling
      # ------------------------------------------------------------------------
      
      mean_sampling_max_temp =
        mean(
          sampling_max_temp,
          na.rm = TRUE
        ),
      
      
      # ------------------------------------------------------------------------
      # Mean trap net wet weight
      #
      # This is retained as a measure of sampling/collection characteristics.
      # ------------------------------------------------------------------------
      
      mean_weight =
        mean(
          net_wet_weight,
          na.rm = TRUE
        ),
      
      
      # ------------------------------------------------------------------------
      # Total sampling duration
      #
      # The duration is calculated from the earliest start time to the latest
      # end time within each locality × year group.
      # ------------------------------------------------------------------------
      
      sampling_duration_hours =
        as.numeric(
          difftime(
            max(
              end_time,
              na.rm = TRUE
            ),
            min(
              start_time,
              na.rm = TRUE
            ),
            units = "hours"
          )
        ),
      
      .groups = "drop"
    )
  
  
  # ============================================================================
  # 29. Create one spatial record per locality
  # ============================================================================
  #
  # The richness calculation is at locality × year level, but the spatial
  # locations are represented by one geometry per locality.
  #
  # `distinct()` therefore retains ONE spatial feature for each locality.
  # ============================================================================
  
  locations <- nina_seminat %>%
    
    dplyr::select(
      locality,
      geometry
    ) %>%
    
    dplyr::distinct(
      locality,
      .keep_all = TRUE
    )
  
  
  # ============================================================================
  # 30. Join richness and sampling information to spatial locations
  # ============================================================================
  #
  # The calculated richness table is joined to the spatial locality points.
  #
  # The resulting object is:
  #
  #     national_monitoring_insect_richness
  #
  # and contains observed insect richness together with the associated
  # sampling variables.
  # ============================================================================
  
  national_monitoring_insect_richness <- locations %>%
    
    dplyr::left_join(
      richness,
      by = "locality"
    ) %>%
    
    dplyr::mutate(
      dataset = "national_insect_monitoring"
    )
  
  
  # ============================================================================
  # 31. Plot observed insect richness
  # ============================================================================
  #
  # This provides a quick spatial diagnostic showing where the monitoring
  # locations are and the observed richness associated with them.
  # ============================================================================
  
  library(ggplot2)
  
  
  ggplot(
    national_monitoring_insect_richness
  ) +
    
    geom_sf(
      aes(
        colour = insect_richness
      ),
      size = 1.5,
      alpha = 0.7
    ) +
    
    scale_colour_viridis_c(
      name = "Insect richness"
    ) +
    
    labs(
      title = "Observed insect richness",
      subtitle = "National insect monitoring"
    ) +
    
    theme_bw()
  
  
  # ============================================================================
  # 32. Prepare final output for Shapefile
  # ============================================================================
  #
  # Several variable names are shortened here to make them easier to work
  # with in subsequent modelling scripts.
  #
  # The resulting names are:
  #
  #     local    = locality
  #     insct_r  = insect_richness
  #     min_temp = mean_sampling_min_temp
  #     max_temp = mean_sampling_max_temp
  #     weight   = mean_weight
  #     samp_dur = sampling_duration_hours
  # ============================================================================
  
  national_monitoring_insect_richness_shp <-
    national_monitoring_insect_richness %>%
    
    dplyr::rename(
      local = locality,
      insct_r = insect_richness,
      min_temp = mean_sampling_min_temp,
      max_temp = mean_sampling_max_temp,
      weight = mean_weight,
      samp_dur = sampling_duration_hours
    )
  
  
  # ============================================================================
  # 33. Save the final dataset as a Shapefile
  # ============================================================================
  
  sf::st_write(
    national_monitoring_insect_richness_shp,
    national_monitoring_file,
    delete_layer = TRUE,
    quiet = FALSE
  )
  
  
  # ============================================================================
  # 34. Reload the saved Shapefile
  # ============================================================================
  #
  # The final object is reloaded from disk so that downstream analyses use
  # exactly the data that were written to the output file.
  # ============================================================================
  
  national_monitoring_insect_richness <- sf::st_read(
    national_monitoring_file,
    quiet = TRUE
  )
  
  
} else {
  
  
  # ============================================================================
  # 35. If the final Shapefile already exists, load it directly
  # ============================================================================
  
  national_monitoring_insect_richness <- sf::st_read(
    national_monitoring_file,
    quiet = TRUE
  )
}

