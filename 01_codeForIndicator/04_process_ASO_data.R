################################################################################
# Calculate observed plant richness at ASO locations
#
# This script:
#
#   1. Loads ASO plant occurrence records.
#
#   2. Converts longitude/latitude coordinates into spatial points.
#
#   3. Reprojects the observations to the analysis CRS.
#
#   4. Extracts standardised plant species names.
#
#   5. Retains only plant species represented in the plant SDM.
#
#   6. Identifies unique ASO sampling locations.
#
#   7. Creates a plant × location detection matrix:
#
#          1 = plant species observed at least once at the location
#          0 = plant species not observed at the location
#
#   8. Calculates observed plant species richness at each ASO location.
#
#   9. Joins the richness values back to the spatial locations.
#
#  10. Saves the resulting locations as a Shapefile.
#
# Final output:
#
#     aso_plant_richness
#
# Each feature represents one unique ASO location and contains:
#
#     - x coordinate
#     - y coordinate
#     - geometry
#     - plant_richness
#     - dataset = "ASO"
#
################################################################################


# ==============================================================================
# 0. Load required packages
# ==============================================================================

library(dplyr)
library(tidyr)
library(sf)


# ==============================================================================
# 1. Define the output file
# ==============================================================================

aso_richness_file <- file.path(
  resultFolder,
  "polygons",
  "aso_location_with_richness.shp"
)


# ==============================================================================
# 2. Create the ASO richness dataset if it does not already exist
# ==============================================================================

if (!file.exists(aso_richness_file)) {
  
  
  # ============================================================================
  # 3. Read and prepare ASO occurrence data
  # ============================================================================
  #
  # The ASO occurrence data are stored as a tab-delimited text file.
  #
  # The original coordinates are:
  #
  #     decimalLongitude
  #     decimalLatitude
  #
  # and are assumed to use WGS84 (EPSG:4326).
  #
  # The observations are then transformed to the analysis CRS used throughout
  # the spatial modelling workflow.
  # ============================================================================
  
  aso <- readr::read_delim(
    file.path(
      resultFolder,
      "modelFit",
      "asoOccurrences.csv"
    ),
    delim = "\t",
    escape_double = FALSE,
    trim_ws = TRUE
  ) %>%
    
    # --------------------------------------------------------------------------
  # Convert longitude/latitude to spatial point features
  # --------------------------------------------------------------------------
  
  sf::st_as_sf(
    coords = c(
      "decimalLongitude",
      "decimalLatitude"
    ),
    crs = 4326
  ) %>%
    
    # --------------------------------------------------------------------------
  # Transform the points to the common analysis CRS
  # --------------------------------------------------------------------------
  
  sf::st_transform(
    crs
  ) %>%
    
    # --------------------------------------------------------------------------
  # Extract a standardised scientific name
  #
  # Only the first two words of the scientific name are retained:
  #
  #     Genus species
  #
  # For example:
  #
  #     "Achillea millefolium L."
  #
  # becomes:
  #
  #     "Achillea millefolium"
  #
  # This is necessary to match the species names used by the plant SDM.
  # --------------------------------------------------------------------------
  
  mutate(
    simpleScientificName = stringr::str_extract(
      scientificName,
      "^[A-Za-z]+\\s+[a-z]+"
    )
  )
  
  
  # =============================================================================
  # 4. Identify plant species represented in the plant SDM
  # =============================================================================
  #
  # The names of the layers in `plants` define the plant species included in
  # the SDM.
  #
  # The observed ASO richness calculation will therefore use exactly the same
  # species set as the plant predictions.
  # =============================================================================
  
  plant_names <- names(
    plants
  )
  
  
  # =============================================================================
  # 5. Keep only ASO records for species represented in the plant SDM
  # =============================================================================
  #
  # ASO may contain observations for plant species that are not included in
  # the current SDM.
  #
  # These records are removed so that observed richness is calculated only
  # for species represented by the model.
  #
  # Spaces are replaced with underscores to match the raster layer naming
  # convention.
  #
  # Example:
  #
  #     Achillea millefolium
  #
  # becomes:
  #
  #     Achillea_millefolium
  # =============================================================================
  
  aso_plants <- aso %>%
    
    mutate(
      plant_name = gsub(
        " ",
        "_",
        simpleScientificName
      )
    ) %>%
    
    filter(
      plant_name %in% plant_names
    )
  
  
  # =============================================================================
  # 6. Identify unique ASO sampling locations
  # =============================================================================
  #
  # Multiple plant observations can occur at the same ASO location.
  #
  # We therefore identify unique locations using their projected x and y
  # coordinates.
  #
  # Only one spatial feature is retained for each unique coordinate pair.
  # =============================================================================
  
  aso_locations <- aso_plants %>%
    
    mutate(
      x = sf::st_coordinates(geometry)[, 1],
      y = sf::st_coordinates(geometry)[, 2]
    ) %>%
    
    distinct(
      x,
      y,
      .keep_all = TRUE
    ) %>%
    
    select(
      x,
      y,
      geometry
    )
  
  
  # =============================================================================
  # 7. Create the plant detection matrix
  # =============================================================================
  #
  # The detection matrix has:
  #
  #     rows    = ASO sampling locations
  #     columns = plant species
  #
  # Each cell contains:
  #
  #     1 = species observed at the location
  #     0 = species not observed at the location
  #
  # If a species is recorded multiple times at the same location, it still
  # contributes only ONE detection.
  # =============================================================================
  
  aso_presence <- aso_plants %>%
    
    # --------------------------------------------------------------------------
  # Extract projected coordinates
  # --------------------------------------------------------------------------
  
  mutate(
    x = sf::st_coordinates(geometry)[, 1],
    y = sf::st_coordinates(geometry)[, 2]
  ) %>%
    
    # --------------------------------------------------------------------------
  # Retain one record per location × species combination
  # --------------------------------------------------------------------------
  
  distinct(
    x,
    y,
    plant_name
  ) %>%
    
    # --------------------------------------------------------------------------
  # Mark the species as detected
  # --------------------------------------------------------------------------
  
  mutate(
    detected = 1L
  ) %>%
    
    # --------------------------------------------------------------------------
  # Convert from long format to a location × species matrix
  #
  # Species that were not detected at a location receive zero.
  # --------------------------------------------------------------------------
  
  pivot_wider(
    id_cols = c(x, y),
    names_from = plant_name,
    values_from = detected,
    values_fill = 0
  )
  
  
  # =============================================================================
  # 8. Add plant species that were never observed in ASO
  # =============================================================================
  #
  # A plant species may be present in the plant SDM but absent from all ASO
  # records.
  #
  # Such species will not automatically appear as columns in `aso_presence`.
  #
  # They are explicitly added with zero detection at every ASO location.
  #
  # This ensures that the detection matrix contains the complete set of
  # plant species represented in the SDM.
  # =============================================================================
  
  missing_plants <- setdiff(
    plant_names,
    names(aso_presence)
  )
  
  
  if (length(missing_plants) > 0) {
    
    aso_presence[missing_plants] <- 0L
  }
  
  
  # =============================================================================
  # 9. Ensure plant columns follow the same order as the plant raster layers
  # =============================================================================
  #
  # This ensures that:
  #
  #     ASO species columns
  #
  # are ordered in exactly the same way as:
  #
  #     plant raster layers
  #
  # Coordinate columns remain at the beginning.
  # =============================================================================
  
  aso_presence <- aso_presence %>%
    
    select(
      x,
      y,
      all_of(plant_names)
    )
  
  
  # =============================================================================
  # 10. Calculate observed plant species richness
  # =============================================================================
  #
  # For each ASO location:
  #
  #     plant richness =
  #         sum of detected plant species
  #
  # Because every species is represented by either 0 or 1, this is simply the
  # number of modelled plant species observed at that location.
  # =============================================================================
  
  aso_presence <- aso_presence %>%
    
    mutate(
      plant_richness = rowSums(
        across(
          all_of(plant_names)
        )
      )
    )
  
  
  # =============================================================================
  # 11. Add richness values to the spatial ASO locations
  # =============================================================================
  #
  # The calculated richness is joined back to the spatial point dataset using
  # the projected x and y coordinates.
  #
  # The resulting object contains one point per unique ASO location.
  # =============================================================================
  
  aso_plant_richness <- aso_locations %>%
    
    select(
      x,
      y,
      geometry
    ) %>%
    
    # --------------------------------------------------------------------------
  # Add observed plant richness
  # --------------------------------------------------------------------------
  
  left_join(
    aso_presence %>%
      select(
        x,
        y,
        plant_richness
      ),
    by = c(
      "x",
      "y"
    )
  ) %>%
    
    # --------------------------------------------------------------------------
  # Record the source dataset
  # --------------------------------------------------------------------------
  
  mutate(
    dataset = "ASO"
  )
  
  
  # =============================================================================
  # 12. Save the ASO richness dataset
  # =============================================================================
  
  sf::st_write(
    aso_plant_richness,
    aso_richness_file,
    delete_layer = TRUE,
    quiet = FALSE
  )
  
  
  # =============================================================================
  # 13. Reload the saved Shapefile
  # =============================================================================
  #
  # This ensures that the object used downstream corresponds to the spatial
  # dataset actually written to disk.
  # =============================================================================
  
  aso_plant_richness <- sf::st_read(
    aso_richness_file,
    quiet = TRUE
  )
  
  
} else {
  
  
  # =============================================================================
  # 14. If the output already exists, load it directly
  # =============================================================================
  
  aso_plant_richness <- sf::st_read(
    aso_richness_file,
    quiet = TRUE
  )
}
