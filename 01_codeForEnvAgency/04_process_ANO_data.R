################################################################################
# Calculate observed plant richness at ANO locations
#
# This script:
#
#   1. Loads ANO occurrence records.
#
#   2. Converts the longitude/latitude coordinates into spatial points.
#
#   3. Reprojects the observations to the analysis CRS.
#
#   4. Extracts plant species names and retains only species represented in
#      the plant SDM predictions.
#
#   5. Identifies unique ANO sampling locations.
#
#   6. Creates a plant × location detection matrix:
#
#          1 = plant species observed at least once at the location
#          0 = plant species not observed at the location
#
#   7. Calculates observed plant species richness at each ANO location.
#
#   8. Saves the resulting spatial dataset as a Shapefile.
#
# Final output:
#
#     ano_plant_richness
#
# Each feature represents one unique ANO location and contains:
#
#     - x coordinate
#     - y coordinate
#     - geometry
#     - plant_richness
#     - dataset = "ANO"
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

ano_richness_file <- file.path(
  resultFolder,
  "polygons",
  "ano_location_with_richness.shp"
)


# ==============================================================================
# 2. Create the ANO richness dataset if it does not already exist
# ==============================================================================

if (!file.exists(ano_richness_file)) {
  
  
  # ============================================================================
  # 3. Read and prepare ANO occurrence data
  # ============================================================================
  #
  # The ANO occurrence data are stored as a tab-delimited text file.
  #
  # The longitude and latitude fields are initially in geographic coordinates
  # using WGS84 (EPSG:4326).
  #
  # They are then transformed to the analysis CRS used by the rest of the
  # workflow.
  # ============================================================================
  
  ano <- readr::read_delim(
    file.path(
      resultFolder,
      "modelFit",
      "anoOccurrences.csv"
    ),
    delim = "\t",
    escape_double = FALSE,
    trim_ws = TRUE
  ) %>%
    
    # --------------------------------------------------------------------------
  # Convert occurrence records into spatial point features
  # --------------------------------------------------------------------------
  
  sf::st_as_sf(
    coords = c(
      "decimalLongitude",
      "decimalLatitude"
    ),
    crs = 4326
  ) %>%
    
    # --------------------------------------------------------------------------
  # Transform the points to the analysis CRS
  #
  # `crs` is assumed to contain the coordinate reference system used by
  # the wider analysis.
  # --------------------------------------------------------------------------
  
  sf::st_transform(
    crs
  ) %>%
    
    # --------------------------------------------------------------------------
  # Extract a standardised scientific name
  #
  # Some occurrence records may contain additional taxonomic information
  # after the species name.
  #
  # The regular expression:
  #
  #     ^[A-Za-z]+\\s+[a-z]+
  #
  # extracts the first two words:
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
  # The names of the layers in `plants` define the species that are included
  # in the plant SDM.
  #
  # Only ANO observations belonging to these species will be used to calculate
  # richness.
  # =============================================================================
  
  plant_names <- names(
    plants
  )
  
  
  # =============================================================================
  # 5. Keep only ANO records for species represented in the plant SDM
  # =============================================================================
  #
  # ANO may contain observations for species that are not included in the
  # current plant SDM.
  #
  # Those observations are excluded so that observed richness is calculated
  # using exactly the same species set as the plant predictions.
  #
  # Spaces in species names are replaced with underscores to match the naming
  # convention used by the raster layers.
  #
  # Example:
  #
  #     "Achillea millefolium"
  #
  # becomes:
  #
  #     "Achillea_millefolium"
  # =============================================================================
  
  ano_plants <- ano %>%
    
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
  # 6. Identify unique ANO sampling locations
  # =============================================================================
  #
  # The occurrence data may contain multiple observations of the same species
  # and/or multiple species at the same coordinates.
  #
  # We therefore identify unique spatial locations based on their projected
  # x and y coordinates.
  #
  # Only one feature is retained for each unique coordinate pair.
  # =============================================================================
  
  ano_locations <- ano_plants %>%
    
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
  # The goal is to create a matrix where:
  #
  #     rows    = ANO sampling locations
  #     columns = plant species
  #
  # Each cell contains:
  #
  #     1 = species was observed at that location
  #     0 = species was not observed at that location
  #
  # Multiple observations of the same species at the same location are
  # collapsed into a single detection.
  # =============================================================================
  
  ano_presence <- ano_plants %>%
    
    # --------------------------------------------------------------------------
  # Extract the projected coordinates
  # --------------------------------------------------------------------------
  
  mutate(
    x = sf::st_coordinates(geometry)[, 1],
    y = sf::st_coordinates(geometry)[, 2]
  ) %>%
    
    # --------------------------------------------------------------------------
  # Keep one record per location × species combination
  #
  # This prevents repeated observations of the same species from affecting
  # the richness calculation.
  # --------------------------------------------------------------------------
  
  distinct(
    x,
    y,
    plant_name
  ) %>%
    
    # --------------------------------------------------------------------------
  # Mark every observed species as detected
  # --------------------------------------------------------------------------
  
  mutate(
    detected = 1L
  ) %>%
    
    # --------------------------------------------------------------------------
  # Convert the long-format occurrence data into a wide detection matrix
  #
  # Missing combinations of location × species are assigned zero.
  # --------------------------------------------------------------------------
  
  pivot_wider(
    id_cols = c(x, y),
    names_from = plant_name,
    values_from = detected,
    values_fill = 0
  )
  
  
  # =============================================================================
  # 8. Add plant species that were never observed in the ANO data
  # =============================================================================
  #
  # A plant species may be included in the plant SDM but have no observations
  # in the ANO dataset.
  #
  # Such species will not appear as columns in `ano_presence`.
  #
  # We explicitly add them with zero detection at every location so that the
  # richness calculation uses the complete plant SDM species set.
  # =============================================================================
  
  missing_plants <- setdiff(
    plant_names,
    names(ano_presence)
  )
  
  
  if (length(missing_plants) > 0) {
    
    ano_presence[missing_plants] <- 0L
  }
  
  
  # =============================================================================
  # 9. Ensure plant columns are in the same order as the plant SDM
  # =============================================================================
  #
  # This makes the detection matrix consistent with the ordering of the
  # plant prediction raster layers.
  #
  # The coordinate columns remain first, followed by plant species.
  # =============================================================================
  
  ano_presence <- ano_presence %>%
    
    select(
      x,
      y,
      all_of(plant_names)
    )
  
  
  # =============================================================================
  # 10. Calculate observed plant species richness
  # =============================================================================
  #
  # For each ANO location:
  #
  #     plant richness =
  #         sum of detected plant species
  #
  # Because each plant species is represented by either 0 or 1, this is simply
  # the number of plant species observed at the location.
  # =============================================================================
  
  ano_presence <- ano_presence %>%
    
    mutate(
      plant_richness = rowSums(
        across(
          all_of(plant_names)
        )
      )
    )
  
  
  # =============================================================================
  # 11. Add richness values to the spatial ANO locations
  # =============================================================================
  #
  # The richness values are joined back onto the spatial point dataset using
  # the projected x and y coordinates.
  #
  # The resulting object contains one spatial point per unique ANO location.
  # ==============================================================================
  
  ano_plant_richness <- ano_locations %>%
    
    select(
      x,
      y,
      geometry
    ) %>%
    
    # --------------------------------------------------------------------------
  # Add observed plant richness
  # --------------------------------------------------------------------------
  
  left_join(
    ano_presence %>%
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
  # Identify the source dataset
  # --------------------------------------------------------------------------
  
  mutate(
    dataset = "ANO"
  )
  
  
  # =============================================================================
  # 12. Save the ANO plant-richness locations
  # =============================================================================
  
  sf::st_write(
    ano_plant_richness,
    ano_richness_file,
    delete_layer = TRUE,
    quiet = FALSE
  )
  
  
  # =============================================================================
  # 13. Reload the saved Shapefile
  # =============================================================================
  #
  # Reloading ensures that the object used downstream corresponds to the
  # dataset that has actually been written to disk.
  # =============================================================================
  
  ano_plant_richness <- sf::st_read(
    ano_richness_file,
    quiet = TRUE
  )
  
  
} else {
  
  
  # =============================================================================
  # 14. If the output already exists, load it directly
  # =============================================================================
  
  ano_plant_richness <- sf::st_read(
    ano_richness_file,
    quiet = TRUE
  )
}
