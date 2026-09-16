library(dplyr)
library(tidyr)
library(sf)

if(!file.exists(file.path(
  resultFolder, "polygons",
  paste0("ninaInsect_with_richness_for_",var_of_interest,".shp")
))){
  ## ========================================================
  ## 1. Read and prepare ASO data
  ## ========================================================
  
  insectData <- readr::read_delim(
    file.path(
      resultFolder,
      "modelFit/ninaInsectData.csv"
    ),
    delim = "\t",
    escape_double = FALSE,
    trim_ws = TRUE
  ) %>%
    sf::st_as_sf(
      coords = c(
        "decimalLongitude",
        "decimalLatitude"
      ),
      crs = 4326
    ) %>%
    st_transform(
      crs
    ) %>%
    mutate(
      simpleScientificName = coalesce(
        str_extract(
          scientificName,
          "^[A-Za-z]+\\s+[a-z]+"
        )
      )
    )
  
  
  ## ========================================================
  ## 2. Insect species represented in the SDM
  ## ========================================================
  
  species_names
  
  
  ## ========================================================
  ## 3. Keep only ASO records for species in plants
  ## ========================================================
  
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
  
  
  ## ========================================================
  ## 4. Extract unique ASO locations
  ## ========================================================
  
  nina_insects_locations <- nina_insects %>%
    mutate(
      x = st_coordinates(geometry)[, 1],
      y = st_coordinates(geometry)[, 2]
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
  
  
  ## ========================================================
  ## 5. Create detection matrix
  ##
  ## 1 = plant observed at least once at location
  ## 0 = plant not observed at location
  ## ========================================================
  
  insect_presence <- nina_insects %>%
    mutate(
      x = st_coordinates(geometry)[, 1],
      y = st_coordinates(geometry)[, 2]
    ) %>%
    distinct(
      x,
      y,
      insect_name
    ) %>%
    mutate(
      detected = 1L
    ) %>%
    pivot_wider(
      id_cols = c(x, y),
      names_from = insect_name,
      values_from = detected,
      values_fill = 0
    )
  
  
  ## ========================================================
  ## 6. Add plant species that were never observed
  ## ========================================================
  
  missing_plants <- setdiff(
    species_names,
    names(insect_presence)
  )
  
  if (length(missing_plants) > 0) {
    insect_presence[missing_plants] <- 0L
  }
  
  
  ## ========================================================
  ## 7. Ensure plant columns follow raster order
  ## ========================================================
  
  insect_presence <- insect_presence %>%
    select(
      x,
      y,
      all_of(species_names)
    )
  
  
  ## ========================================================
  ## 8. Calculate observed plant species richness
  ## ========================================================
  
  insect_presence <- insect_presence %>%
    mutate(
      insect_richness = rowSums(
        across(
          all_of(species_names)
        )
      )
    )
  
  
  ## ========================================================
  ## 9. Add richness to ASO locations
  ## ========================================================
  
  nina_insect_richness <- nina_insects_locations %>%
    select(
      x,
      y,
      geometry
    ) %>%
    left_join(
      insect_presence %>%
        select(
          x,
          y,
          insect_richness
        ),
      by = c("x", "y")
    ) %>%
    dplyr::mutate(dataset = "NINA_insect")
  
  
  ## ========================================================
  ## Result
  ## ========================================================
  
  shapefile_path <- file.path(
    resultFolder, "polygons",
    paste0("ninaInsect_with_richness_for_",var_of_interest,".shp")
  )
  
  sf::st_write(
    nina_insect_richness,
    shapefile_path,
    delete_layer = TRUE,
    quiet = FALSE
  )
  
  nina_insect_richness <- sf::st_read(
    file.path(
      resultFolder, "polygons",
      paste0("ninaInsect_with_richness_for_",var_of_interest,".shp")
    ),
    quiet = TRUE
  )
} else {
  nina_insect_richness <- sf::st_read(
    file.path(
      resultFolder, "polygons",
      paste0("ninaInsect_with_richness_for_",var_of_interest,".shp")
    ),
    quiet = TRUE
  )
}