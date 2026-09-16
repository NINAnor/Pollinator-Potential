# ==============================================================================
# POLLINATOR INDICATOR ANALYSIS PIPELINE
# ==============================================================================
#
# This script runs the complete pollinator indicator analysis.
#
# The analysis is organised as a pipeline and calls supporting R scripts from
# the functions/ and pipeline/ directories.
#
# In most cases, the required data have already been downloaded and formatted.
# This script does not automatically re-download or re-format the data.
#
# To force the data to be downloaded and formatted again, manually remove the
# corresponding formatted data. The relevant pipeline scripts will then
# recreate them.
#
# ==============================================================================


# ==============================================================================
# 1. SET WORKING ENVIRONMENT
# ==============================================================================

## ------------------------------------------------------------------
## Choose execution environment
## ------------------------------------------------------------------

ninaServer <- FALSE

if (ninaServer) {
  
  setwd(
    "/data/P-Prosjekter2/412421_okologisk_tilstand_2024/pollinatorIndicators/oldScripts"
  )
  
} else {
  
  setwd(
    "P:/412421_okologisk_tilstand_2024/pollinatorIndicators/oldScripts"
  )
}


# ==============================================================================
# 2. LOAD FUNCTIONS
# ==============================================================================


options(encoding = "UTF-8")
Sys.setenv(LANG = "en_US.UTF-8")

## Source all functions used by the analysis
invisible(
  lapply(
    list.files(
      "functions",
      pattern = "\\.R$",
      full.names = TRUE
    ),
    source
  )
)


# ==============================================================================
# 3. LOAD PACKAGES AND SET TERRA OPTIONS
# ==============================================================================

library(terra)
library(lme4)
library(emmeans)
library(ggplot2)
library(dplyr)

terraOptions(
  memfrac  = 0.5,       # Fraction of available RAM
  memmax   = 4,          # Maximum memory allocation in GB
  threads  = 2,          # Number of processing threads
  parallel = TRUE,
  todisk   = TRUE,       # Prefer disk-based processing
  tempdir  = "C:/terra_tmp",
  datatype = "FLT4S",
  progress = 3
)

## --------------------------------------------------------
## Encoding
## --------------------------------------------------------

options(encoding = "UTF-8")

# ==============================================================================
# 4. ANALYSIS PARAMETERS
# ==============================================================================

## ------------------------------------------------------------------
## Input and output directories
## ------------------------------------------------------------------

dataFolder <- "Data"

resultFolder <- "C:/ninaProjects/pollinatorIndicator/kwakuAugust2026/kwakuAugust2026"


## ------------------------------------------------------------------
## Model and plotting options
## ------------------------------------------------------------------

modelPlots <- FALSE
prediction_year = 2024

## ------------------------------------------------------------------
## Open-lowland threshold
## ------------------------------------------------------------------

openlowland_threshold <- 1


## ------------------------------------------------------------------
## Coordinate reference systems
## ------------------------------------------------------------------

## EPSG code used for the spatial analysis
crs <- 25833

## CRS used for defining the study region
newCrs <-
  "+proj=utm +zone=33 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=km +no_defs"


# ==============================================================================
# 5. DEFINE STUDY REGION
# ==============================================================================

regionGeometry <- defineRegion(
  level  = "country",
  region = "Norway",
  crs    = newCrs
)


# ==============================================================================
# 6. LOAD SPECIES AND TRAIT INFORMATION
# ==============================================================================

## ------------------------------------------------------------------
## Load functions used for interaction/model plots
## ------------------------------------------------------------------

source(
  "pipeline/webForInterractions/webPlotsForModels.R"
)


## ------------------------------------------------------------------
## Load biotraits if not already available
## ------------------------------------------------------------------

if (!exists("biotraits")) {
  
  source(
    "pipeline/predictions/pollinatorZoneGroupings.R"
  )
}


## ------------------------------------------------------------------
## Create species classification table
## ------------------------------------------------------------------

speciesForTaxon <- biotraits %>%
  
  dplyr::select(
    acceptedScientificName,
    Taxon,
    zone
  ) %>%
  
  dplyr::group_by(
    Taxon,
    zone
  ) %>%
  
  dplyr::distinct() %>%
  
  dplyr::ungroup() %>%
  
  dplyr::mutate(
    
    simpleScientificName = stringr::str_extract(
      acceptedScientificName,
      "^[A-Za-z]+\\s+[a-z]+"
    )
    
  ) %>%
  
  dplyr::filter(
    !is.na(simpleScientificName)
  )


## ------------------------------------------------------------------
## Define groups for which indicators will be calculated
## ------------------------------------------------------------------

vars_of_interest <- c(
  "all",
  unique(speciesForTaxon$Taxon),
  unique(speciesForTaxon$zone)

)


# ========================================
# Get caliberation raster
# =======================================

source("01_codeForEnvAgency/00D_get_calibration_raster.R")

# ==============================================================================
# 8. CREATE OPEN-LOWLAND MASK
# ==============================================================================
useCLCForLowlandBoundary <- FALSE
useCLCplusForLowlandBoundary <- TRUE
source("01_codeForEnvAgency/00E_get_openlowland_boundary.R")


# ==============================================================================
# 9. LOAD REFERENCE LOCATIONS
# ==============================================================================
source("01_codeForEnvAgency/00F_get_reference_locations.R")

## -----------------------
## ASO and ANO data locations
## -----------------------

source("01_codeForEnvAgency/04_process_ASO_data.R")
source("01_codeForEnvAgency/04_process_ANO_data.R")

evaluation_locs <- bind_rows(aso_plant_richness,
                             ano_plant_richness)

dataset_colours <- c(
  "ASO" = "#0072B2",
  "ANO" = "#D55E00"
)

p <- ggplot(
  evaluation_locs
) +
  
  geom_sf(
    aes(colour = dataset),
    size = 0.7,
    alpha = 0.7
  ) +
  
  scale_colour_manual(
    values = dataset_colours,
    name = "Dataset"
  ) +
  
  labs(
    x = "Easting",
    y = "Northing",
  ) +
  
  theme_bw() +
  
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 10)
  )

p

## --------------------------------------------------------
## Save plot
## --------------------------------------------------------

plot_file <- file.path(
  resultFolder, "figures",
  "evaluation_locations_with_richness"
)

## PDF
ggsave(
  filename = paste0(plot_file, ".pdf"),
  plot = p,
  width = 10,
  height = 8,
  units = "in"
)

## PNG
ggsave(
  filename = paste0(plot_file, ".png"),
  plot = p,
  width = 10,
  height = 8,
  units = "in",
  dpi = 600
)

# ==============================================================================
# 10. DEFINE BIOGEOGRAPHIC CLASSIFICATION REGIONS
# ==============================================================================

source("01_codeForEnvAgency/00G_get_nature_regions.R", echo = TRUE)


# ==============================================================================
# 7. LOAD SPECIES DISTRIBUTION RASTERS
# ==============================================================================

## ------------------------------------------------------------------
## Insect occurrence probabilities
## ------------------------------------------------------------------

insects <- terra::rast(
  file.path(
    resultFolder, "modelFit",
    "insectProbabilities.tiff"
  )
)


## ------------------------------------------------------------------
## Plant occurrence probabilities
## ------------------------------------------------------------------

plants <- terra::rast(
  file.path(
    resultFolder, "modelFit",
    "plantProbabilities.tiff"
  )
)

## ========================================================
## Expected plant species richness
## ========================================================
if(!file.exists(  filename = file.path(
  resultFolder,
  "plant_expected_richness_50m.tif"
))){
plant_expected_richness <- terra::app(
  plants,
  fun = sum,
  na.rm = TRUE,
  cores = 1,
  filename = file.path(
    resultFolder,
    "plant_expected_richness.tif"
  ),
  overwrite = TRUE
)

names(plant_expected_richness) <- "plant_expected_richness"

## Use the first calibration layer as the exact target grid
template50 <- allCalibrationRaster[[1]]

plant_expected_richness_50m <- terra::project(
  plant_expected_richness,
  template50,
  method = "bilinear",
  filename = file.path(
    resultFolder,
    "plant_expected_richness_50m.tif"
  ),
  overwrite = TRUE
)
}

plant_expected_richness_50m <- terra::rast(file.path(
  resultFolder,
  "plant_expected_richness_50m.tif"
))

# ==============================================================================
# 11. CALCULATE INTERACTION-SUPPORTED POLLINATOR RICHNESS
# ==============================================================================

if(!file.exists(file.path(
  resultFolder,
  "weighted_insects_probability.tif"
))){
source("01_codeForEnvAgency/01_richnessEstimation.R")

terra::writeRaster(weighted_insects_probability$insect_values,
            filename = file.path(
              resultFolder,
              "weighted_insects_probability.tif"
            ),
            overwrite = TRUE)
}

weighted_insects_probability <- terra::rast(file.path(
  resultFolder,
  "weighted_insects_probability.tif"
))

# ==============================================================================
# 12. CALCULATE INDICATORS FOR EACH SPECIES GROUP
# ==============================================================================

n_groups <- length(vars_of_interest)

message(
  "\nStarting indicator calculations for ",
  n_groups,
  " groups...\n"
)


## ------------------------------------------------------------------
## Initialise progress bar
## ------------------------------------------------------------------

progress_bar <- utils::txtProgressBar(
  min = 0,
  max = n_groups,
  style = 3
)


#for (i in seq_along(vars_of_interest)) {
  i <- 1
  var_of_interest <- vars_of_interest[i]
  
  
  ## ----------------------------------------------------------------
  ## Identify current group
  ## ----------------------------------------------------------------
  
  message(
    "\n====================================================\n",
    "Processing group ",
    i,
    " of ",
    n_groups,
    ": ",
    var_of_interest,
    "\n===================================================="
  )
  
  
  # ------------------------------------------------------------------
  # Select species belonging to the current group
  # ------------------------------------------------------------------
  
  if (
    var_of_interest %in%
    c(
      "Butterflies",
      "Hoverflies",
      "Bees"
    )
  ) {
    
    species_names <- intersect(
      
      gsub(
        " ",
        "_",
        speciesForTaxon$simpleScientificName[
          speciesForTaxon$Taxon == var_of_interest
        ]
      ),
      
      names(insects)
    )
    
    
  } else if (
    var_of_interest %in%
    c(
      "alpine",
      "temperate",
      "boreal"
    )
  ) {
    
    species_names <- intersect(
      
      gsub(
        " ",
        "_",
        speciesForTaxon$simpleScientificName[
          speciesForTaxon$zone == var_of_interest
        ]
      ),
      
      names(insects)
    )
    
    
  } else if (var_of_interest == "all") {
    
    species_names <- names(insects)
    
    
  } else {
    
    warning(
      "Unknown group: ",
      var_of_interest,
      ". Skipping."
    )
    
    utils::setTxtProgressBar(
      progress_bar,
      i
    )
    
    next
  }
  
  
  # ------------------------------------------------------------------
  # Report number of species included
  # ------------------------------------------------------------------
  
  message(
    "Species included: ",
    length(species_names)
  )
  
  
  # ------------------------------------------------------------------
  # Create output directory
  # ------------------------------------------------------------------
  
  plotFolder <- file.path(
    resultFolder,
    var_of_interest
  )
  
  if (!dir.exists(plotFolder)) {
    
    dir.create(
      plotFolder,
      recursive = TRUE
    )
  }
  
  
  # ------------------------------------------------------------------
  # Calculate scaled indicator
  # ------------------------------------------------------------------
  rerunCalibrationModel <- FALSE
  rerunIndicators <- FALSE
  if(rerunIndicators){
  source(
    "01_codeForEnvAgency/02_scaledIndicator.R"
  )
  }
  
  
  # ------------------------------------------------------------------
  # Calculate indices and generate plots
  # ------------------------------------------------------------------
  
  source(
    "01_codeForEnvAgency/03b_get_indices_plot.R"
  )
  
  # ------------------------------------------------------------------
  # Evaluate the indices
  # ------------------------------------------------------------------
  
  source("01_codeForEnvAgency/05_plant_evaluation.R")
  
  source("01_codeForEnvAgency/05_insect_evaluation.R")
  
  # ------------------------------------------------------------------
  # Update progress bar
  # ------------------------------------------------------------------
  
  utils::setTxtProgressBar(
    progress_bar,
    i
  )
#}


# ==============================================================================
# 13. CLOSE PROGRESS BAR
# ==============================================================================

close(progress_bar)


# ==============================================================================
# 14. COMPLETION MESSAGE
# ==============================================================================

message(
  "\n\n====================================================\n",
  "Pollinator indicator analysis completed successfully.\n",
  "Processed groups: ",
  n_groups,
  "\n====================================================\n"
)
