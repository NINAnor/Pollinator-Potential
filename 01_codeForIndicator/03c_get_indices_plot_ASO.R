################################################################################
# Indicator estimates at ASO survey locations
# Requires: plotFolder, nature_regions, aso_plant_richness (sf points)
# Input raster/layer are identical to 03b (see 03_survey_indicator_helper.R)
################################################################################

source("01_codeForEnvAgency/03_survey_indicator_helper.R")

aso_results <- run_survey_indicator(aso_plant_richness, out_prefix = "aso")
aso_results
