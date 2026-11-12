################################################################################
# Indicator estimates at ANO survey locations
# Requires: plotFolder, nature_regions, ano_plant_richness (sf points)
# Input raster/layer are identical to 03b (see 03_survey_indicator_helper.R)
################################################################################

source("01_codeForEnvAgency/03_survey_indicator_helper.R")

ano_results <- run_survey_indicator(ano_plant_richness, out_prefix = "ano")
ano_results
