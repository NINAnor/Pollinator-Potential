################################################################################
# Indicator estimates at ASO and ANO locations pooled
# Locations are pooled as-is (not de-duplicated).
################################################################################

source("01_codeForEnvAgency/03_survey_indicator_helper.R")

aso_ano_locations <- dplyr::bind_rows(
  sf::st_sf(geometry = sf::st_geometry(aso_plant_richness)),
  sf::st_sf(geometry = sf::st_geometry(ano_plant_richness))
)

aso_ano_results <- run_survey_indicator(aso_ano_locations, out_prefix = "aso_ano")
aso_ano_results
