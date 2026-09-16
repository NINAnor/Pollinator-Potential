
## ========================================================
## 1. Get evaluation locations
## ========================================================

# Plant evaluation locs
plant_eval_locs <- evaluation_locs %>%
  dplyr::rename(obs_richness = plnt_rc) %>%
  dplyr::mutate(taxa = "plants")


# Insect evaluation locs
source("01_codeForEnvAgency/04_process_NINA_insect_data.R")
source("01_codeForEnvAgency/04_process_national_insect_monitoring_data.R")

insect_eval_locs <- bind_rows(nina_insect_richness,
                               national_monitoring_insect_richness) %>%
  dplyr::rename(obs_richness = insct_r) %>%
  dplyr::mutate(taxa = "insects") %>%
  dplyr::select(all_of(names(plant_eval_locs)))

# Put all evaluation locs together
all_eval_locs <- bind_rows(plant_eval_locs,
                           insect_eval_locs)%>%
  dplyr::mutate(
    x = sf::st_coordinates(geometry)[, 1],
    y = sf::st_coordinates(geometry)[, 2]
  )


# ======================================
# Extract estimated richness and indicators
# ======================================

## Extract expected plant richness
indicator_vals <- terra::extract(
  plant_expected_richness,
  terra::vect(all_eval_locs) 
) %>%
  select(-ID) %>%
  bind_cols(terra::extract(
    indicator_results$allRasts,
    terra::vect(all_eval_locs)
  )) %>%
  select(-ID) %>%
  bind_cols(all_eval_locs)

write.csv(
  sf::st_drop_geometry(indicator_vals),
  file.path(plotFolder, "all_eval_locs.csv"),
  row.names = FALSE
)


# ======================================
# 
# ======================================

################################################################################
### Load rasters
################################################################################

raster_for_val_dir <- paste0(
  "C:/Users/kwaku.adjei/OneDrive - NINA/",
  "222608 - Indikatorutvikling for \u00E5pent lavland 2026 - CLC for kwaku"
)

short_project_dir <- system(
  paste0(
    'cmd /c for %I in ("',
    raster_for_val_dir,
    '") do @echo %~sI'
  ),
  intern = TRUE
)

list.files(short_project_dir)

# raster_for_val_dir <- "C:/Users/kwaku.adjei/OneDrive - NINA/222608 - Indikatorutvikling for åpent lavland 2026 - CLC for kwaku"
# dir.exists(raster_for_val_dir)

r.CLC100.class.occ.projected.freqs <- terra::rast(file.path(short_project_dir, "r CLC100 1km class freq.tif"))
Landscape.heterogeneity.CLC <- terra::rast(file.path(short_project_dir, "r CLC100 1km Landscape.heterogeneity.tif"))
names(Landscape.heterogeneity.CLC) <- "CLC_heterogeneity"

r.CLC10.class.occ.projected.freqs <- terra::rast(file.path(short_project_dir, "r CLC10m 1km class freq.tif"))
Landscape.heterogeneity.CLCplus <- terra::rast(file.path(short_project_dir, "r CLC10Plus 1km Landscape.heterogeneity.tif"))
names(Landscape.heterogeneity.CLCplus) <- "CLCplus_heterogeneity"

landScapeDiff <- scale(Landscape.heterogeneity.CLCplus)-scale(Landscape.heterogeneity.CLC)
names(landScapeDiff) <- "lScapeDiff"

landScapeRaster <- c(Landscape.heterogeneity.CLC,
                     Landscape.heterogeneity.CLCplus,
                     landScapeDiff)

# ==============================================================================
# 1. Load CLC rasters
# ==============================================================================

r.CLC100 <- r.CLC100.class.occ.projected.freqs
r.CLC10  <- r.CLC10.class.occ.projected.freqs

# Landscape heterogeneity
LDirCLC100 <- Landscape.heterogeneity.CLC
LDirCLC10  <- Landscape.heterogeneity.CLCplus

names(LDirCLC100) <- "LDirCLC100"
names(LDirCLC10)  <- "LDirCLCplus10"


# ==============================================================================
# 2. CLC100 -> CLC+ Backbone crosswalk
# ==============================================================================

CLC100_to_CLCplus <- c(
  `1`  = 1,
  `2`  = 1,
  `3`  = 1,
  `4`  = 1,
  `5`  = 1,
  `6`  = 1,
  `7`  = 9,
  `8`  = 9,
  `9`  = 9,
  `10` = 6,
  `11` = 6,
  `12` = 7,
  `18` = 6,
  `20` = 7,
  `21` = 6,
  `23` = 3,
  `24` = 2,
  `25` = NA,   # Mixed forest
  `26` = 6,
  `27` = 5,
  `29` = 5,
  `30` = 9,
  `31` = 9,
  `32` = 9,
  `33` = 9,
  `34` = 11,
  `35` = 6,
  `36` = 6,
  `39` = 9,
  `40` = 10,
  `41` = 10,
  `43` = 10,
  `44` = 10
)


# ==============================================================================
# 3. Aggregate CLC100 class frequencies to CLC+ backbone classes
# ==============================================================================

CLC100.class <- as.integer(
  sub("^CLC100_", "", names(r.CLC100))
)

CLCplus.class <- unname(
  CLC100_to_CLCplus[as.character(CLC100.class)]
)

keep <- !is.na(CLCplus.class)

target.classes <- sort(unique(CLCplus.class[keep]))

r.CLC100.to.CLCplus <- tapp(
  r.CLC100[[keep]],
  index = CLCplus.class[keep],
  fun = sum,
  na.rm = TRUE
)

names(r.CLC100.to.CLCplus) <- paste0(
  "CLCplus_",
  target.classes
)


# ==============================================================================
# 4. Check that the two 1-km products have matching geometry
# ==============================================================================

stopifnot(
  compareGeom(
    r.CLC100.to.CLCplus,
    r.CLC10,
    stopOnError = FALSE
  )
)


# ==============================================================================
# 5. Calculate CLC100 -> CLC+ differences
# ==============================================================================

# Extract class numbers from CLC+10 layer names
CLC10.class <- as.integer(
  sub("^CLCPlus10_", "", names(r.CLC10))
)

# Keep only CLC+10 classes that also exist in aggregated CLC100
keep <- CLC10.class %in% target.classes

CLC10.class <- CLC10.class[keep]

# Numeric layer indices
CLC10_idx <- which(keep)

# Corresponding CLC100 layers
CLC100_idx <- match(
  CLC10.class,
  target.classes
)

# Extract layers
CLC10_common <- r.CLC10[[CLC10_idx]]
CLC100_common <- r.CLC100.to.CLCplus[[CLC100_idx]]

# Calculate differences
CLC_diffs <- CLC10_common - CLC100_common

# Give the difference layers informative names
names(CLC_diffs) <- paste0(
  "CLC_",
  CLC10.class,
  "_diff"
)

CLC_diffs

allCalibrationRaster <- c(landScapeRaster,
                          CLC_diffs)

# CLC100_forest <- 
#   r.CLC100.class.occ.projected.freqs[["CLC100_23"]] +
#   r.CLC100.class.occ.projected.freqs[["CLC100_24"]] +
#   r.CLC100.class.occ.projected.freqs[["CLC100_25"]]
# 
# CLC100_urban <-
#   r.CLC100.class.occ.projected.freqs[["CLC100_1"]] +
#   r.CLC100.class.occ.projected.freqs[["CLC100_2"]] +
#   r.CLC100.class.occ.projected.freqs[["CLC100_3"]] +
#   r.CLC100.class.occ.projected.freqs[["CLC100_4"]] +
#   r.CLC100.class.occ.projected.freqs[["CLC100_5"]] +
#   r.CLC100.class.occ.projected.freqs[["CLC100_6"]] +
#   r.CLC100.class.occ.projected.freqs[["CLC100_7"]] +
#   r.CLC100.class.occ.projected.freqs[["CLC100_8"]] +
#   r.CLC100.class.occ.projected.freqs[["CLC100_9"]] +
#   r.CLC100.class.occ.projected.freqs[["CLC100_10"]] +
#   r.CLC100.class.occ.projected.freqs[["CLC100_11"]]
# 
# CLC100_grassland <-
#   r.CLC100.class.occ.projected.freqs[["CLC100_18"]] +
#   r.CLC100.class.occ.projected.freqs[["CLC100_26"]]
# 
# CLC100_shrub <-
#   r.CLC100.class.occ.projected.freqs[["CLC100_27"]] +
#   r.CLC100.class.occ.projected.freqs[["CLC100_29"]]
# 
# CLC100_water <-
#   r.CLC100.class.occ.projected.freqs[["CLC100_40"]] +
#   r.CLC100.class.occ.projected.freqs[["CLC100_41"]] +
#   r.CLC100.class.occ.projected.freqs[["CLC100_43"]] +
#   r.CLC100.class.occ.projected.freqs[["CLC100_44"]]
# 
# 
# CLC10_forest <-
#   r.CLC10.class.occ.projected.freqs[["CLCPlus10_2"]] +
#   r.CLC10.class.occ.projected.freqs[["CLCPlus10_3"]] 
# 
# CLC10_urban <-
#   r.CLC10.class.occ.projected.freqs[["CLCPlus10_1"]]
# 
# CLC10_shrub <-
#   r.CLC10.class.occ.projected.freqs[["CLCPlus10_5"]]
# 
# CLC10_grassland <-
#   r.CLC10.class.occ.projected.freqs[["CLCPlus10_6"]] +
#   r.CLC10.class.occ.projected.freqs[["CLCPlus10_7"]]
# 
# CLC10_water <-
#   r.CLC10.class.occ.projected.freqs[["CLCPlus10_10"]]
# 
# forest_difference <- CLC10_forest - CLC100_forest
# 
# urban_difference <- CLC10_urban - CLC100_urban
# 
# grassland_difference <- CLC10_grassland - CLC100_grassland
# 
# shrub_difference <- CLC10_shrub - CLC100_shrub
# 
# water_difference <- CLC10_water - CLC100_water
# 
# difference_covariates <- c(
#   forest_difference,
#   urban_difference,
#   grassland_difference,
#   shrub_difference,
#   water_difference
# )
# 
# names(difference_covariates) <- c(
#   "forest_difference",
#   "urban_difference",
#   "grassland_difference",
#   "shrub_difference",
#   "water_difference"
# )
