################################################################################
# Define NiN reference locations
#
# This script:
#
#   1. Loads the Norwegian NiN habitat polygons from the Miljødirektoratet
#      File Geodatabase.
#
#   2. Selects polygons representing high-quality semi-natural grassland:
#
#          - Main ecosystem = "semi-naturligMark"
#          - Condition      = "god"
#          - Nature type    = one of:
#                "Naturbeitemark"
#                "Slåttemark"
#                "Semi-naturlig eng"
#
#   3. Saves the selected polygons as a Shapefile.
#
#   4. Creates polygon centroids to provide representative point locations.
#
#   5. Produces a diagnostic figure showing:
#          - the original NiN polygons
#          - the corresponding centroids
#
# The selected polygons are subsequently stored in:
#
#     ref_locs
#
################################################################################


# ==============================================================================
# 1. Define the output polygon file
# ==============================================================================

nin_reference_file <- file.path(
  resultFolder,
  "polygons",
  "good_seminatural_grassland_nin.shp"
)


# ==============================================================================
# 2. Extract suitable NiN reference polygons if they do not already exist
# ==============================================================================

if (!file.exists(nin_reference_file)) {
  
  
  # ============================================================================
  # 2.1 Location of the NiN File Geodatabase
  # ============================================================================
  #
  # This is the 2025 NiN dataset from the Norwegian Environment Agency
  # (Miljødirektoratet).
  #
  # The data are stored in EPSG:25833:
  #
  #     ETRS89 / UTM zone 33N
  # ============================================================================
  
  nin_gdb <- paste0(
    "R:/GeoSpatialData/Habitats_biotopes/",
    "Norway_Miljodirektoratet_Naturtyper_nin(instruks)/",
    "Original/verson2025/",
    "Naturtyper - Miljodirektoratets instruks/",
    "Naturtyper_nin_0000_norge_25833_FILEGDB/",
    "Naturtyper_nin_0000_norge_25833_FILEGDB.gdb"
  )
  
  
  # ============================================================================
  # 2.2 Read the NiN polygon layer
  # ============================================================================
  #
  # `naturtyper_nin_omr` contains the mapped NiN nature-type polygons.
  #
  # `quiet = FALSE` allows sf to report information about the layer while it
  # is being read.
  # ============================================================================
  
  nin_polygons <- sf::st_read(
    nin_gdb,
    layer = "naturtyper_nin_omr",
    quiet = FALSE
  )
  
  
  # ============================================================================
  # 2.3 Rename important columns
  # ============================================================================
  #
  # The original NiN dataset contains long or dataset-specific field names.
  # These three fields are renamed for easier use later in the script.
  #
  # IMPORTANT:
  # The current code identifies columns by POSITION:
  #
  #       column 5 -> area_name
  #       column 7 -> main_ecosystem
  #       column 9 -> mapping_year
  #
  # This assumes the NiN data schema does not change.
  #
  # For a production workflow, it is safer to rename columns by their actual
  # names rather than their positions.
  # ============================================================================
  
  names(nin_polygons)[5] <- "area_name"
  names(nin_polygons)[7] <- "main_ecosystem"
  names(nin_polygons)[9] <- "mapping_year"
  
  
  # ============================================================================
  # 2.4 Define the nature types to be used as reference habitats
  # ============================================================================
  #
  # These are the NiN nature types considered to represent semi-natural
  # grassland for the reference dataset.
  # ============================================================================
  
  reference_naturtyper <- c(
    "Naturbeitemark",
    "Slåttemark",
    "Semi-naturlig eng"
  )
  
  
  # ============================================================================
  # 2.5 Select high-quality semi-natural grassland polygons
  # ============================================================================
  #
  # Three criteria are applied:
  #
  #   1. The polygon belongs to the semi-natural ecosystem:
  #
  #          main_ecosystem == "semi-naturligMark"
  #
  #   2. Its mapped ecological condition is good:
  #
  #          tilstand == "god"
  #
  #   3. Its nature type belongs to the predefined reference nature types.
  #
  # Only the nature type and nature-type code are retained in the output.
  #
  # These polygons therefore represent the selected "good" semi-natural
  # grassland reference locations.
  # ============================================================================
  
  good_seminatural_grassland <- nin_polygons %>%
    dplyr::filter(
      main_ecosystem == "semi-naturligMark",
      tilstand == "god",
      naturtype %in% reference_naturtyper
    ) %>%
    dplyr::select(
      naturtype,
      naturtypeKode
    )
  
  
  # ============================================================================
  # 2.6 Store selected polygons as reference locations
  # ============================================================================
  
  ref_locs <- good_seminatural_grassland
  
  
  # ============================================================================
  # 2.7 Save the selected polygons as a Shapefile
  # ============================================================================
  
  sf::st_write(
    good_seminatural_grassland,
    nin_reference_file,
    delete_layer = TRUE,
    quiet = FALSE
  )
  
  
  # ============================================================================
  # 2.8 Reload the saved reference polygons
  # ============================================================================
  #
  # Reloading the file ensures that `ref_locs` represents the actual saved
  # dataset that will be used by subsequent analyses.
  # ============================================================================
  
  ref_locs <- sf::st_read(
    nin_reference_file,
    quiet = TRUE
  )
  
  
} else {
  
  
  # ============================================================================
  # 3. If the reference polygons already exist, simply load them
  # ============================================================================
  
  ref_locs <- sf::st_read(
    nin_reference_file,
    quiet = TRUE
  )
}


# ==============================================================================
# 4. Load packages used for visualisation
# ==============================================================================

library(ggplot2)
library(sf)
library(patchwork)


# ==============================================================================
# 5. Define common colours for the NiN nature types
# ==============================================================================
#
# The same colours are used for both the polygon and centroid maps so that
# the two panels can be compared directly.
# ==============================================================================

nature_type_colours <- c(
  "Naturbeitemark" = "#E76F51",
  "Semi-naturlig eng" = "#00A6A6"
)


# ==============================================================================
# 6. Plot the original NiN polygons
# ==============================================================================

p1 <- ggplot() +
  
  # --------------------------------------------------------------------------
# Draw the NiN reference polygons
# --------------------------------------------------------------------------

geom_sf(
  data = ref_locs,
  aes(fill = natrtyp),
  colour = "black",
  linewidth = 0.1
) +
  
  # --------------------------------------------------------------------------
# Apply the predefined nature-type colours
# --------------------------------------------------------------------------

scale_fill_manual(
  values = nature_type_colours,
  name = "Nature type"
) +
  
  # --------------------------------------------------------------------------
# Use the native spatial extent without additional expansion
# --------------------------------------------------------------------------

coord_sf(
  expand = FALSE
) +
  
  # --------------------------------------------------------------------------
# Axis labels and title
# --------------------------------------------------------------------------

labs(
  x = "Easting",
  y = "Northing",
  title = "NiN polygons"
) +
  
  # --------------------------------------------------------------------------
# Use a clean black-and-white theme
# --------------------------------------------------------------------------

theme_bw() +
  
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    plot.title = element_text(
      face = "bold",
      size = 12
    )
  )


# ==============================================================================
# 7. Calculate polygon centroids
# ==============================================================================
#
# A centroid is a representative point calculated from each polygon.
#
# These points can subsequently be used as reference locations in analyses
# that require point coordinates rather than polygons.
#
# NOTE:
# `st_centroid()` calculates the geometric centroid. For irregular or
# highly concave polygons, the centroid can sometimes fall outside the
# polygon. If the point must always lie inside the polygon, consider using
# `st_point_on_surface()` instead.
# ==============================================================================

ref_points <- ref_locs %>%
  sf::st_centroid()


# ==============================================================================
# 8. Plot the polygon centroids
# ==============================================================================

p2 <- ggplot() +
  
  # --------------------------------------------------------------------------
# Draw one point for each NiN polygon centroid
# --------------------------------------------------------------------------

geom_sf(
  data = ref_points,
  aes(fill = natrtyp),
  shape = 21,
  colour = "black",
  size = 2,
  alpha = 0.8
) +
  
  # --------------------------------------------------------------------------
# Use the same nature-type colours as the polygon map
# --------------------------------------------------------------------------

scale_fill_manual(
  values = nature_type_colours,
  name = "Nature type"
) +
  
  coord_sf(
    expand = FALSE
  ) +
  
  labs(
    x = "Easting",
    y = "Northing",
    title = "Centroids"
  ) +
  
  theme_bw() +
  
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    plot.title = element_text(
      face = "bold",
      size = 12
    )
  )


# ==============================================================================
# 9. Combine the two maps
# ==============================================================================
#
# `patchwork` places the two plots side-by-side.
#
# `guides = "collect"` ensures that only one shared legend is displayed
# rather than having a separate legend for each panel.
# ==============================================================================

p <- patchwork::wrap_plots(
  p1,
  p2,
  ncol = 2,
  guides = "collect"
) &
  theme(
    legend.position = "right"
  )


# Display the combined figure
p


# ==============================================================================
# 10. Save the diagnostic figure
# ==============================================================================

plot_file <- file.path(
  resultFolder,
  "figures",
  "good_seminatural_grassland_nin"
)


# ------------------------------------------------------------------------------
# 10.1 Save PDF
# ------------------------------------------------------------------------------

ggplot2::ggsave(
  filename = paste0(plot_file, ".pdf"),
  plot = p,
  width = 10,
  height = 10,
  units = "in"
)


# ------------------------------------------------------------------------------
# 10.2 Save high-resolution PNG
# ------------------------------------------------------------------------------

ggplot2::ggsave(
  filename = paste0(plot_file, ".png"),
  plot = p,
  width = 10,
  height = 10,
  units = "in",
  dpi = 600
)
