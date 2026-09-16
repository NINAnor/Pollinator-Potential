################################################################################
# Define Norwegian nature regions
#
# This script:
#
#   1. Loads Norwegian county boundaries from the national administrative
#      boundary dataset.
#
#   2. Assigns each county to one of five broader nature regions:
#
#          South
#          East
#          West
#          Central
#          North
#
#   3. Dissolves the county boundaries within each nature region to create
#      one polygon per region.
#
#   4. Produces a map showing the resulting nature regions.
#
# The resulting spatial object is:
#
#     nature_regions
#
# It contains one polygon for each nature region and can subsequently be used
# for spatial overlays, mapping, or assigning regional membership to other
# spatial datasets.
#
################################################################################


# ==============================================================================
# 1. Define the path to the Norwegian county geodatabase
# ==============================================================================
#
# The dataset contains the official Norwegian county (fylke) boundaries.
#
# The data are in EPSG:25833:
#
#     ETRS89 / UTM zone 33N
#
# This is useful because the resulting polygons can be directly compared with
# other spatial data already represented in the same projected coordinate
# system.
# ==============================================================================

gdb <- paste0(
  "R:/GeoSpatialData/AdministrativeUnits/",
  "Norway_AdministrativeUnits/Original/Norway_County/",
  "versjon2025/Administrative enheter fylker FGDB-format/",
  "Basisdata_0000_Norge_25833_Fylker_FGDB/",
  "Basisdata_0000_Norge_25833_Fylker_FGDB.gdb"
)


# ==============================================================================
# 2. Read the county polygons
# ==============================================================================
#
# The `fylke` layer contains the individual Norwegian county polygons.
#
# `sf::st_read()` returns these as an sf object, meaning that the dataset
# contains both:
#
#     - attribute information (e.g. county number)
#     - geometry (the county boundaries)
# ==============================================================================

counties <- sf::st_read(
  gdb,
  layer = "fylke"
)


# ==============================================================================
# 3. Assign each county to a broader nature region
# ==============================================================================
#
# The county number (`fylkesnummer`) is used to assign counties to one of
# five broader geographical/nature regions.
#
# `case_when()` evaluates the conditions from top to bottom.
#
# Counties that do not match any of the specified groups receive NA.
# ==============================================================================

nature_regions <- counties %>%
  
  dplyr::mutate(
    
    # --------------------------------------------------------------------------
    # Convert county number to character
    #
    # This ensures that the county codes can be reliably compared with the
    # character values specified below.
    # --------------------------------------------------------------------------
    
    fylkesnummer = as.character(
      fylkesnummer
    ),
    
    
    # --------------------------------------------------------------------------
    # Assign counties to nature regions
    # --------------------------------------------------------------------------
    
    nature_region = dplyr::case_when(
      
      
      # ========================================================================
      # SOUTH
      # ========================================================================
      #
      # Counties:
      #
      #     40
      #     42
      #
      # ========================================================================
      
      fylkesnummer %in% c(
        "40",
        "42"
      ) ~ "South",
      
      
      # ========================================================================
      # EAST
      # ========================================================================
      #
      # Counties:
      #
      #     31
      #     32
      #     03
      #     39
      #     33
      #     34
      #
      # ========================================================================
      
      fylkesnummer %in% c(
        "31",
        "32",
        "03",
        "39",
        "33",
        "34"
      ) ~ "East",
      
      
      # ========================================================================
      # WEST
      # ========================================================================
      #
      # Counties:
      #
      #     46
      #     11
      #
      # ========================================================================
      
      fylkesnummer %in% c(
        "46",
        "11"
      ) ~ "West",
      
      
      # ========================================================================
      # CENTRAL
      # ========================================================================
      #
      # Counties:
      #
      #     50
      #     15
      #
      # ========================================================================
      
      fylkesnummer %in% c(
        "50",
        "15"
      ) ~ "Central",
      
      
      # ========================================================================
      # NORTH
      # ========================================================================
      #
      # Counties:
      #
      #     18
      #     55
      #     56
      #
      # ========================================================================
      
      fylkesnummer %in% c(
        "18",
        "55",
        "56"
      ) ~ "North",
      
      
      # ------------------------------------------------------------------------
      # Any county not explicitly assigned above receives NA.
      # ------------------------------------------------------------------------
      
      TRUE ~ NA_character_
    )
  ) %>%
  
  
  # =============================================================================
# 4. Remove counties that were not assigned to a nature region
# =============================================================================
#
# This ensures that only the five defined nature regions are retained.
# =============================================================================

dplyr::filter(
  !is.na(nature_region)
) %>%
  
  
  # =============================================================================
# 5. Dissolve counties within each nature region
# =============================================================================
#
# Several counties may belong to the same nature region.
#
# `group_by(nature_region)` followed by `summarise()` dissolves their
# boundaries, producing one spatial polygon per nature region.
#
# For example:
#
#     East
#       ├── County A
#       ├── County B
#       ├── County C
#       └── County D
#
# becomes:
#
#     East = one combined spatial polygon
#
# This is useful for producing regional masks and spatial overlays.
# =============================================================================

dplyr::group_by(
  nature_region
) %>%
  
  dplyr::summarise(
    .groups = "drop"
  )


# ==============================================================================
# 6. Plot the nature regions
# ==============================================================================
#
# Each region is displayed using a different fill colour automatically
# determined by ggplot2.
#
# The black boundaries make the regional borders easier to see.
# ==============================================================================

p <- ggplot2::ggplot() +
  
  ggplot2::geom_sf(
    data = nature_regions,
    aes(fill = nature_region),
    colour = "black",
    linewidth = 0.4
  ) +
  
  # --------------------------------------------------------------------------
# Optional region labels
#
# Uncomment this section if you want the name of each region displayed
# directly on the map.
# --------------------------------------------------------------------------

# ggplot2::geom_sf_text(
#   data = nature_regions,
#   aes(label = nature_region),
#   fontface = "bold",
#   size = 4
# ) +

# --------------------------------------------------------------------------
# Legend title
# --------------------------------------------------------------------------

ggplot2::labs(
  fill = "Nature region"
) +
  
  # --------------------------------------------------------------------------
# Use a simple map-friendly theme
# --------------------------------------------------------------------------

ggplot2::theme_bw() +
  
  ggplot2::theme(
    panel.grid = element_blank(),
    legend.position = "right"
  )


# Display the map
p


# ==============================================================================
# 7. Define output path
# ==============================================================================

plot_file <- file.path(
  resultFolder,
  "figures",
  "nature_regions"
)


# ==============================================================================
# 8. Save PDF version
# ==============================================================================

ggplot2::ggsave(
  filename = paste0(plot_file, ".pdf"),
  plot = p,
  width = 8,
  height = 8,
  units = "in"
)


# ==============================================================================
# 9. Save high-resolution PNG version
# ==============================================================================

ggplot2::ggsave(
  filename = paste0(plot_file, ".png"),
  plot = p,
  width = 8,
  height = 8,
  units = "in",
  dpi = 300
)
