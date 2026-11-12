
library(sf)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(scales)

# ============================================================
# 1. Prepare spatial data
# ============================================================
regional_values <- sf::st_read(file.path(plotFolder, "regional_values_ano.gpkg"))
insect_indicator_reference_values <- readRDS(file.path(plotFolder, "insect_indicator_reference_values.rds"))

plot_data <- regional_values |>
  mutate(
    x = sf::st_coordinates(geometry)[, 1],
    y = sf::st_coordinates(geometry)[, 2]
  ) |>
  sf::st_drop_geometry()

# Ensure the required variables are present
required_columns <- c(
  "x", "y", "Richness",
  "Calibrated_richness",
  "Calibrated regional indicator"
)

missing_columns <- setdiff(required_columns, names(plot_data))

if (length(missing_columns) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

# Shared colour limits for observed and calibrated richness
richness_limits <- range(
  c(
    plot_data$Richness,
    plot_data$Calibrated_richness
  ),
  na.rm = TRUE
)

# ============================================================
# 2. Top row: observed and calibrated richness
# ============================================================

richness_long <- plot_data |>
  select(
    x, y,
    Richness,
    Calibrated_richness
  ) |>
  pivot_longer(
    cols = c(Richness, Calibrated_richness),
    names_to = "variable",
    values_to = "value"
  ) |>
  mutate(
    variable = recode(
      variable,
      Richness = "Estimated richness",
      Calibrated_richness = "Calibrated richness"
    ),
    variable = factor(
      variable,
      levels = c(
        "Estimated richness",
        "Calibrated richness"
      )
    )
  )

p_richness <- ggplot(
  richness_long,
  aes(x = x, y = y, colour = value)
) +
  geom_point(size = 1.2, na.rm = TRUE) +
  facet_wrap(~variable, nrow = 1) +
  scale_colour_viridis_c(
    name = "Species richness",
    limits = richness_limits,
    oob = scales::squish
  ) +
  coord_equal() +
  labs(
    title = "Observed and calibrated richness"
  ) +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey90"),
    legend.position = "right"
  )

# ============================================================
# 3. Middle row: calibrated regional indicator
# ============================================================

p_indicator <- ggplot(
  plot_data,
  aes(
    x = x,
    y = y,
    colour = `Calibrated regional indicator`
  )
) +
  geom_point(size = 1.2, na.rm = TRUE) +
  scale_colour_viridis_c(
    name = "Scaled indicator",
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    labels = scales::label_number(accuracy = 0.1),
    oob = scales::squish
  ) +
  coord_equal() +
  labs(
    title = "Calibrated regional indicator"
  ) +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right"
  )

# ============================================================
# 4. Prepare regional reference values
# ============================================================

required_reference_columns <- c(
  "nature_region",
  "richness", "min_richness", "max_richness",
  "cal_richness", "min_cal_richness", "max_cal_richness"
)

missing_reference_columns <- setdiff(
  required_reference_columns,
  names(insect_indicator_reference_values)
)

if (length(missing_reference_columns) > 0) {
  stop(
    "Missing reference columns: ",
    paste(missing_reference_columns, collapse = ", ")
  )
}

reference_summary <- insect_indicator_reference_values |>
  select(
    nature_region,
    richness,
    min_richness,
    max_richness,
    cal_richness,
    min_cal_richness,
    max_cal_richness
  ) |>
  pivot_longer(
    cols = -nature_region,
    names_to = "statistic",
    values_to = "value"
  ) |>
  mutate(
    # The richness and cal_richness columns are the reference medians.
    measure = case_when(
      statistic %in% c(
        "richness",
        "min_richness",
        "max_richness"
      ) ~ "Observed richness reference",

      statistic %in% c(
        "cal_richness",
        "min_cal_richness",
        "max_cal_richness"
      ) ~ "Calibrated richness reference"
    ),

    statistic = case_when(
      statistic %in% c(
        "richness",
        "cal_richness"
      ) ~ "Median",

      statistic %in% c(
        "min_richness",
        "min_cal_richness"
      ) ~ "Minimum",

      statistic %in% c(
        "max_richness",
        "max_cal_richness"
      ) ~ "Maximum"
    )
  ) |>
  select(
    nature_region,
    measure,
    statistic,
    value
  ) |>
  pivot_wider(
    names_from = statistic,
    values_from = value
  ) |>
  mutate(
    nature_region = factor(
      nature_region,
      levels = c(
        "Central", "East", "North", "South", "West"
      )
    ),
    measure = factor(
      measure,
      levels = c(
        "Observed richness reference",
        "Calibrated richness reference"
      )
    )
  )

# Inspect the summary to verify that medians are populated
print(reference_summary)

# ============================================================
# 5. Bottom row: reference medians and ranges
# ============================================================

p_reference <- ggplot(reference_summary, aes(x = nature_region)) +
  geom_linerange(
    aes(ymin = Minimum, ymax = Maximum),
    linewidth = 0.8,
    colour = "grey50",
    na.rm = TRUE
  ) +
  # Mark the lower and upper tips of each linerange
  geom_point(
    aes(y = Minimum),
    shape = 25, size = 3,
    fill = "#0072B2",
    colour = "#0072B2",
    na.rm = TRUE
  ) +
  geom_point(
    aes(y = Maximum),
    shape = 24, size = 3,
    fill = "#D55E00",
    colour = "#D55E00",
    na.rm = TRUE
  ) +
  # Show the median as a circle
  geom_point(
    aes(y = Median),
    shape = 21, size = 3,
    fill = "black",
    colour = "black",
    na.rm = TRUE
  ) +
  facet_wrap(~measure, scales = "free_y", nrow = 1) +
  labs(
    title = "Regional reference values used for scaling",
    subtitle = "Triangles mark the minimum and maximum tips; circles mark medians",
    x = "Nature region",
    y = "Reference richness"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    strip.background = element_rect(fill = "grey90"),
    panel.grid.minor = element_blank()
  )

# ============================================================
# 6. Combine all three rows
# ============================================================

final_plot <- (
  p_richness /
  p_indicator /
  p_reference
) +
  plot_layout(
    heights = c(1.1, 1, 1)
  ) +
  plot_annotation(
    title = paste(
      "Spatial comparison of richness,",
      "calibrated indicator and reference conditions"
    ),
    subtitle = paste(
      "ANO sampling locations.",
      "The richness panels share a colour scale;",
      "the indicator uses a 0–1 scale."
    )
  )

# Display the final figure
print(final_plot)

# Optional: save the figure
ggsave(
  filename = file.path(
    plotFolder,
    "ANO_richness_indicator_reference_plot.png"
  ),
  plot = final_plot,
  width = 15,
  height = 16,
  units = "in",
  dpi = 300,
  bg = "white"
)


# All rasts at NiN Polygons
input_raster <- file.path(plotFolder, "all_indicator_raster.tif")
all_rasts <- terra::rast(input_raster)
ref_locs
insect_indicator_region_lookup <- readRDS(file.path(plotFolder, "insect_indicator_region_lookup.rds")) %>%
  select(region_id, nature_region)
  

library(terra)
library(sf)
library(dplyr)
library(tidyr)
library(ggplot2)

# 1. Convert polygons to SpatVector
ref_locs_vect <- terra::vect(ref_locs)

# 2. Extract raster values within each polygon
# Exact = FALSE includes raster cells touched by each polygon
extracted <- terra::extract(
  all_rasts,
  ref_locs_vect,
  cells = FALSE,
  ID = TRUE
) 

  reference_polygon_values <- extracted%>%
  dplyr::select(calibrated_richness, richness, region_id) %>%
    na.omit()

library(dplyr)
library(tidyr)
library(ggplot2)

# Add original region names
reference_polygon_values_named <- reference_polygon_values %>%
  left_join(
    insect_indicator_region_lookup,
    by = "region_id"
  ) %>%
  mutate(
    nature_region = factor(
      nature_region,
      levels = c("Central", "East", "North", "South", "West")
    )
  )

# Reshape for plotting
reference_polygon_values_long <- reference_polygon_values_named %>%
  pivot_longer(
    cols = c(calibrated_richness, richness),
    names_to = "indicator",
    values_to = "value"
  ) %>%
  mutate(
    indicator = recode(
      indicator,
      calibrated_richness = "Calibrated richness",
      richness = "Observed richness"
    )
  )

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

# Add original region names
reference_polygon_values_named <- reference_polygon_values %>%
  left_join(insect_indicator_region_lookup, by = "region_id") %>%
  mutate(
    nature_region = factor(
      nature_region,
      levels = c("Central", "East", "North", "South", "West")
    )
  )

# Convert richness variables to long format
reference_polygon_values_long <- reference_polygon_values_named %>%
  pivot_longer(
    cols = c(calibrated_richness, richness),
    names_to = "indicator",
    values_to = "value"
  ) %>%
  mutate(
    indicator = recode(
      indicator,
      calibrated_richness = "Calibrated richness",
      richness = "Estimated richness"
    )
  )

# Row 1: Boxplots
p_boxplot <- ggplot(
  reference_polygon_values_long,
  aes(x = nature_region, y = value)
) +
  geom_boxplot(outlier.size = 0.6, na.rm = TRUE) +
  facet_wrap(~indicator, scales = "free_y", nrow = 1) +
  labs(
    title = "Distribution of richness values by region",
    x = NULL,
    y = "Richness"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    strip.background = element_rect(fill = "grey90"),
    panel.grid.minor = element_blank()
  )

# Calculate requested quantiles
# Calculate quantiles and summary statistics by region and indicator
quantile_summary <- reference_polygon_values_long %>%
  group_by(nature_region, indicator) %>%
  summarise(
    Mean = mean(value, na.rm = TRUE),
    Minimum = min(value, na.rm = TRUE),
    `50th` = quantile(value, 0.50, na.rm = TRUE),
    `75th` = quantile(value, 0.75, na.rm = TRUE),
    `80th` = quantile(value, 0.80, na.rm = TRUE),
    `90th` = quantile(value, 0.90, na.rm = TRUE),
    `95th` = quantile(value, 0.95, na.rm = TRUE),
    `99th` = quantile(value, 0.99, na.rm = TRUE),
    Maximum = max(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c(
      Mean, Minimum, `50th`, `75th`, `80th`,
      `90th`, `95th`, `99th`, Maximum
    ),
    names_to = "statistic",
    values_to = "value"
  ) %>%
  mutate(
    statistic = factor(
      statistic,
      levels = c(
        "Minimum", "Mean", "50th", "75th", "80th",
        "90th", "95th", "99th", "Maximum"
      )
    )
  )

# Row 2: Quantiles and summary statistics
p_quantiles <- ggplot(
  quantile_summary,
  aes(
    x = nature_region,
    y = value,
    colour = statistic,
    shape = statistic,
    linetype = statistic,
    group = statistic
  )
) +
  geom_line(
    #position = position_dodge(width = 0.65),
    linewidth = 0.8
  ) +
  geom_point(
    #position = position_dodge(width = 0.65),
    size = 3
  ) +
  facet_wrap(~indicator, scales = "free_y", nrow = 1) +
  scale_colour_manual(
    name = "Statistic",
    values = c(
      "Minimum" = "#000000",
      "Mean"    = "#E69F00",
      "50th"    = "#0072B2",
      "75th"    = "#009E73",
      "80th"    = "#CC79A7",
      "90th"    = "#D55E00",
      "95th"    = "#56B4E9",
      "99th"    = "#6A3D9A",
      "Maximum" = "#666666"
    )
  ) +
  scale_shape_manual(
    name = "Statistic",
    values = c(
      "Minimum" = 4,
      "Mean"    = 16,
      "50th"    = 17,
      "75th"    = 15,
      "80th"    = 18,
      "90th"    = 8,
      "95th"    = 3,
      "99th"    = 7,
      "Maximum" = 0
    )
  ) +
  scale_linetype_manual(
    name = "Statistic",
    values = c(
      "Minimum" = "solid",
      "Mean"    = "dashed",
      "50th"    = "solid",
      "75th"    = "dashed",
      "80th"    = "dotted",
      "90th"    = "solid",
      "95th"    = "dashed",
      "99th"    = "dotted",
      "Maximum" = "solid"
    )
  ) +
  labs(
    title = "Richness quantiles and summary statistics by region",
    x = "Nature region",
    y = "Richness"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    strip.background = element_rect(fill = "grey90"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.box = "horizontal"
  ) +
  guides(
    colour = guide_legend(nrow = 3, byrow = TRUE),
    shape = guide_legend(nrow = 3, byrow = TRUE),
    linetype = guide_legend(nrow = 3, byrow = TRUE)
  )

# Combine with the existing boxplot row
final_plot <- p_boxplot / p_quantiles +
  plot_layout(heights = c(1, 1.3)) +
  plot_annotation(
    title = "Regional richness distributions and summary statistics within the NiN polygons"
  )

print(final_plot)
