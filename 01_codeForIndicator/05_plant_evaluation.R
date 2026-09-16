
################################################################################
# Plant indicator evaluation
#
# Purpose:
#   1. Extract indicator predictions and expected plant richness at evaluation
#      locations.
#   2. Compare observed plant richness with the plant and pollinator indicators.
#   3. Quantify correlations between expected plant richness and each indicator.
#   4. Evaluate the effect of calibration on the regional and global indicators.
#
# Main outputs:
#   - plant_indicator_evaluation.pdf
#   - plant_indicator_evaluation.png
#   - plant_calibration_difference_tests.csv
#   - plant_calibration_difference_effect.pdf
#   - plant_calibration_difference_effect.png
################################################################################


# ==============================================================================
# 0. Load required packages and define plotting colours
# ==============================================================================

library(terra)
library(sf)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(lme4)
library(emmeans)
library(readr)


# Colours used to distinguish the two datasets throughout the analysis.
dataset_colours <- c(
  "ASO" = "#0072B2",
  "ANO" = "#D55E00"
)


# ==============================================================================
# 1. Extract indicator predictions at evaluation locations
# ==============================================================================

# ------------------------------------------------------------------------------
# 1.1 Load the raster containing all indicator predictions
# ------------------------------------------------------------------------------

allRasts <- terra::rast(
  file.path(
    plotFolder,
    "all_indicator_raster.tif"
  )
)


# ------------------------------------------------------------------------------
# 1.2 Extract indicator values at evaluation locations
# ------------------------------------------------------------------------------

# Convert evaluation locations from sf to SpatVector and extract the values
# from every indicator layer.
indicator_values <- terra::extract(
  allRasts,
  terra::vect(evaluation_locs)
)


# ------------------------------------------------------------------------------
# 1.3 Extract expected plant richness
# ------------------------------------------------------------------------------

# Extract the expected plant richness surface at the same evaluation
# locations. This provides the reference quantity against which the
# indicators are evaluated.
plant_richness_values <- terra::extract(
  plant_expected_richness_50m,
  terra::vect(evaluation_locs)
)


# ------------------------------------------------------------------------------
# 1.4 Assign nature regions to evaluation locations
# ------------------------------------------------------------------------------

# Spatially join each evaluation location to its corresponding nature region.
# Only the nature_region variable is retained because the remaining spatial
# attributes are not required for the subsequent analyses.
regions_eval_locs <- sf::st_join(
  evaluation_locs,
  nature_regions %>%
    dplyr::select(nature_region),
  join = sf::st_within
)


# ==============================================================================
# 2. Construct the evaluation dataset
# ==============================================================================

# Combine:
#   - indicator predictions,
#   - expected plant richness,
#   - evaluation-location metadata, and
#   - nature-region information.
#
# Only variables required for the analyses below are retained.

rt <- indicator_values %>%
  
  # Add expected plant richness.
  bind_cols(
    plant_richness_values %>%
      select(-ID)
  ) %>%
  
  # Add evaluation-location metadata and nature region.
  bind_cols(
    regions_eval_locs
  ) %>%
  
  # Retain only variables required for evaluation.
  select(
    richness,
    calibrated_richness,
    richness_without_interaction,
    regional_clamped,
    calibrated_regional_clamped,
    regional_clamped_without_interaction,
    global_clamped,
    calibrated_global_clamped,
    global_clamped_without_interaction,
    plnt_rc,
    plant_expected_richness,
    dataset,
    nature_region
  ) %>%
  
  # Remove locations with missing values from the evaluation dataset.
  na.omit()


# ==============================================================================
# 3. Prepare data for indicator–plant richness comparisons
# ==============================================================================

# Remove geometry if rt is an sf object. Geometry is not required for the
# correlation analysis or scatter plots.
plot_data <- rt %>%
  sf::st_drop_geometry() %>%
  
  select(
    plant_expected_richness,
    plnt_rc,
    calibrated_regional_clamped,
    calibrated_richness,
    richness,
    dataset
  ) %>%
  
  filter(
    !is.na(plant_expected_richness),
    !is.na(dataset)
  )


# ==============================================================================
# 4. Calculate correlations between plant richness and indicators
# ==============================================================================

# ------------------------------------------------------------------------------
# Function: get_correlation()
#
# Calculates Pearson's correlation coefficient between two variables.
#
# Arguments:
#
#   data  - data frame containing the variables
#   x_var - name of the x-axis variable
#   y_var - name of the y-axis variable
#
# Returns:
#
#   A formatted text label containing Pearson's correlation coefficient.
# ------------------------------------------------------------------------------

get_correlation <- function(
    data,
    x_var,
    y_var
) {
  
  # Keep only complete observations for the two variables.
  data %>%
    filter(
      !is.na(.data[[x_var]]),
      !is.na(.data[[y_var]])
    ) %>%
    
    # Calculate the correlation independently for ASO and ANO.
    group_by(dataset) %>%
    
    summarise(
      r = cor(
        .data[[x_var]],
        .data[[y_var]],
        method = "pearson"
      ),
      .groups = "drop"
    ) %>%
    
    # Create a label for displaying the correlation on the plot.
    mutate(
      label = paste0(
        dataset,
        ": r = ",
        sprintf("%.2f", r)
      )
    )
}



# ==============================================================================
# 5. Function for generating indicator scatter plots
# ==============================================================================

# ------------------------------------------------------------------------------
# Function: make_scatter()
#
# Produces a standardised scatter plot showing the relationship between
# expected plant richness and an indicator.
#
# The plot includes:
#   - observations coloured by dataset,
#   - dataset-specific linear regression lines,
#   - Pearson correlation coefficients, and
#   - a panel label (A–D).
# ------------------------------------------------------------------------------

make_scatter <- function(
    data,
    x_var,
    y_var,
    x_label,
    y_label,
    panel_label
) {
  
  # Calculate dataset-specific Pearson correlations.
  correlations <- get_correlation(
    data = data,
    x_var = x_var,
    y_var = y_var
  )
  
  ggplot(
    data,
    aes(
      x = .data[[x_var]],,
      y = .data[[y_var]],
      colour = dataset
    )
  ) +
    
    # --------------------------------------------------------------------------
  # Observed/predicted locations
  # --------------------------------------------------------------------------
  geom_point(
    alpha = 0.6,
    size = 1.8
  ) +
    
    # --------------------------------------------------------------------------
  # Dataset-specific linear relationships
  # --------------------------------------------------------------------------
  geom_smooth(
    method = "lm",
    se = FALSE,
    linewidth = 0.8
  ) +
    
    # --------------------------------------------------------------------------
  # Pearson correlation coefficients
  # --------------------------------------------------------------------------
  #
  # Inf places the labels in the upper-right corner. The vertical position
  # is varied so that ASO and ANO labels do not overlap.
  geom_text(
    data = correlations,
    aes(
      x = Inf,
      y = Inf,
      label = label,
      colour = dataset
    ),
    hjust = 1.1,
    vjust = seq(
      from = 1.2,
      to = 2.4,
      length.out = nrow(correlations)
    ),
    inherit.aes = FALSE,
    fontface = "bold",
    size = 3.5
  ) +
    
    # --------------------------------------------------------------------------
  # Panel label
  # --------------------------------------------------------------------------
  annotate(
    "text",
    x = -Inf,
    y = Inf,
    label = panel_label,
    hjust = -0.5,
    vjust = 1.5,
    fontface = "bold",
    size = 5
  ) +
    
    # --------------------------------------------------------------------------
  # Dataset colours and axis labels
  # --------------------------------------------------------------------------
  scale_colour_manual(
    values = dataset_colours
  ) +
    
    labs(
      x = x_label,
      y = y_label,
      colour = "Dataset"
    ) +
    
    # --------------------------------------------------------------------------
  # Plot theme
  # --------------------------------------------------------------------------
  theme_classic(
    base_size = 12
  ) +
    
    theme(
      legend.position = "bottom",
      axis.title = element_text(face = "bold"),
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        linewidth = 0.5
      )
    )
}


# ==============================================================================
# 6. Generate the four indicator evaluation plots
# ==============================================================================

# A: Comparison with observed plant richness.
# pA <- make_scatter(
#   plot_data,
#   y_var = "plnt_rc",
#   y_label = "Observed plant richness",
#   panel_label = "A"
# )


pA <- make_scatter(
  data = plot_data,
  x_var = "calibrated_richness",
  y_var = "plnt_rc",
  x_label = "Pollinator potential",
  y_label = "Observed plant richness",
  panel_label = "A"
)

# B: Relationship with the calibrated regional indicator.
pB <- make_scatter(
  data = plot_data,
  x_var = "calibrated_regional_clamped",
  y_var = "plnt_rc",
  x_label = "Scaled indicator",
  y_label = " ",
  panel_label = "B"
)


# C: Relationship with pollination potential.
pC <- make_scatter(
  data = plot_data,
  x_var = "richness",
  y_var = "plnt_rc",
  x_label = "Expected insect richness",
  y_label = " ",
  panel_label = "C"
)

# D: Relationship with expected pollinator richness.
# pD <- make_scatter(
#   plot_data,
#   y_var = "richness",
#   y_label = "Expected pollinator richness",
#   panel_label = "D"
# )


# ==============================================================================
# 7. Combine indicator evaluation plots
# ==============================================================================

# Arrange the four panels in a 2 × 2 layout.
#
# guides = "collect" ensures that only one shared legend is displayed.
combined_plot <- (pA | pB | pC) + #/
  # (pC | pD) +
  plot_layout(
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom"
  )


# Display the combined figure.
combined_plot


# ==============================================================================
# 8. Save indicator evaluation figure
# ==============================================================================

# ------------------------------------------------------------------------------
# PDF
# ------------------------------------------------------------------------------

ggsave(
  filename = file.path(
    plotFolder,
    "plant_indicator_evaluation.pdf"
  ),
  plot = combined_plot,
  width = 12,
  height = 12,
  units = "in",
  device = "pdf"
)


# ------------------------------------------------------------------------------
# PNG
# ------------------------------------------------------------------------------

ggsave(
  filename = file.path(
    plotFolder,
    "plant_indicator_evaluation.png"
  ),
  plot = combined_plot,
  width = 12,
  height = 12,
  units = "in",
  dpi = 600,
  device = "png"
)


# ==============================================================================
# 9. Evaluate the effect of calibration
# ==============================================================================

# The following analysis compares the uncalibrated and calibrated versions
# of the regional and global indicators.
#
# Positive difference:
#     uncalibrated indicator > calibrated indicator
#
# Negative difference:
#     calibrated indicator > uncalibrated indicator
#
# A difference of zero indicates no effect of calibration.


# ------------------------------------------------------------------------------
# 9.1 Calculate calibration-induced differences
# ------------------------------------------------------------------------------

indicator_difference <- rt %>%
  
  transmute(
    nature_region,
    dataset,
    
    # Regional indicator:
    # uncalibrated - calibrated
    regional_difference =
      regional_clamped -
      calibrated_regional_clamped,
    
    # Global indicator:
    # uncalibrated - calibrated
    global_difference =
      global_clamped -
      calibrated_global_clamped
  ) %>%
  
  # Convert regional and global differences into a single long-format column.
  pivot_longer(
    cols = c(
      regional_difference,
      global_difference
    ),
    names_to = "indicator_type",
    values_to = "difference"
  ) %>%
  
  # Replace technical variable names with readable labels.
  mutate(
    indicator_type = recode(
      indicator_type,
      regional_difference = "Regional",
      global_difference = "National"
    ),
    
    # Explicitly define the desired plotting/model order.
    indicator_type = factor(
      indicator_type,
      levels = c(
        "Regional",
        "National"
      )
    )
  )


# ==============================================================================
# 10. Statistical tests of calibration effects
# ==============================================================================

# For each indicator type × nature-region combination:
#
#   1. Calculate sample size.
#   2. Calculate the mean calibration-induced difference.
#   3. Calculate the standard deviation.
#   4. Test whether the mean difference differs from zero.
#
# P-values are subsequently adjusted using the Benjamini–Hochberg procedure
# to account for multiple comparisons.

calibration_tests <- indicator_difference %>%
  
  group_by(
    indicator_type,
    nature_region
  ) %>%
  
  summarise(
    
    # Number of non-missing observations.
    n = sum(
      !is.na(difference)
    ),
    
    # Mean difference between uncalibrated and calibrated indicators.
    mean_difference = mean(
      difference,
      na.rm = TRUE
    ),
    
    # Standard deviation of the difference.
    sd_difference = sd(
      difference,
      na.rm = TRUE
    ),
    
    # One-sample t-test:
    # H0: mean difference = 0
    p_value = t.test(
      difference,
      mu = 0
    )$p.value,
    
    .groups = "drop"
  ) %>%
  
  # Control the false discovery rate across the statistical tests.
  mutate(
    p_adjusted = p.adjust(
      p_value,
      method = "BH"
    ),
    
    # Identify statistically significant calibration effects.
    significant = p_adjusted < 0.05
  )


# ------------------------------------------------------------------------------
# 10.1 Save statistical results
# ------------------------------------------------------------------------------

write_csv(
  calibration_tests,
  file.path(
    plotFolder,
    "plant_calibration_difference_tests.csv"
  )
)


# ==============================================================================
# 11. Mixed-effects model of calibration effects
# ==============================================================================

# Model:
#
# difference ~ indicator type × nature region + random effect of dataset
#
# Fixed effects:
#   indicator_type
#       Tests whether calibration effects differ between regional and global
#       indicators.
#
#   nature_region
#       Tests whether calibration effects differ among nature regions.
#
#   indicator_type × nature_region
#       Tests whether the difference between regional and global calibration
#       effects depends on nature region.
#
# Random effect:
#   dataset
#       Accounts for systematic differences between the ASO and ANO datasets.

difference_model <- lmer(
  difference ~
    indicator_type * factor(nature_region) +
    (1 | dataset),
  data = indicator_difference,
  REML = TRUE
)


# ==============================================================================
# 12. Estimated marginal means
# ==============================================================================

# Calculate the model-estimated mean calibration effect for each combination
# of indicator type and nature region.
emm <- emmeans(
  difference_model,
  ~ indicator_type | nature_region
)


# Convert the emmeans object into a data frame and identify whether the
# confidence interval excludes zero.
emm_df <- as.data.frame(emm) %>%
  
  mutate(
    
    # Significant effect if the confidence interval does not overlap zero.
    significant =
      asymp.LCL > 0 |
      asymp.UCL < 0,
    
    # Human-readable significance label.
    effect_label = if_else(
      significant,
      "Significant",
      "Not significant"
    )
  )


# ==============================================================================
# 13. Plot estimated calibration effects
# ==============================================================================

p_effects <- ggplot(
  emm_df,
  aes(
    x = factor(nature_region),
    y = emmean,
    colour = indicator_type
  )
) +
  
  # Reference line representing zero calibration effect.
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  
  # 95% confidence intervals around the estimated marginal means.
  geom_errorbar(
    aes(
      ymin = asymp.LCL,
      ymax = asymp.UCL
    ),
    position = position_dodge(
      width = 0.6
    ),
    width = 0.15,
    linewidth = 0.7
  ) +
  
  # Estimated marginal means.
  geom_point(
    position = position_dodge(
      width = 0.6
    ),
    size = 3
  ) +
  
  labs(
    x = "Nature region",
    y = "Uncalibrated − calibrated indicator",
    colour = "Indicator"
  ) +
  
  theme_bw() +
  
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom",
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  )


# Display the calibration-effect plot.
p_effects


# ==============================================================================
# 14. Save calibration-effect figure
# ==============================================================================

# ------------------------------------------------------------------------------
# PDF
# ------------------------------------------------------------------------------

ggsave(
  filename = file.path(
    plotFolder,
    "plant_calibration_difference_effect.pdf"
  ),
  plot = p_effects,
  width = 10,
  height = 7,
  units = "in",
  device = "pdf"
)


# ------------------------------------------------------------------------------
# PNG
# ------------------------------------------------------------------------------

ggsave(
  filename = file.path(
    plotFolder,
    "plant_calibration_difference_effect.png"
  ),
  plot = p_effects,
  width = 10,
  height = 7,
  units = "in",
  dpi = 600,
  device = "png"
)


################################################################################
# End of script
################################################################################

