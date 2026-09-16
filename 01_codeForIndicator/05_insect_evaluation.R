################################################################################
# Insect indicator evaluation
#
# Purpose:
#
#   1. Prepare observed insect richness from the insect monitoring data.
#   2. Extract indicator predictions at insect monitoring locations.
#   3. Compare indicator values with:
#        - expected plant richness,
#        - observed insect richness, and
#        - expected pollinator richness.
#   4. Quantify the effect of calibration on regional and global indicators.
#
# Main outputs:
#
#   - richness_calibration_relationships.pdf
#   - richness_calibration_relationships.png
#   - insect_calibration_difference_tests.csv
#   - insect_calibration_difference_effect.pdf
#   - insect_calibration_difference_effect.png
################################################################################


# ==============================================================================
# 0. Load required packages
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


# ==============================================================================
# 1. Define dataset colours
# ==============================================================================

# Colours used to distinguish the two insect datasets in figures.
dataset_colours <- c(
  "NINA_insect" = "#0072B2",
  "national_insect_monitoring" = "#D55E00"
)


# ==============================================================================
# 2. Load and prepare observed insect richness
# ==============================================================================

# ------------------------------------------------------------------------------
# 2.1 Load insect monitoring data
# ------------------------------------------------------------------------------

# NINA insect monitoring data can also be loaded here if required.
#
# source("01_codeForEnvAgency/04_process_NINA_insect_data.R")

# Load the national insect monitoring data.
source(
  "01_codeForEnvAgency/04_process_national_insect_monitoring_data.R"
)


# ------------------------------------------------------------------------------
# 2.2 Calculate mean observed insect richness per location
# ------------------------------------------------------------------------------

# Multiple observations may occur at the same local location.
#
# We therefore calculate the mean observed insect richness for each location
# before extracting raster-based indicator predictions.

observed_richness <- bind_rows(
  national_monitoring_insect_richness
) %>%
  
  group_by(local) %>%
  
  summarise(
    insct_r = mean(insct_r),
    .groups = "drop"
  )


# ==============================================================================
# 3. Extract indicator predictions at insect monitoring locations
# ==============================================================================

# ------------------------------------------------------------------------------
# 3.1 Extract all indicator predictions
# ------------------------------------------------------------------------------

# Extract the values of every indicator raster at the observed insect
# monitoring locations.
#
# The resulting object contains one row per monitoring location and one
# column for each indicator raster layer.

indicator_values <- terra::extract(
  allRasts,
  terra::vect(observed_richness)
)


# ------------------------------------------------------------------------------
# 3.2 Extract expected plant richness
# ------------------------------------------------------------------------------

# Extract expected plant richness from the 50-m prediction raster at the same
# insect monitoring locations.
#
# This allows us to evaluate whether plant-based information is associated
# with observed insect richness and the resulting pollination indicators.

plant_richness_values <- terra::extract(
  plant_expected_richness_50m,
  terra::vect(observed_richness)
)


# ------------------------------------------------------------------------------
# 3.3 Assign nature regions to monitoring locations
# ------------------------------------------------------------------------------

# Spatially join each monitoring location to its corresponding nature region.
#
# Only the nature_region variable is retained because this is the regional
# grouping required for the calibration-effect analysis.

regions_eval_locs <- sf::st_join(
  observed_richness,
  
  nature_regions %>%
    dplyr::select(nature_region),
  
  join = sf::st_within
)


# ==============================================================================
# 4. Construct the insect indicator evaluation dataset
# ==============================================================================

# Combine:
#
#   - indicator predictions,
#   - expected plant richness,
#   - observed insect richness, and
#   - nature-region information.
#
# The resulting object, `rt`, is the main analysis dataset used throughout
# the remainder of the script.

rt <- indicator_values %>%
  
  # Add expected plant richness.
  bind_cols(
    plant_richness_values %>%
      select(-ID)
  ) %>%
  
  # Add monitoring-location information and nature region.
  bind_cols(
    regions_eval_locs
  ) %>%
  
  # Retain only variables required for the evaluation.
  select(
    
    # Expected pollinator richness.
    richness,
    
    # Calibrated expected pollinator richness.
    calibrated_richness,
    
    # Expected pollinator richness without interaction.
    richness_without_interaction,
    
    # Regional indicator.
    regional_clamped,
    
    # Calibrated regional indicator.
    calibrated_regional_clamped,
    
    # Regional indicator without interaction.
    regional_clamped_without_interaction,
    
    # Global indicator.
    global_clamped,
    
    # Calibrated global indicator.
    calibrated_global_clamped,
    
    # Global indicator without interaction.
    global_clamped_without_interaction,
    
    # Observed insect richness.
    insct_r,
    
    # Expected plant richness.
    plant_expected_richness,
    
    # Nature region.
    nature_region
  ) %>%
  
  # Remove locations containing missing values.
  na.omit()


# ==============================================================================
# 5. Prepare data for richness relationship plots
# ==============================================================================

# The four scatter plots evaluate the following relationships:
#
#   A. Pollination potential vs expected plant richness
#   B. Pollination potential vs observed insect richness
#   C. Scaled indicator vs observed insect richness
#   D. Pollination potential vs expected pollinator richness
#
# In each case, the variable being evaluated is plotted against the relevant
# ecological richness measure.

plot_data <- rt %>%
  
  # Geometry is not required for these statistical/graphical comparisons.
  sf::st_drop_geometry() %>%
  
  select(
    plant_expected_richness,
    calibrated_richness,
    insct_r,
    calibrated_regional_clamped,
    richness
  )


# ==============================================================================
# 6. Calculate Pearson correlations
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
  dat <- data %>%
    filter(
      !is.na(.data[[x_var]]),
      !is.na(.data[[y_var]])
    )
  
  
  # Calculate Pearson's correlation coefficient.
  cor_value <- cor(
    dat[[x_var]],
    dat[[y_var]],
    method = "pearson"
  )
  
  
  # Format the correlation for display on the figure.
  paste0(
    "Pearson r = ",
    sprintf("%.2f", cor_value)
  )
}


# ==============================================================================
# 7. Function for creating standardised scatter plots
# ==============================================================================

# ------------------------------------------------------------------------------
# Function: make_scatter()
#
# Creates a standardised scatter plot containing:
#
#   - individual observations,
#   - a linear regression line,
#   - its 95% confidence interval,
#   - Pearson's correlation coefficient, and
#   - a panel label.
#
# Arguments:
#
#   data        - plotting data
#   x_var       - variable plotted on the x-axis
#   y_var       - variable plotted on the y-axis
#   x_label     - x-axis label
#   y_label     - y-axis label
#   panel_label - panel identifier (A, B, C, D)
# ------------------------------------------------------------------------------

make_scatter <- function(
    data,
    x_var,
    y_var,
    x_label,
    y_label,
    panel_label
) {
  
  # Calculate Pearson's correlation for this relationship.
  cor_label <- get_correlation(
    data = data,
    x_var = x_var,
    y_var = y_var
  )
  
  
  ggplot(
    data,
    aes(
      x = .data[[x_var]],
      y = .data[[y_var]]
    )
  ) +
    
    # --------------------------------------------------------------------------
  # Monitoring locations
  # --------------------------------------------------------------------------
  geom_point(
    alpha = 0.6,
    size = 2
  ) +
    
    # --------------------------------------------------------------------------
  # Linear relationship
  #
  # The shaded region represents the 95% confidence interval around the
  # estimated regression line.
  # --------------------------------------------------------------------------
  geom_smooth(
    method = "lm",
    se = TRUE,
    linewidth = 0.8
  ) +
    
    # --------------------------------------------------------------------------
  # Pearson correlation
  #
  # Position the correlation coefficient in the upper-right corner.
  # --------------------------------------------------------------------------
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = cor_label,
    hjust = 1.1,
    vjust = 1.5,
    fontface = "bold",
    size = 4
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
  # Axis labels
  # --------------------------------------------------------------------------
  labs(
    x = x_label,
    y = y_label
  ) +
    
    # --------------------------------------------------------------------------
  # Plot theme
  # --------------------------------------------------------------------------
  theme_classic(
    base_size = 12
  ) +
    
    theme(
      
      axis.title = element_text(
        face = "bold"
      ),
      
      axis.text = element_text(
        colour = "black"
      ),
      
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        linewidth = 0.5
      )
    )
}


# ==============================================================================
# 8. Generate the four richness relationship plots
# ==============================================================================

# ------------------------------------------------------------------------------
# Panel A
#
# Question:
#   Is pollination potential associated with expected plant richness?
# ------------------------------------------------------------------------------

pA <- make_scatter(
  data = plot_data,
  x_var = "calibrated_richness",
  y_var = "plant_expected_richness",
  x_label = "Pollinator potential",
  y_label = "Expected plant richness",
  panel_label = "A"
)


# ------------------------------------------------------------------------------
# Panel B
#
# Question:
#   Is pollination potential associated with observed insect richness?
# ------------------------------------------------------------------------------

pB <- make_scatter(
  data = plot_data,
  x_var = "calibrated_richness",
  y_var = "insct_r",
  x_label = "Pollinator potential",
  y_label = "Observed insect richness",
  panel_label = "A"
)


# ------------------------------------------------------------------------------
# Panel C
#
# Question:
#   Is the calibrated regional indicator associated with observed insect
#   richness?
# ------------------------------------------------------------------------------

pC <- make_scatter(
  data = plot_data,
  x_var = "calibrated_regional_clamped",
  y_var = "insct_r",
  x_label = "Scaled indicator",
  y_label = " ",
  panel_label = "B"
)


# ------------------------------------------------------------------------------
# Panel D
#
# Question:
#   How strongly does pollination potential track the underlying expected
#   pollinator richness surface?
# ------------------------------------------------------------------------------

pD <- make_scatter(
  data = plot_data,
  x_var = "richness",
  y_var = "insct_r",
  x_label = "Expected richness",
  y_label = " ",
  panel_label = "C"
)


# ==============================================================================
# 9. Combine the four plots
# ==============================================================================

# Arrange the four panels in a 2 × 2 layout.
# p_combined <- (
#   pA | pB
# ) /
#   (
#     pC | pD
#   )

p_combined <- (
pB | pC | pD
)


# Display the combined figure.
p_combined


# ==============================================================================
# 10. Save richness relationship figure
# ==============================================================================

# ------------------------------------------------------------------------------
# PDF
# ------------------------------------------------------------------------------

ggsave(
  filename = file.path(
    plotFolder,
    "richness_calibration_relationships.pdf"
  ),
  plot = p_combined,
  width = 12,
  height = 10,
  units = "in"
)


# ------------------------------------------------------------------------------
# High-resolution PNG
# ------------------------------------------------------------------------------

ggsave(
  filename = file.path(
    plotFolder,
    "richness_calibration_relationships.png"
  ),
  plot = p_combined,
  width = 12,
  height = 10,
  units = "in",
  dpi = 600
)


# ==============================================================================
# 11. Evaluate the effect of calibration
# ==============================================================================

# The second part of the analysis evaluates how much calibration changes
# the regional and global indicators.
#
# Difference is defined as:
#
#     uncalibrated indicator - calibrated indicator
#
# Therefore:
#
#   positive difference  -> calibration decreases the indicator
#   negative difference  -> calibration increases the indicator
#   difference = 0       -> no effect of calibration


# ------------------------------------------------------------------------------
# 11.1 Calculate calibration-induced differences
# ------------------------------------------------------------------------------

indicator_difference <- rt %>%
  
  transmute(
    nature_region,
    
    # --------------------------------------------------------------------------
    # Regional indicator
    # --------------------------------------------------------------------------
    regional_difference =
      regional_clamped -
      calibrated_regional_clamped,
    
    # --------------------------------------------------------------------------
    # Global indicator
    # --------------------------------------------------------------------------
    global_difference =
      global_clamped -
      calibrated_global_clamped
  ) %>%
  
  # Convert regional and global differences into long format.
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
    
    # Define the plotting/model order.
    indicator_type = factor(
      indicator_type,
      levels = c(
        "Regional",
        "National"
      )
    )
  )


# ==============================================================================
# 12. Statistical tests of calibration effects
# ==============================================================================

# For every indicator type × nature-region combination, calculate:
#
#   - sample size,
#   - mean calibration-induced difference,
#   - standard deviation,
#   - one-sample t-test against zero.
#
# The null hypothesis is:
#
#     H0: mean calibration difference = 0
#
# P-values are adjusted using the Benjamini–Hochberg procedure to control the
# false discovery rate across the multiple regional comparisons.

calibration_tests <- indicator_difference %>%
  
  group_by(
    indicator_type,
    nature_region
  ) %>%
  
  summarise(
    
    # Number of observations.
    n = sum(
      !is.na(difference)
    ),
    
    # Mean difference between uncalibrated and calibrated indicators.
    mean_difference = mean(
      difference,
      na.rm = TRUE
    ),
    
    # Standard deviation of the calibration-induced difference.
    sd_difference = sd(
      difference,
      na.rm = TRUE
    ),
    
    # Test whether the mean difference differs from zero.
    p_value = t.test(
      difference,
      mu = 0
    )$p.value,
    
    .groups = "drop"
  ) %>%
  
  # Adjust p-values for multiple comparisons.
  mutate(
    p_adjusted = p.adjust(
      p_value,
      method = "BH"
    ),
    
    # Flag statistically significant calibration effects.
    significant = p_adjusted < 0.05
  )


# ==============================================================================
# 13. Save calibration test results
# ==============================================================================

write_csv(
  calibration_tests,
  file.path(
    plotFolder,
    "insect_calibration_difference_tests.csv"
  )
)


# ==============================================================================
# 14. Model calibration effects
# ==============================================================================

# Fit a linear model to test whether calibration effects differ according to:
#
#   1. indicator type (Regional vs Global),
#   2. nature region, and
#   3. their interaction.
#
# Model:
#
#   difference ~ indicator_type × nature_region
#
# The interaction tests whether the difference between Regional and Global
# calibration effects changes among nature regions.

difference_model <- lm(
  difference ~
    indicator_type *
    factor(nature_region),
  data = indicator_difference
)


# ==============================================================================
# 15. Calculate estimated marginal means
# ==============================================================================

# Estimate the model-adjusted mean calibration effect for each
# indicator type × nature-region combination.

emm <- emmeans(
  difference_model,
  ~ indicator_type | nature_region
)


# ==============================================================================
# 16. Prepare estimated effects for plotting
# ==============================================================================

emm_df <- as.data.frame(emm) %>%
  
  mutate(
    
    # A calibration effect is considered significant when its confidence
    # interval does not include zero.
    significant =
      lower.CL > 0 |
      upper.CL < 0,
    
    # Create a readable significance label.
    effect_label = if_else(
      significant,
      "Significant",
      "Not significant"
    )
  )


# Inspect estimated marginal means.
emm_df


# ==============================================================================
# 17. Plot estimated calibration effects
# ==============================================================================

p_effects <- ggplot(
  emm_df,
  aes(
    x = factor(nature_region),
    y = emmean,
    colour = indicator_type
  )
) +
  
  # Zero represents no change due to calibration.
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  
  # 95% confidence intervals around the estimated marginal means.
  geom_errorbar(
    aes(
      ymin = lower.CL,
      ymax = upper.CL
    ),
    position = position_dodge(
      width = 0.6
    ),
    width = 0.15,
    linewidth = 0.7
  ) +
  
  # Estimated mean calibration effects.
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
    axis.title = element_text(
      face = "bold"
    ),
    legend.title = element_text(
      face = "bold"
    )
  )


# Display the calibration-effect plot.
p_effects


# ==============================================================================
# 18. Save calibration-effect figure
# ==============================================================================

# ------------------------------------------------------------------------------
# PDF
# ------------------------------------------------------------------------------

ggsave(
  filename = file.path(
    plotFolder,
    "insect_calibration_difference_effect.pdf"
  ),
  plot = p_effects,
  width = 10,
  height = 7,
  units = "in",
  device = "pdf"
)


# ------------------------------------------------------------------------------
# High-resolution PNG
# ------------------------------------------------------------------------------

ggsave(
  filename = file.path(
    plotFolder,
    "insect_calibration_difference_effect.png"
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

