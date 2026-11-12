library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(openxlsx)


observed_data <- reference_polygon_values

# ============================================================
# 1. Always calculate regional parameters
# ============================================================
global_min <- min(observed_data$calibrated_reference_richness, na.rm = TRUE)

regional_parameters <- observed_data %>%
  group_by(nature_region) %>%
  summarise(
    x_min = 1,
    x_max = max(calibrated_reference_richness, na.rm = TRUE),
    n_observed = sum(is.finite(calibrated_reference_richness)),
   X_0 = 0,
    #X_0 = global_min,
    X_60 = median(calibrated_reference_richness, na.rm = TRUE),
    X_100 = max(calibrated_reference_richness, na.rm = TRUE),
    #X_100 = quantile(calibrated_reference_richness, 0.9, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    mu = (log(x_min) + log(x_max)) / 2,
    sigma = (log(x_max) - log(x_min)) / (2 * 1.96),
    region_id = as.numeric(factor(nature_region))
  )

# These are always available, including when checks = FALSE
reference_values <- regional_parameters %>%
  select(nature_region, region_id, X_0, X_60, X_100)

# ============================================================
# 2. Run simulations and diagnostics only when requested
# ============================================================

if (checks) {

  n_replicates <- 1000
  n_density_points <- 512

  set.seed(123)

  # ----------------------------------------------------------
  # Simulate n_replicates datasets per region
  # ----------------------------------------------------------

  simulated_data <- purrr::pmap_dfr(
    regional_parameters %>%
      select(nature_region, n_observed, mu, sigma),
    function(nature_region, n_observed, mu, sigma) {

      tidyr::expand_grid(
        simulation_id = seq_len(n_replicates),
        observation_id = seq_len(n_observed)
      ) %>%
        mutate(
          nature_region = nature_region,
          richness = rlnorm(
            n(),
            meanlog = mu,
            sdlog = sigma
          )
        )
    }
  ) %>%
    select(nature_region, simulation_id, observation_id, richness)

  # ----------------------------------------------------------
  # Validate simulations
  # ----------------------------------------------------------

  simulation_check <- simulated_data %>%
    count(nature_region, simulation_id) %>%
    group_by(nature_region) %>%
    summarise(
      n_replicates = n_distinct(simulation_id),
      min_sample_size = min(n),
      max_sample_size = max(n),
      .groups = "drop"
    )

  print(simulation_check)

  # ----------------------------------------------------------
  # Prepare observed values
  # ----------------------------------------------------------

  observed_plot_data <- observed_data %>%
    filter(is.finite(calibrated_reference_richness)) %>%
    group_by(nature_region) %>%
    mutate(
      observation_id = row_number(),
      richness = calibrated_reference_richness
    ) %>%
    ungroup() %>%
    select(nature_region, observation_id, richness)

  # ----------------------------------------------------------
  # Estimate density for each replicate
  # ----------------------------------------------------------

  x_grid <- seq(
    min(c(observed_plot_data$richness, simulated_data$richness)),
    max(c(observed_plot_data$richness, simulated_data$richness)),
    length.out = n_density_points
  )

  replicate_densities <- simulated_data %>%
    group_by(nature_region, simulation_id) %>%
    group_modify(~ {
      d <- density(
        .x$richness,
        from = min(x_grid),
        to = max(x_grid),
        n = n_density_points
      )

      tibble(richness = d$x, density = d$y)
    }) %>%
    ungroup()

  mean_density <- replicate_densities %>%
    group_by(nature_region, richness) %>%
    summarise(
      density = mean(density),
      .groups = "drop"
    )

  # ----------------------------------------------------------
  # Density plot
  # ----------------------------------------------------------

  density_plot <- ggplot() +
    geom_line(
      data = replicate_densities,
      aes(x = richness, y = density, group = simulation_id),
      colour = "steelblue",
      alpha = 0.06,
      linewidth = 0.35
    ) +
    geom_line(
      data = mean_density,
      aes(x = richness, y = density),
      colour = "darkred",
      linewidth = 1.3
    ) +
    geom_density(
      data = observed_plot_data,
      aes(x = richness),
      colour = "darkgreen",
      linewidth = 1.2,
      linetype = "dashed"
    ) +
    facet_wrap(~ nature_region, scales = "free") +
    labs(
      title = "Observed and simulated richness distributions",
      subtitle = paste(n_replicates, "replicates per region"),
      x = "Calibrated reference richness",
      y = "Density",
      caption = paste(
        "Blue: individual simulated densities;",
        "red: mean simulated density;",
        "green dashed: observed density."
      )
    ) +
    theme_minimal(base_size = 13) +
    theme(
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold")
    )

  print(density_plot)

  # ----------------------------------------------------------
  # Prepare boxplot data
  # ----------------------------------------------------------

  simulated_median <- simulated_data %>%
    group_by(nature_region, observation_id) %>%
    summarise(
      richness = median(richness),
      .groups = "drop"
    ) %>%
    mutate(type = "Median of 1,000 simulations")

  simulated_first <- simulated_data %>%
    filter(simulation_id == 1) %>%
    transmute(
      nature_region,
      observation_id,
      richness,
      type = "Simulation 1"
    )

  boxplot_data <- bind_rows(
    observed_plot_data %>%
      mutate(type = "Observed"),
    simulated_median,
    simulated_first
  ) %>%
    mutate(
      type = factor(
        type,
        levels = c(
          "Observed",
          "Median of 1,000 simulations",
          "Simulation 1"
        )
      )
    )

  # ----------------------------------------------------------
  # Boxplot
  # ----------------------------------------------------------

  boxplot_plot <- ggplot(
    boxplot_data,
    aes(x = nature_region, y = richness, fill = type)
  ) +
    geom_boxplot(
      position = position_dodge(width = 0.85),
      width = 0.7,
      alpha = 0.8
    ) +
    labs(
      title = "Observed and simulated richness by region",
      x = "Nature region",
      y = "Calibrated reference richness",
      fill = "Dataset"
    ) +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold"))

  print(boxplot_plot)

}


# Regional parameters to save

regional_parameters_to_save <- observed_data %>%
  group_by(nature_region) %>%
  summarise(
    x_min = min(calibrated_reference_richness, na.rm = TRUE),
    x_max = max(calibrated_reference_richness, na.rm = TRUE),
    n_observed = sum(is.finite(calibrated_reference_richness)),

    X_0   = 0,
    X_50  = median(calibrated_reference_richness, na.rm = TRUE),
    X_90  = quantile(calibrated_reference_richness, 0.90, na.rm = TRUE),
    X_95  = quantile(calibrated_reference_richness, 0.95, na.rm = TRUE),
    X_99  = quantile(calibrated_reference_richness, 0.99, na.rm = TRUE),
    X_100 = max(calibrated_reference_richness, na.rm = TRUE),

    .groups = "drop"
  )



wb <- createWorkbook()

addWorksheet(wb, "Regional thresholds")

writeDataTable(
  wb,
  sheet = "Regional thresholds",
  x = regional_parameters_to_save,
  tableStyle = "TableStyleMedium2"
)

setColWidths(
  wb,
  sheet = "Regional thresholds",
  cols = 1:ncol(regional_parameters_to_save),
  widths = "auto"
)

freezePane(
  wb,
  sheet = "Regional thresholds",
  firstRow = TRUE
)

saveWorkbook(
  wb,
  file = file.path(plotFolder ,"regional_parameters.xlsx"),
  overwrite = TRUE
)
