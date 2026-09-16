library(terra)
library(sf)
library(dplyr)
library(tidyr)
library(ggplot2)
library(grid)

## ========================================================
## 1. Extract indicator layers
## ========================================================

allRasts <- terra::rast(file.path(
  plotFolder,
  "all_indicator_raster.tif"
))

regional_clamped <-
  allRasts[[
    "regional_clamped"
  ]]

calibrated_regional_clamped <-
  allRasts[[
    "calibrated_regional_clamped"
  ]]

regional_clamped_without_interaction <-
  allRasts[[
    "regional_clamped_without_interaction"
  ]]

global_clamped <-
  allRasts[[
    "global_clamped"
  ]]

calibrated_global_clamped <-
  allRasts[[
    "calibrated_global_clamped"
  ]]

global_clamped_without_interaction <-
  allRasts[[
    "global_clamped_without_interaction"
  ]]

region_id <-
  allRasts[[
    "region_id"
  ]]


## ========================================================
## 2. Convert REGIONAL indicators to data frame
## ========================================================

regional_df <- terra::as.data.frame(
  c(
    regional_clamped,
    calibrated_regional_clamped,
    regional_clamped_without_interaction
  ),
  xy = TRUE,
  na.rm = TRUE
)

names(regional_df)[3:5] <- c(
  "Regional indicator",
  "Calibrated regional indicator",
  "Regional indicator without interaction"
)


## ========================================================
## 3. Convert GLOBAL indicators to data frame
## ========================================================

global_df <- terra::as.data.frame(
  c(
    global_clamped,
    calibrated_global_clamped,
    global_clamped_without_interaction
  ),
  xy = TRUE,
  na.rm = TRUE
)

names(global_df)[3:5] <- c(
  "Global indicator",
  "Calibrated global indicator",
  "Global indicator without interaction"
)


## ========================================================
## 4. Combine all indicators
## ========================================================

indicator_df <- bind_rows(
  regional_df,
  global_df
)


## ========================================================
## 5. Convert indicators to long format
## ========================================================

indicator_df <- indicator_df %>%
  pivot_longer(
    cols = -c(x, y),
    names_to = "indicator_type",
    values_to = "indicator"
  )


## ========================================================
## 6. Set indicator order
## ========================================================

indicator_df$indicator_type <- factor(
  indicator_df$indicator_type,
  levels = c(
    "Regional indicator",
    "Calibrated regional indicator",
    "Regional indicator without interaction",
    "Global indicator",
    "Calibrated global indicator",
    "Global indicator without interaction"
  )
)


## ========================================================
## 7. Add region ID
## ========================================================

region_df <- terra::as.data.frame(
  region_id,
  xy = TRUE,
  na.rm = FALSE
)

names(region_df)[3] <- "region_id"


indicator_df <- indicator_df %>%
  left_join(
    region_df,
    by = c("x", "y")
  )


## ========================================================
## 8. Add nature-region names
## ========================================================

indicator_df <- indicator_df %>%
  left_join(
    indicator_results$region_lookup %>%
      select(
        region_id,
        nature_region
      ),
    by = "region_id"
  )


## ========================================================
## 9. Calculate mean indicator for each region
## ========================================================

mean_by_region <- indicator_df %>%
  group_by(
    indicator_type,
    region_id,
    nature_region
  ) %>%
  summarise(
    mean_indicator = mean(
      indicator,
      na.rm = TRUE
    ),
    .groups = "drop"
  )


## ========================================================
## 10. Create regional mean table
## ========================================================

region_mean_labels <- mean_by_region %>%
  filter(
    !is.na(nature_region)
  ) %>%
  select(
    nature_region,
    indicator_type,
    mean_indicator
  ) %>%
  pivot_wider(
    names_from = indicator_type,
    values_from = mean_indicator
  ) %>%
  mutate(
    region_label = paste0(
      nature_region,
      "\n",
      
      "R: ",
      sprintf(
        "%.2f",
        `Regional indicator`
      ),
      
      "   R-C: ",
      sprintf(
        "%.2f",
        `Calibrated regional indicator`
      ),
      
      "   R-NI: ",
      sprintf(
        "%.2f",
        `Regional indicator without interaction`
      ),
      
      "\n",
      
      "G: ",
      sprintf(
        "%.2f",
        `Global indicator`
      ),
      
      "   G-C: ",
      sprintf(
        "%.2f",
        `Calibrated global indicator`
      ),
      
      "   G-NI: ",
      sprintf(
        "%.2f",
        `Global indicator without interaction`
      )
    )
  )


## ========================================================
## 11. Named vector for region legend
## ========================================================

region_labels <- region_mean_labels$region_label

names(region_labels) <-
  region_mean_labels$nature_region


## ========================================================
## 12. Convert regions to polygons
## ========================================================

region_poly <- terra::as.polygons(
  region_id,
  dissolve = TRUE,
  na.rm = TRUE
)

region_poly <- sf::st_as_sf(
  region_poly
)

names(region_poly)[
  names(region_poly) == "region_id"
] <- "region_id"


## ========================================================
## 13. Add region names
## ========================================================

region_poly <- region_poly %>%
  left_join(
    indicator_results$region_lookup %>%
      select(
        region_id,
        nature_region
      ),
    by = "region_id"
  )


## ========================================================
## 14. Plot
## ========================================================

p <- ggplot() +
  
  ## ======================================================
## Indicator raster
## ======================================================

geom_raster(
  data = indicator_df,
  aes(
    x = x,
    y = y,
    fill = indicator
  )
) +
  
  
  ## ======================================================
## Nature-region boundaries
## ======================================================

geom_sf(
  data = region_poly,
  aes(
    colour = nature_region
  ),
  fill = NA,
  linewidth = 0.8
) +
  
  
  ## ======================================================
## Scaled indicator colour scale
## ======================================================

scale_fill_viridis_c(
  name = "Scaled indicator",
  limits = c(0, 1),
  breaks = seq(
    0,
    1,
    0.25
  ),
  na.value = NA,
  guide = guide_colourbar(
    direction = "horizontal",
    title.position = "top",
    title.hjust = 0.5,
    barwidth = unit(4, "cm"),
    barheight = unit(0.4, "cm")
  )
) +
  
  
  ## ======================================================
## Nature-region legend
## ======================================================

scale_colour_discrete(
  name = paste0(
    "Nature region\n",
    "R = Regional | R-C = Calibrated regional | ",
    "R-NI = Regional without interaction\n",
    "G = Global | G-C = Calibrated global | ",
    "G-NI = Global without interaction"
  ),
  breaks = region_mean_labels$nature_region,
  labels = region_labels,
  guide = guide_legend(
    nrow = 3,
    byrow = TRUE,
    title.position = "top",
    title.hjust = 0,
    keywidth = unit(0.5, "cm"),
    keyheight = unit(0.5, "cm")
  )
) +
  
  
  ## ======================================================
## Facets
## ======================================================

facet_wrap(
  ~ indicator_type,
  ncol = 2,
  labeller = as_labeller(
    c(
      
      "Regional indicator" =
        "Regional indicator",
      
      "Calibrated regional indicator" =
        "Calibrated regional indicator",
      
      "Regional indicator without interaction" =
        "Regional indicator\n(no interaction)",
      
      "Global indicator" =
        "Global indicator",
      
      "Calibrated global indicator" =
        "Calibrated global indicator",
      
      "Global indicator without interaction" =
        "Global indicator\n(no interaction)"
    )
  )
) +
  
  
  ## ======================================================
## Coordinates
## ======================================================

coord_sf(
  expand = FALSE
) +
  
  
  ## ======================================================
## Axis labels
## ======================================================

labs(
  x = "Easting",
  y = "Northing"
) +
  
  
  ## ======================================================
## Theme
## ======================================================

theme_bw() +
  
  theme(
    
    panel.grid = element_blank(),
    
    ## ----------------------------------------------------
    ## Legends
    ## ----------------------------------------------------
    
    legend.position = "right",
    
    legend.direction = "horizontal",
    
    legend.box = "vertical",
    
    legend.box.just = "left",
    
    legend.spacing.y = unit(
      0.3,
      "cm"
    ),
    
    legend.text = element_text(
      size = 8
    ),
    
    legend.title = element_text(
      face = "bold",
      size = 9
    ),
    
    ## ----------------------------------------------------
    ## Facet strips
    ## ----------------------------------------------------
    
    strip.text = element_text(
      face = "bold",
      size = 10,
      lineheight = 1.1
    ),
    
    strip.background = element_rect(
      fill = "grey90",
      colour = "black"
    ),
    
    ## ----------------------------------------------------
    ## Panel spacing
    ## ----------------------------------------------------
    
    panel.spacing = unit(
      0.8,
      "lines"
    )
  )


p
