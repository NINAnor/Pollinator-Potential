source("01_codeForEnvAgency/03f_getCellCout.R")

region_map <- readRDS(file.path(
      plotFolder,
      "region_results.rds"
    )) %>%
  dplyr::mutate(dataset = "All open lowland boundary")


aso_map <- readRDS(file.path(
      plotFolder,
      "aso_region_results.rds"
    )) %>%
  dplyr::mutate(dataset = "ASO")

ano_map <- readRDS(file.path(
      plotFolder,
      "ano_region_results.rds"
    )) %>%
  dplyr::mutate(dataset = "ANO")


aso_ano_map <- readRDS(file.path(
      plotFolder,
      "aso_ano_region_results.rds"
    )) %>%
  dplyr::mutate(dataset = "ASO + ANO")

all_results <- bind_rows(region_map,
aso_map,
ano_map,
aso_ano_map) %>%
  dplyr::filter(indicator %in% "Calibrated regional indicator") %>%
  dplyr::left_join(., allCellCounts, by =c("nature_region", "dataset")) %>%
  dplyr::select(dataset, nature_region, bootstrap_mean,
  bootstrap_sd, lower_quantile, upper_quantile, n_non_NA) %>%
  dplyr::mutate(nature_region = if_else(nature_region %in% "national", "National", nature_region))

colnames(all_results) <- c("Dataset", "Region", "Mean", "Standard deviation", "Lower quantile", "Upper quantile", "Cell counts")

# Sasve workbook


# ============================================================
# Save all_results as an Excel worksheet
# ============================================================

# Install openxlsx if it is not already installed
if (!requireNamespace("openxlsx", quietly = TRUE)) {
  install.packages("openxlsx")
}

library(openxlsx)

# Define output file
output_file <- file.path(
  plotFolder,
  "all_results.xlsx"
)

# Create workbook
wb <- createWorkbook()

# Add worksheet
addWorksheet(
  wb,
  sheetName = "all_results"
)

# Write results to worksheet
writeData(
  wb,
  sheet = "all_results",
  x = all_results,
  withFilter = TRUE
)

# Format header
header_style <- createStyle(
  textDecoration = "bold",
  fgFill = "#D9EAF7",
  border = "Bottom"
)

addStyle(
  wb,
  sheet = "all_results",
  style = header_style,
  rows = 1,
  cols = 1:ncol(all_results),
  gridExpand = TRUE
)

# Format numeric columns to five decimal places
numeric_style <- createStyle(
  numFmt = "0.00000"
)

addStyle(
  wb,
  sheet = "all_results",
  style = numeric_style,
  rows = 2:(nrow(all_results) + 1),
  cols = 3:6,
  gridExpand = TRUE
)

# Freeze the header row
freezePane(
  wb,
  sheet = "all_results",
  firstRow = TRUE
)

# Automatically adjust column widths
setColWidths(
  wb,
  sheet = "all_results",
  cols = 1:ncol(all_results),
  widths = "auto"
)

# Save workbook
saveWorkbook(
  wb,
  file = output_file,
  overwrite = TRUE
)

# Confirm output
message("Excel file saved to: ", output_file)

# Verify that the file exists
file.exists(output_file)

