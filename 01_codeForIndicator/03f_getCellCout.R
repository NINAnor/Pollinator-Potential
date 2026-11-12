################################################################################
# Counts of valid observations per dataset and region, taken from the results
# files written by 03b (all open-lowland cells), 03c (ASO), 03d (ANO) and
# 03g (ASO + ANO). Run those scripts first.
#
# For 03b, n_sampled_cells is the number of valid raster cells in the region
# (all cells are used, not a subsample). For the survey datasets it is the
# number of survey locations with a valid indicator value.
################################################################################

output_file <- file.path(plotFolder, "allCellCounts.csv")

result_files <- c(
  "All open lowland boundary" = "region_results.rds",
  "ASO" = "aso_region_results.rds",
  "ANO" = "ano_region_results.rds",
  "ASO + ANO" = "aso_ano_region_results.rds"
)

missing_files <- result_files[!file.exists(file.path(plotFolder, result_files))]
if (length(missing_files)) {
  stop("Missing results file(s); run 03b/03c/03d/03g first: ",
       paste(missing_files, collapse = ", "))
}

allCellCounts <- purrr::imap_dfr(result_files, function(f, dataset) {
  readRDS(file.path(plotFolder, f)) |>
    dplyr::distinct(nature_region, n_sampled_cells) |>
    dplyr::transmute(
      nature_region,
      n_non_NA = n_sampled_cells,
      dataset = dataset
    )
})

write.csv(allCellCounts, output_file, row.names = FALSE)
