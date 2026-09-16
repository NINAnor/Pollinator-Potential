# Create a function that estimates covariate values for regions outside the defined
# region, due to the mesh extensions
covariateMissingValues <- function(where, index, layerName, parameters) {
  # Extract the values
  v <- eval_spatial(parameters[[index]], where, layer = layerName)
  # Fill in missing values
  if (any(is.na(v))) {
    v <- bru_fill_missing(parameters[[index]], where, v)
  }
  return(v)
}
