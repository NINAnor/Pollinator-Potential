#' @title \emph{createPredictionData}: Creates a prediction field for the pipeline run to project projections into

#' @description We need to define the resolution and projection of the predictions that are outputted from the intSDM. This function creates that.
#'
#' @param resolution A numeric vector giving the size of the pixels which should be used in the prediction projection.
#' @param regionGeometry An sf multipolygon object giving the region that the projection should cover.
#' @param proj The projection's crs.
#' 
#' 
#' @return An sf object as a collection of points covering the specified region
#'





createPredictionData <- function(resolution, regionGeometry, proj = "Proj.4 +proj=utm +zone=32 +datum=WGS84 +units=km +no_defs +type=crs", environmentalCovariates = NULL) {
  if(is.null(environmentalCovariates)){
  Tron <- st_transform(regionGeometry, proj)
  predictionData <- st_make_grid(Tron, c(resolution[1], resolution[2]))
  intersectData <- sapply(st_intersects(predictionData, Tron), function(z) if (length(z)==0) FALSE else TRUE)
  predictionData <- predictionData[intersectData] #if need be, convert to points and select every 5th?
  predictionData <- st_as_sf(predictionData)#, data = data.frame())
  PredictionData <- predictionData %>% st_cast('POINT')
  PredictionData <- PredictionData[seq(1,length(PredictionData), 5)]
  st_geometry(PredictionData) <- 'geometry'
  } else {
    # Transform regionGeometry
    Tron <- st_transform(regionGeometry, proj)
    
    # Create a prediction raster
    predRast <- rast(ext(environmentalCovariates), 
                     res = c(resolution[1], resolution[2]), 
                     crs = proj)
    
    # Create a prediction grid using one of the environmental covariates
    predGrid <- project(environmentalCovariates$bio1, predRast)
    
    predGrid <- Tron %>%
      st_transform(., proj)%>%
      vect( )%>%
      mask(predGrid, .)
    
    # define geometries to combine with prediction 
    geometries <- xyFromCell(predGrid, seq(ncell(predGrid))) %>% 
      as.data.frame() %>% 
      st_as_sf(coords = c("x", "y"), crs = proj)
   
    # Create the prediction data 
    PredictionData <- predGrid %>% 
      as.data.frame(na.rm = FALSE)%>% 
      reduce(cbind) %>% 
      bind_cols(geometries) %>% 
      st_sf()%>%
      st_transform(., proj)%>% # transform the prediction data to the units used in model runs
      filter(rowSums(is.na(.)) != (ncol(.)-1))%>% # Now we take all the rows with sum of NAs equal to the number of columns of the dataframe minus the 
      pull(geometry)%>%
      st_as_sf()
    
    # Set the geometry of the prediction data
    st_geometry(PredictionData) <- "geometry"
  }
  return(PredictionData)
}
