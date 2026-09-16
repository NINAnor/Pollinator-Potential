# Function to estimate indices
estimateIndicesFnx <- function(predGrid, 
                               regionGeometry,
                               indicator){
  
  # Extract the values for each location
  vals <- terra::extract(predGrid, vect(regionGeometry))%>% #extract the covariates within the prediction grid
    dplyr::select(bio1, Latitude, landCover) %>% #select only the geo-climatic variables
    cbind(., terra::extract(indicator, vect(regionGeometry), xy = TRUE))%>% #extract the diversity values within the prediction grid
    na.omit() 
  
  # Estimate the maximum diversity at each zone and save them
  maxVals <- vals %>%
    dplyr::mutate(climaticZones = ifelse(Latitude > 0, 1,2))%>%
    dplyr::group_by(climaticZones)%>%
    dplyr::summarise(maxReference = max(mean),
                     minReference = min(mean))
  
  #scale the diversity by their geo-climatic maximum values
  indicators <- vals %>%
    select(Latitude, x, y, mean)%>%
    rename(all = mean)%>%
    dplyr::mutate(climaticZones = ifelse(Latitude > 0, 1, 2))%>%
    dplyr::mutate(allReference = ifelse(climaticZones == 1, all/maxVals$maxReference[1], (all/maxVals$maxReference[2])))%>%
    sf::st_as_sf(., coords = c("x", "y"), crs = newCrs)%>% # convert to sf object
    dplyr::select(climaticZones, 
                  all,
                  allReference) 
  
  # convert the results to a raster to plot
  spPred <- rasterize(indicators, 
                      predRast, 
                      c("climaticZones", 
                        "all",
                        "allReference"))%>%
    terra::crop(., regionGeometry)
  
  # return the results needed
  ret <- list(maxVals = maxVals, 
              indicator = indicators, 
              spPred = spPred)
  
  return(ret)
}