# Estimate diversity indices based on their taxon: either Bees, Butterflies or Hoverflies
# The function also takes into consideration whether you want to return the diversity estimate of a particular pollinator species
# If pollinatorSpeciesNames = NULL, then the function considers all the pollinator species within a zone or taxon.
# speciesForTaxon is the dataframe that contains the plant species prediction results. This is formatted from the fitted ISDMS to the ASO plant species
# taxa is the taxon we are interested in returning results for. It takes the values of either "allTaxa", "Butterflies", "Bees", "Hoverflies" 
# predRast is the prediction raster created in the masterScript
# silent indicates whether we want to print the species names or not



estimateTaxonDiversityWithoutInterraction <- function(pollinatorSpeciesNames = NULL, #Names of the pollinator species (should be the simplescientific name)
                                   speciesForTaxon, # The data with the plant species results
                                   taxa, 
                                   predRast, # The prediction raster
                                   silent = FALSE){
  
  if(!taxa %in% c("allTaxa", "Butterflies", "Bees", "Hoverflies" )) stop("Taxa must either be 'allTaxa', 'Butterflies', 'Bees', 'Hoverflies ")
  # Filter out the pollinator species we want
  if(!taxa %in% "allTaxa"){
    dat <- speciesForTaxon%>%
      dplyr::filter(Taxon %in% taxa)%>%
      na.omit()%>%
      filter(genus %in% colnames(interractionMatrix))
  } else if (taxa == "allTaxa"){
    
    dat <- speciesForTaxon%>%
      dplyr::filter(Taxon %in% c("Bees", "Butterflies", "Hoverflies"))%>%
      na.omit()  %>%
      filter(genus %in% colnames(interractionMatrix))
  }
  
  returnAllPrediction <- is.null(pollinatorSpeciesNames )
  
  if(is.null(pollinatorSpeciesNames )){
    pollinatorSpeciesNames <- dat$simpleScientificName
  } else {
    dat <- dat %>%
      dplyr::filter(simpleScientificName %in% pollinatorSpeciesNames)
  }
  
  if(!silent) print(pollinatorSpeciesNames)
  
  
  #indSpeciesDiversity <- lapply(seq_along(pollinatorSpeciesNames), function(x){
  sp <- pollinatorSpeciesNames#[x]
  genera <- dat$genus
  
  
  #get unique taxa
  focalTaxa <- dat$Taxon
  #alltaxa <- c("alpine", "boreal", "temperate")
  zones <- unique(dat$zone)
  
    #for now, we ensure that dat has only one pollinator group
    #if(length(unique(dat$Taxon))>1) stop("Select one zone only")
    
    # Load the intensity of the species from the zone
    spPred <- list()
    for(j in 1:length(sp)){
      print(paste("Estimating diversity indices for ", sp[j]))
      path <- paste(c(dataFolder, "modelOutputs", "pollinators", focalTaxa[j], sp[j]), collapse = "/")
      if(!file.exists(path)) next
      results <- readRDS(file.path(path, paste0("predictions", ".rds")))
      results$intensity <- round(exp(results$mean), 1)*10
      
      # Transform the intensity to the original scale
      diversity <-  results
      
      
      # define geometries to combine with prediction 
      geometries <- xyFromCell(diversity, seq(ncell(diversity))) %>% 
        as.data.frame() %>% 
        st_as_sf(coords = c("x", "y"), crs = newCrs)
      
      # Create the prediction data 
      PredictionData <- diversity %>% 
        as.data.frame(na.rm = FALSE)%>%
        bind_cols(geometries) %>% 
        tidyr::drop_na() %>%
        st_sf()%>%
        st_transform(., newCrs)
      
      PredictionData <- sf::st_intersection(PredictionData, 
                                            regionGeometry)
      
      # PredictionData <- PredictionData %>%
      #   mutate(richness = vegan::specnumber(across(all_of(plantSpecies))),
      #          simpson = vegan::diversity(across(all_of(plantSpecies)), "invsimpson"),
      #          shannon= vegan::diversity(across(all_of(plantSpecies))))
      #mutate(averageIntensity = sum(across(plantSpecies), na.rm = TRUE)/length(plantSpecies))%>%
      # reduce(cbind) %>% 
      
      #%>% # transform the prediction data to the units used in model runs
      #filter(rowSums(is.na(.)) != (ncol(.)-1)) # Now we take all the rows with sum of NAs equal to the number of columns of the dataframe minus the 
      
      #names(PredictionData) <- c(plantSpecies, "geometry")#, "speciesNumber", "simpson", "shannon")
      
      PredictionData <- PredictionData%>%
        rename(averageIntensity = intensity)%>%
        select(all_of(c("averageIntensity","geometry")))
      
      names(PredictionData) <- c(sp[j], "geometry")
      spPred[[j]] <- PredictionData #rasterize(PredictionData, predRast, names(PredictionData)[-length(names(PredictionData))])
      
      rm(list = c("PredictionData", "diversity", "results", "geometries"))
     # gc()
      gc()
    }
    
    sf_combined <- spPred[!sapply(spPred, is.null)]%>%
      do.call(cbind, .)%>% select(-starts_with("geometry.")) %>%
      dplyr::mutate(richness = vegan::specnumber(across(.cols = !starts_with("geometry"))),
                    simpson = vegan::diversity(across(.cols = !starts_with("geometry")), "invsimpson"),
                    shannon= vegan::diversity(across(.cols = !starts_with("geometry"))))%>%
      dplyr::select(richness, shannon, simpson)
  

    
    # path <- paste0(dataFolder, "/modelOutputs/diversity/", sp)
    # # spPred <- rasterize(spPred, predRast, names(spPred))
    # if(!file.exists(path)){
    #   dir.create(path)
    # }
    
    #print("Saving Results")
    
    #st_write(sf_combined, file.path(path, paste0("predsWithoutInterraction.shp")))
    #spPred <- rasterize(sf_combined , predRast, c("richness", "shannon", "simpson"))
  
    return(sf_combined)
  
}  