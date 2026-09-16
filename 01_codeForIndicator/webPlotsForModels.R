# Interraction web
# Import and format the pollinator datasets
if(!file.exists(paste0(dataFolder, "/interractionsData/interractionMatrix.csv"))){
  FinalSelectionOfLepidopteraAndBeesAndHoverflies <- allPollinatorsDataset(dataFolder, 
                                                                           print = FALSE)
  Bio1Traits <- readr::read_csv("Data/pollinatorDataFolder/bioTraits.csv")%>%
    mutate(
      simpleScientificName = coalesce(
        #redList$species[match(acceptedScientificName, redList$GBIFName)],  # Match redList species
        str_extract(acceptedScientificName, "^[A-Za-z]+\\s+[a-z]+")        # Extract binomial name
      )
    )%>%
    mutate(genus = simpleScientificName)%>%
    separate(genus, into = c("genus", "other"), sep = " ")
  
  #NiNIndicators <- read.csv(paste0(dataFolder, "/interractionsData/Indikatorarter NiN eng.csv"), sep = ";")
  #NiNIndicatorsGenera <- do.call(rbind,lapply(unique(NiNIndicators$Genus),function(x)NiNIndicators[NiNIndicators$Genus %in% x,][1,]))
  #FinalSelectionOfLepidopteraAndBeesAndHoverflies$PlantGenusPollinator <- with(FinalSelectionOfLepidopteraAndBeesAndHoverflies, paste(PlantGenus,ValidSpeciesName))
  # OneInteractionPerPlantGenusPollinator[OneInteractionPerPlantGenusPollinator$PlantGenus %in% NiNIndicatorsGenera$Genus,]%>%
  #OneInteractionPerPlantGenusPollinator <- do.call(rbind,lapply(unique(FinalSelectionOfLepidopteraAndBeesAndHoverflies$PlantGenusPollinator),function(x)FinalSelectionOfLepidopteraAndBeesAndHoverflies[FinalSelectionOfLepidopteraAndBeesAndHoverflies$PlantGenusPollinator%in% x,][1,]))
  
  
  rr <- data.frame(scientificName = unique(FinalSelectionOfLepidopteraAndBeesAndHoverflies$ValidSpeciesName),
                   acceptedScientificName = sapply(unique(FinalSelectionOfLepidopteraAndBeesAndHoverflies$ValidSpeciesName), function(x) findGBIFName(x)))
  
  InteractionsInIdealMeadow <-  FinalSelectionOfLepidopteraAndBeesAndHoverflies%>%
    group_by(PlantGenus, ValidSpeciesName)%>%
    dplyr::distinct( )%>%
    ungroup() %>%
    mutate(genus = ValidSpeciesName)%>%
    separate(genus, into = c("genus", "other"), sep = " ")%>%
    select(PlantGenus, Taxon, genus)%>%
    dplyr::mutate(Taxon = ifelse(Taxon == "Butterfly", "Butterflies", Taxon))
  
  otherPollinatorSp <- unique(Bio1Traits$genus[!Bio1Traits$genus %in% InteractionsInIdealMeadow$genus])  
  otherPollinatorSp <- otherPollinatorSp[!is.na(otherPollinatorSp)]

  
  countries <- ne_countries(scale = "medium", returnclass = "sf")
  
  # Filter Europe
  scand <- countries[countries$name %in% c("Norway", "Sweden"), ]
  poly <- st_transform(scand, 4326)
  bbox <- st_bbox(poly)
  
  # Pollinator interactions (example: bees)
interactions <-  lapply(otherPollinatorSp, function(x){
  interactions <- rglobi::get_interactions_by_taxa(
    sourcetaxon = x,
    interactiontype = "pollinates",
    bbox = bbox
  )


 ret <-  interactions %>%
    select(source_taxon_name, target_taxon_name, study_citation, study_source_citation)
 
 return(ret)
  })

# Keep only those with at least one row
interactions1 <- interactions[sapply(interactions, nrow) > 0] %>%
  do.call("rbind", .) %>%
  dplyr::select(source_taxon_name, target_taxon_name) %>%
  separate(source_taxon_name, into = c("genus", "other"), sep = " ") %>%
  separate(target_taxon_name, into = c("PlantGenus", "otherSp"), sep = " ") %>%
  filter(PlantGenus %in% InteractionsInIdealMeadow$PlantGenus) %>%
  dplyr::left_join(., 
                   Bio1Traits %>% select(genus, Taxon),
                   by = "genus",
                   keep = FALSE) %>%
  dplyr::select(colnames(InteractionsInIdealMeadow))
  
allInterractions <- bind_rows(InteractionsInIdealMeadow,
                              interactions1)

  
  # InteractionsInIdealMeadowWithTraits <- Bio1Traits %>%
  #   select(genus, Taxon)%>%
  #   dplyr::full_join(allInterractions, 
  #                    .,
  #                    by = c("genus", "Taxon"))

  #dplyr::select("ValidSpeciesName", "Taxon", "acceptedScientificName", "PlantGenus")
  
  # tpr <- InteractionsInIdealMeadow %>%
  #   filter(Taxon == "Butterflies")
  # tpr%>%
  #   filter(genus %in% unique(Bio1Traits$genus))
  # 
  # table(tpr$genus)
  # tpr$genus[!is.na(tpr$PlantGenus)]  
  
  # FInd out the species in the aso data that match up to the specific genus
  if(!exists("asoDatasf")) asoDatasf <- readRDS(file = paste0(dataFolder, "/plantDataFolder/formattedData/asoData.RDS"))
  
  # Merge the asoData genus to the species data
  asoDatasf <- asoDatasf %>%
    mutate(
      simpleScientificName = coalesce(
        str_extract(acceptedScientificName, "^[A-Za-z]+\\s+[a-z]+")
      ),
      genus   = str_extract(simpleScientificName, "^[A-Za-z]+"),
      otherSp = str_extract(simpleScientificName, "(?<=\\s)[a-z]+")
    ) %>%
    st_drop_geometry() %>%
    distinct(simpleScientificName, genus)
  
  interractionData <- merge(allInterractions, 
                            as.data.frame(asoDatasf[, c("simpleScientificName", "genus")])[, 1:2],
                            by.x = "PlantGenus",
                            by.y = "genus")%>%
    group_by(genus, Taxon)%>%
    distinct()
  
  # save the interractions in ideal meadow
  # write.csv(InteractionsInIdealMeadowWithTraits, 
  #           paste0(dataFolder, "/interractionsData/InteractionsInIdealMeadowWithTraitsAlt.csv"))
  
  write.csv(interractionData, 
            paste0(dataFolder, "/interractionsData/interractionDataAlt.csv"))
  
  # Create the plant pollinator for network
  #PlantPollinatorsForNetwork <- InteractionsInIdealMeadowWithTraits[c("ValidSpeciesName","PlantGenus")]
  PlantPollinatorsForNetwork <- interractionData[c("genus","simpleScientificName")]
  names(PlantPollinatorsForNetwork)[1:2] <- c("higher","lower")
  PlantPollinatorsForNetwork$freq <- 1
  PlantPollinatorsForNetwork$webID <- 1
  dfForWeb <- PlantPollinatorsForNetwork[c("lower","higher","webID","freq")]
  
  # Create webs from dataframe
  WebReady <- bipartite::frame2webs(dfForWeb, varnames = c("lower", "higher", "webID", "freq"))
  
  # save for later use
  save(WebReady,
       file = paste0(dataFolder, "/interractionsData/webPlotAlt.RData"))
  
  if(modelPlots){
    plotweb(WebReady[[1]])
  }
  
  
  # Calculte the interraction matrix
  interractionMatrix <- as.data.frame(WebReady[[1]]) 
  interractionMatrix$species <- rownames(interractionMatrix)
  rownames(interractionMatrix) <- NULL
  # save the results
  write.csv(interractionMatrix, 
            paste0(dataFolder, "/interractionsData/interractionMatrixAlt.csv"))
  
} else {
  load(paste0(dataFolder,"/interractionsData/webPlotAlt.RData"))
  interractionMatrix <- readr::read_csv(paste0(dataFolder,"/interractionsData/interractionMatrixAlt.csv"))
  
  if(modelPlots){
    plotweb(WebReady[[1]])
  }
}

interractionProb <- interractionMatrix

# # Calculate the
# taxNames <- names(interractionMatrix[, 2:4])
interractionProb <- interractionMatrix%>%
  #                     as.matrix()%>%
  #                     #prop.table(., margin = 1)%>%
  #                     as.data.frame(., col.names = taxnames)%>%
  #                     dplyr::mutate_at(vars(taxNames), .funs = function(x) ifelse(x >0, 1, 0))%>%
  dplyr::mutate(species = interractionMatrix$species)%>%
  mutate(
    simpleScientificName = coalesce(       #redList$species[match(acceptedScientificName, redList$GBIFName)],  # Match redList species
      str_extract(species, "^[A-Za-z]+\\s+[a-z]+")        # Extract binomial name     
    ),
    # Replace space with underscore in simpleScientificName
    species = gsub("-", "", gsub("","", gsub(" ", "_", simpleScientificName)))
  ) 
# 
# colnames(interractionProb)[2] <- "Butterflies"
