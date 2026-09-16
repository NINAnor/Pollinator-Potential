#########################################################################################################
## Combine Hoverfly:PlantGenus with Bee and butterfly
#########################################################################################################

formatHoverflyDataset <- function(dataFolder, ValidNorwegianPollinators){
  HoverFlyPlantFull <- read.csv(paste0(dataFolder,
                                       "/Hoverflies results/results.csv"))
  
  HoverflyPlantSubSet <- HoverFlyPlantFull[c("Pollinator.Species","Plant.Species")]
  HoverflyPlantSubSet$PlantGenus <- do.call(c,
                                            lapply(strsplit(HoverflyPlantSubSet$Plant.Species," "),
                                                   function(x)x[1]))
  names(HoverflyPlantSubSet) <- c("Hoverfly",
                                  "PlantSpecies",
                                  "PlantGenus")
  
  # Exclude hoverflies not found in Norway during the current climate period (1991-2020)
  HoverFlyOccurrences <- read.csv(paste0(dataFolder, "/Artskart exports/Syrphidae/Syrphidae added LatLon.csv"), sep = ";")
  HoverFlyOccurrences <- HoverFlyOccurrences[grepl(" ",HoverFlyOccurrences$Species),]
  HoverFlyOccurrences$Species[HoverFlyOccurrences$Species %in% "Eriozona (Megasyrphus) erratica"] <- "Eriozona erratica"
  HoverflyPlantSubSetOccursInNorway <- HoverflyPlantSubSet[HoverflyPlantSubSet$Hoverfly %in% HoverFlyOccurrences$Species,]
  HoverflyPlantSubSetOccursInNorway$Taxon <- "Hoverflies"
  
  return(HoverflyPlantSubSetOccursInNorway)
}
