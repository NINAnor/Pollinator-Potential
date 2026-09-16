
hoverfliesSpatialData <- function(dataFolder,
         allPollinatorsDataset,
         crs = crs,
         print = FALSE){
  HoverFlyOccurrences <- read.csv(paste0(dataFolder, "/Artskart exports/Syrphidae/Syrphidae added LatLon.csv"), sep = ";")
  HoverflySelectedSpeciesForVec<- HoverFlyOccurrences[HoverFlyOccurrences$Species %in% allPollinatorsDataset$ValidSpeciesName,]
  HoverflyPoints <- vect(HoverflySelectedSpeciesForVec, 
                         geom=c("Lon", "Lat"), 
                         crs=crs, 
                         keepgeom=FALSE)
  names(HoverflyPoints) <- "validScientificName"
  if(print) plot(HoverflyPoints)
  HoverflyPoints$Taxon <- "Hoverflies"
  return(HoverflyPoints)
}


