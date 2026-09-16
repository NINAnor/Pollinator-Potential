
butterflySpatialData <- function(dataFolder,
                                   allPollinatorsDataset,
                                   crs = crs,
                                   print = FALSE){
butterfliesForVect <- read.csv("Data/Artskart exports/Selected lepidoptera families/Selected lepidoptera families.csv", sep = ";")
butterfliesForVect <- butterfliesForVect[butterfliesForVect$individualCount > 0,]
butterflySelectedSpeciesForVec <- butterfliesForVect[butterfliesForVect$validScientificName %in% allPollinatorsDataset$ValidSpeciesName,]
butterflyPoints <- vect(butterflySelectedSpeciesForVec, geom=c("longitude.1", "latitude.1"), crs= crs, keepgeom=FALSE)
butterflyPoints <- butterflyPoints[,"validScientificName"]
if(print) plot(butterflyPoints)
butterflyPoints$Taxon <- "Butterflies"
return(butterflyPoints)
}