
beesSpatialData <- function(dataFolder,
                                  allPollinatorsDataset,
                                  crs = crs,
                                  print = FALSE){
  
BeesForVect <- read.csv(paste0(dataFolder, "/Artskart exports/bees/Bees from ADB removed 0 individuals.csv"), sep = ";")
BeesForVect <- BeesForVect[!BeesForVect$individualCount %in% 0,]
BeesForVect <- BeesForVect[c("validScientificName","Lon", "Lat")]
BumbleBeesForVec <- BeesForVect[grepl("Bombus",BeesForVect$validScientificName),]
SolitaryBeesForVec <- BeesForVect[!grepl("Bombus",BeesForVect$validScientificName),]
BumbleBeesForVecSensuStricto <- BumbleBeesForVec[grepl("(Bombus s. str.)",BumbleBeesForVec$validScientificName),]
BumbleBeesForVecNonSensuStricto <- BumbleBeesForVec[!grepl("(Bombus s. str.)",BumbleBeesForVec$validScientificName),]
BumbleBeesForVecNonSensuStricto$validScientificName<- do.call(c,lapply(strsplit(BumbleBeesForVecNonSensuStricto$validScientificName," "),function(x)paste(x[1],x[3])))
BumbleBeesForVecSensuStricto$validScientificName<- do.call(c,lapply(strsplit(BumbleBeesForVecSensuStricto$validScientificName," "),function(x)paste(x[1],x[5])))
BindedBeesForVec <- rbind(BumbleBeesForVecNonSensuStricto,
                          BumbleBeesForVecSensuStricto,
                          SolitaryBeesForVec )
FilteredBindedBeesForVec <- BindedBeesForVec[BindedBeesForVec$validScientificName %in% allPollinatorsDataset$ValidSpeciesName,]
BeePoints <- vect(FilteredBindedBeesForVec, 
                  geom=c("Lon", "Lat"),
                  crs=projCRS, keepgeom=FALSE)

if(print) plot(BeePoints)

BeePoints$Taxon <- "Bees"
return(BeePoints)
}
