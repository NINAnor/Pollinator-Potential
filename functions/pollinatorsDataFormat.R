#########################################################################
## Load list of Norwegian bees and butterflies
#########################################################################


formatPollinatorDataset <- function(dataFolder){
  ValidNorwegianPollinators <- read.csv(paste0(dataFolder,"/Artsnavnebase/Lepidoptera and bee names.csv"), sep = ";")
  ValidNorwegianPollinators$ValidSpeciesName <- with(ValidNorwegianPollinators, paste(Slekt,Art))
  ValidNorwegianPollinators[grepl("Bombus",ValidNorwegianPollinators$ValidSpeciesName),]
  
  #select the following pollinators with genera
  ValidNorwegianPollinators <- ValidNorwegianPollinators[!ValidNorwegianPollinators$Slekt %in% c("Nomada",
                                                                                                 "Stelis",
                                                                                                 "Sphecodes",
                                                                                                 "Coelioxys",
                                                                                                 "Biastes",
                                                                                                 "Epeolus"),]
  
  return(ValidNorwegianPollinators)
}