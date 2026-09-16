
# Load the required packages
source("functions/loadPackages.R")


#########################################################################
## Load list of Norwegian bees and butterflies
#########################################################################
formatBeeDataset <- function(dataFolder, validNorwegianPollinators, print = TRUE){
ValidNorwegianPollinators <- read.csv("Data/Artsnavnebase/Lepidoptera and bee names.csv", sep = ";")
ValidNorwegianPollinators$ValidSpeciesName <- with(ValidNorwegianPollinators, paste(Slekt,Art))
ValidNorwegianPollinators[grepl("Bombus",ValidNorwegianPollinators$ValidSpeciesName),]

#select the following pollinators with genera
ValidNorwegianPollinators <- ValidNorwegianPollinators[!ValidNorwegianPollinators$Slekt %in% c("Nomada",
                                                                                               "Stelis",
                                                                                               "Sphecodes",
                                                                                               "Coelioxys",
                                                                                               "Biastes",
                                                                                               "Epeolus"),]

#########################################################################
## Load Bee hostplants
#########################################################################
BeesAndPlantGenera <- read.csv("Data/Bees hosts/Bee species and plant genera.csv", sep = ";")
BeesAndPlantGenera$ValidBeeSpecies <- with(BeesAndPlantGenera, 
                                           paste(BeeGenus,BeeSpecies))
BeesAndPlantGenera$Taxon <- "Bees"


# Now, we merge the Norwegian pollinators with the plants
MergedNorwegianBeesAndHosts <- merge(ValidNorwegianPollinators,
                                     BeesAndPlantGenera,
                                     by.x = "ValidSpeciesName",
                                     by.y = "ValidBeeSpecies")

MergedNorwegianBeesAndHosts <- MergedNorwegianBeesAndHosts[c("ValidSpeciesName",
                                                             "PlantGenus",
                                                             "Taxon")]
MergedNorwegianBeesAndHosts <- MergedNorwegianBeesAndHosts[!MergedNorwegianBeesAndHosts$ValidSpeciesName %in% c("Apis mellifera",
                                                                                                                "Bombus bohemicus",
                                                                                                                "Bombus campestris",
                                                                                                                "Bombus cryptarum",
                                                                                                                "Bombus lucorum",
                                                                                                                "Bombus magnus",
                                                                                                                "Bombus terrestris",
                                                                                                                "Bombus norvegicus",
                                                                                                                "Bombus rupestris",
                                                                                                                "Bombus sylvestris"),]
if(print){ 
  message("These are the list of Norwegian Bees and Host")
  MergedNorwegianBeesAndHosts
}
return(MergedNorwegianBeesAndHosts)
}








