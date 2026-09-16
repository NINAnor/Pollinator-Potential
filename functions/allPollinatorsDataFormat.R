#Load the dataset

#Load the valid norwegian pollinators
allPollinatorsDataset <- function(dataFolder, print = FALSE){
  
ValidNorwegianPollinators <- formatPollinatorDataset(dataFolder)

# Format the bee datasets and merge with the valid pollinators
MergedNorwegianBeesAndHosts <- formatBeeDataset(dataFolder, ValidNorwegianPollinators,
                                                print = FALSE)

# Format the butterfly dataset 
NorLepidopteraToCombine <- formatButterflyDataset(dataFolder, ValidNorwegianPollinators, print = FALSE)

# Set the names of the butterfly dataset to be the same as that of the bee
# dataset and combine the two datasets
names(NorLepidopteraToCombine) <- names(MergedNorwegianBeesAndHosts)
CombinedLidopteraAndBees <- rbind(NorLepidopteraToCombine,MergedNorwegianBeesAndHosts)
#head(CombinedLidopteraAndBees)
#sort(unique(CombinedLidopteraAndBees$PlantGenus))

# format the hoverfly dataset and combine it to the butterfly and bees dataset
HoverflyPlantSubSetOccursInNorway <- formatHoverflyDataset(dataFolder,ValidNorwegianPollinators )
names(HoverflyPlantSubSetOccursInNorway)[c(1,3, 4)] <- names(CombinedLidopteraAndBees)
#names(HoverflyPlantSubSetOccursInNorway)
CombinedLidopteraAndBeesAndHoverflies <- rbind(CombinedLidopteraAndBees,HoverflyPlantSubSetOccursInNorway[names(CombinedLidopteraAndBees)])


####################################################################################
## Exclude plant genera not found in Norway
####################################################################################
PlantNamesFromADB <- read.csv(paste0(dataFolder, "/Artsnavnebase/Plantae for R.csv"), sep = ";")
CombinedLidopteraAndBeesAndHoverfliesFilteredPlantGenera <- CombinedLidopteraAndBeesAndHoverflies[CombinedLidopteraAndBeesAndHoverflies$PlantGenus %in% PlantNamesFromADB$Slekt,]
CombinedLidopteraAndBeesAndHoverfliesFilteredPlantGenera


##########################################################################################################
## Exclude woody plant genera 
##########################################################################################################

LedaWood <- read.csv(paste0(dataFolder, "/LedaWoodiness/leda wood.csv"), sep = ";",header = TRUE)


fn.getPlantGenusName <- function(x){
  if(length(x) > 2){
    Genus = x[2] 
  } else {
    Genus = x[1]
  }
  return(Genus)}

LedaWood$PlantGenus <- do.call(c,lapply(strsplit(LedaWood$SBS.name, " "),fn.getPlantGenusName))

# First check how many plant genera are found in the Leda data base #

table(unique(CombinedLidopteraAndBeesAndHoverfliesFilteredPlantGenera$PlantGenus) %in% unique(LedaWood$PlantGenus))
# 341 plant genera occur in the Leda data base, 23 do not.
data.frame(PlantGenus = unique(CombinedLidopteraAndBeesAndHoverfliesFilteredPlantGenera$PlantGenus[!CombinedLidopteraAndBeesAndHoverfliesFilteredPlantGenera$PlantGenus %in% unique(LedaWood$PlantGenus)]))


# Any plant genus whoose members are typically woody or semi-woody are removed (i.e. their sum > sum of non-woody)

SelectedLedaGenera <- LedaWood[LedaWood$PlantGenus %in% CombinedLidopteraAndBeesAndHoverfliesFilteredPlantGenera$PlantGenus,]
#head(SelectedLedaGenera)
FreqWoodiness <- as.data.frame.matrix(with(SelectedLedaGenera, table(PlantGenus,woodiness)))
FreqWoodiness$TypicallyWoody <- rowSums(FreqWoodiness[,2:3]) > FreqWoodiness[,1]
#head(FreqWoodiness)


# Now load plant genera found in Norway but not in Leda


## Now exclude all plant genera that are not typically woody
## Manually classified in excel so that predominantly herbaceous plants will be retained
NonLedaWood <- read.csv(paste0(dataFolder,"/LedaWoodiness/Plant genera not in Leda.csv"), sep = ";",header = TRUE)
NonLedaGeneraToInclude <- NonLedaWood[NonLedaWood$TypicallyWoody %in% TRUE,"Plant"]

## Select lepidoptera, bees, and plant-pollintor links
FinalSelectionOfLepidopteraAndBeesAndHoverflies <- CombinedLidopteraAndBeesAndHoverfliesFilteredPlantGenera[CombinedLidopteraAndBeesAndHoverfliesFilteredPlantGenera$PlantGenus %in%sort(c(rownames(FreqWoodiness[FreqWoodiness$TypicallyWoody %in% FALSE,]),NonLedaGeneraToInclude)),]
#head(FinalSelectionOfLepidopteraAndBeesAndHoverflies)
return(FinalSelectionOfLepidopteraAndBeesAndHoverflies)
}
