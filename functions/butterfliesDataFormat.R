#########################################################################
## Load Lepidoptera hostplants
#########################################################################

formatButterflyDataset <- function(dataFolder,ValidNorwegianPollinators, print = FALSE){
  
  #Load EU butterfly Indicators
  EUbutterflyIndicators <- c("Lasiommata megera",
                             "Coenonympha pamphilus",
                             "Lycaena phlaeas",
                             "Ochlodes sylvanus",
                             "Polyommatus icarus",
                             "Thymelicus acteon",
                             "Anthocharis cardamines",
                             "Cupido minimus",
                             "Cyaniris semiargus",
                             "Erynnis tages",
                             "Lysandra bellargus",
                             "Lysandra coridon",
                             "Maniola jurtina",
                             "Euphydryas aurinia",
                             "Phengaris arion",
                             "Phengaris nausithous",
                             "Spialia sertorius")
  
  EuropeLepidopteraHosts <- read.csv(paste0(dataFolder, "/Lepidoptera hosts/resource.csv"))
  EuropeLepidopteraHosts$ValidSpeciesName <- with(EuropeLepidopteraHosts, paste(Insect.Genus,Insect.Species))
  
  if(print) head(EuropeLepidopteraHosts)
  
  EuropeanButterfliesWithHostsInNorway <- EuropeLepidopteraHosts[EuropeLepidopteraHosts$ValidSpeciesName %in% ValidNorwegianPollinators$ValidSpeciesName,]
  
  NorwegianButterfliesWithHosts <- merge(ValidNorwegianPollinators,EuropeLepidopteraHosts,by = "ValidSpeciesName")
  # Only select the families:
  
  EngSommerFuglerNorge <- NorwegianButterfliesWithHosts[NorwegianButterfliesWithHosts$Familie %in% c("Geometridae","Zygaenidae","Hesperiidae","Papilionidae", "Pieridae", "Lycaenidae","Nymphalidae") ,]
  
  # check for typos and remove species with no information on which genera they depend on
  if(print) sort(unique(EngSommerFuglerNorge$Hostplant.Species))
  
  EngSommerFuglerNorge[EngSommerFuglerNorge$Hostplant.Species %in% c("arvense","arvensis"),] #ok
  EngSommerFuglerNorge[EngSommerFuglerNorge$Hostplant.Species %in% c("dioica", "dioicus"),] #ok
  EngSommerFuglerNorge[EngSommerFuglerNorge$Hostplant.Species %in% c("lupulina" , "lupulus"),] #ok
  EngSommerFuglerNorge[EngSommerFuglerNorge$Hostplant.Species %in% c("silvaticum", "sylvaticum"),] # typo
  EngSommerFuglerNorge[EngSommerFuglerNorge$Hostplant.Species %in% c("vulgare","vulgaris"),]#ok
  EngSommerFuglerNorge$Hostplant.Species[EngSommerFuglerNorge$Hostplant.Species %in% "silvaticum"]  <- "sylvaticum"
  EngSommerFuglerNorge$HostPlantSpecies <- with(EngSommerFuglerNorge,paste(Hostplant.Genus,Hostplant.Species))
  EngSommerFuglerNorge <- EngSommerFuglerNorge[!EngSommerFuglerNorge$Hostplant.Genus %in% "",]
  
  EngSommerFuglerNorge[EngSommerFuglerNorge$Hostplant.Species %in% "",]$Hostplant.Genus
  
  #########################################################################################################
  ## combine bees and butterflies
  #########################################################################################################
  NorLepidopteraToCombine <- EngSommerFuglerNorge[c("ValidSpeciesName","Hostplant.Genus")]
  
  NorLepidopteraToCombine$Taxon <- "Butterfly"
  return(NorLepidopteraToCombine)
}
