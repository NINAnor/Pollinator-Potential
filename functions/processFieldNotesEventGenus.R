
#' @title \emph{processFieldNotesEvent}: Turns a presence only dataset into a presence absence dataset for surveyed data.

#' @description Some datasets have been reduced by GBIF to presence-only, despite the fact they are in fact presence-absence datastes. This function downloads these datasets directly from the source and adds the absences back in.
#'
#' @param focalEndpoint An endpoint through which the original dataset can be downloaded.
#' @param dataFolder A directory in which to save the data downloaded directly from the source.
#' @param datasetName The name of the dataset to be downloaded.
#' @param regionGeometry An sf object encompassing our region of study, as produced by defineRegion.
#' @param focalTaxon A dataframe giving the key and names of each taxonomic group we are downloading.
#' 
#' @return A new dataset with absences added back in.
#'
#' @import sf
#' 
#' 
processFieldNotesEventGenus <- function(focalEndpoint, tempFolder, datasetName, regionGeometry, taxonKey, newCrs, projcrs, focalData) {
  
  # Get the relevant endpoint
  
  # Download and unzip file in temp folder
  options(timeout=100)
  shortDataName <- str_replace_all(datasetName, " ", "")
  download.file(focalEndpoint, paste0(tempFolder,"/",  shortDataName ,".zip"), mode = "wb")
 # outputPath <- paste0(getwd(), "/", file.path(tempFolderName, shortDataName))
  outputPath <- file.path(tempdir(), shortDataName)
 # outputPath <- file.path(getwd(), paste0(tempFolderName, "/", shortDataName, "1"))
 # dir.create(outputPath)
 # dir.exists(tempFolderName)
 # dir.exists(outputPath)
 # TAF::mkdir(outputPath)
  #dir.create(outputPath, recursive = FALSE, mode = "wb")
  #dirname(outputPath)
 # outputPath <- paste0(tempFolder,"/",shortDataName)
 #dir.create(outputPath)
  #outputDir <- paste0("Data\\")
  unzip(paste0(tempFolder,"/",shortDataName ,".zip"), exdir =  outputPath, overwrite = TRUE)#paste0(tempFolderName,"/",   shortDataName))
  
  # Species records
  surveySpecies <- data.frame(surveyedSpecies = c("Diptera", "Lepidoptera", "Hymenoptera"),
                              taxonKey = c(811, 797, 1457))  
  # Function to change file names
  foo = function(x){
    paste(toupper(substring(x, 1, 1)),
          tolower(substring(x, 2, nchar(x))),
          sep = "")
  }
  
  # Load in event and occurrence data
  events <- read.delim(paste0(outputPath ,"/event.txt"))
  occurrence <- read.delim(paste0(outputPath,"/occurrence.txt"))
  
  # Get all species surveyed
  if("order" %in% colnames(occurrence)){
  surveyedSpecies <-  unique(occurrence$scientificName)
  } else {
   occurrence <- occurrence %>%
      dplyr::filter(occurrenceID %in% focalData$occurrenceID)
    occurrence$order <- focalData$order[match(occurrence$occurrenceID, focalData$occurrenceID)]
    surveyedSpecies <-  unique(occurrence$scientificName)
  }
  
  #surveyedSpecies <- surveyedSpecies[surveyedSpecies %in% c("Diptera", "Lepidoptera", "Hymenoptera")]
  
  # Find only eventIDs within our regionGeometry
  eventLocations <- events %>%
    filter(!is.na(decimalLatitude) & !is.na(decimalLongitude)) %>%
    dplyr::select(decimalLatitude, decimalLongitude, eventID, coordinateUncertaintyInMeters) %>%
    distinct()
  eventLocationsSF <- st_as_sf(eventLocations,                         
                               coords = c("decimalLongitude", "decimalLatitude"),
                               crs = projcrs)%>%
    st_transform(., newCrs)
  
  eventLocationsSF <- st_intersection(eventLocationsSF, regionGeometry) %>%
    filter(coordinateUncertaintyInMeters <= 200)
  
  # At this point we may find that there are no relevant points from this dataset available - 
  # in this case we want to finish the function early
  if (nrow(eventLocationsSF) == 0) {
    return(NULL)
  }
  
  # Make sure we have a 'year' column
  if (!("year" %in% colnames(events))) {
    if ("eventDate" %in% colnames(occurrence)) {
      eventTable <- distinct(occurrence[,c("eventID", "eventDate")])
      eventTable$year <- format(as.Date(eventTable$eventDate), "%Y") 
      events$year <- eventTable$year[match(events$eventID, eventTable$eventID)]
    } else {
      events$year <- NA
    } }
  
  # Get a dates table to match years to events
  eventDates <- events %>%
    filter(eventID %in% eventLocationsSF$eventID) %>%
    dplyr::select(year, eventID) %>%
    distinct()
  
  # Build a species table
  # speciesLegend <- data.frame(surveyedSpecies = surveyedSpecies, 
  #                             acceptedScientificName = sapply(surveyedSpecies, FUN = findGBIFName),
  #                             taxonKey = sapply(surveyedSpecies, FUN = function(x) {taxaCheck(x, taxonKey)})) %>%
  #   filter(!is.na(taxonKey))
  taxonKey <- unique(focalData$taxonKeyProject[match(occurrence$occurrenceID, focalData$occurrenceID)])
  genusLegend <- data.frame(surveyedSpecies = surveyedSpecies,
                            taxonKey = taxonKey[!is.na(taxonKey)])
  
  # Create table with all data combinations that we can match to
  allSurveyedSpecies <- focalData$acceptedScientificName #c("Diptera", "Lepidoptera", "Hymenoptera")
  eventTable <- expand.grid(order = allSurveyedSpecies, eventID = unique(eventLocationsSF$eventID))
  eventTable <- merge(eventTable, eventDates, all.x = TRUE, by = "eventID")
  
  # Create an individual count 
  occurrence$individualCount <- 1
  
  # Add in occurrence data, an NA in individualCount column means the species was NOT found in the survey
  eventTableWithOccurrences <- merge(eventTable, occurrence[,c("eventID", "order", "individualCount")], all.x = TRUE,
                                     by.x = c("order", "eventID"), by.y = c("order", "eventID"))
  eventTableWithOccurrences$individualCount[is.na(eventTableWithOccurrences$individualCount)] <- 0
  
  eventTableWithOccurrences <- eventTableWithOccurrences%>%
    group_by(order, eventID)%>%
    dplyr::summarise(individualCount = sum(individualCount))%>%
    dplyr::mutate(year = NA)%>%
    dplyr::ungroup()
  
  eventTableWithOccurrences$geometry <- eventLocationsSF$geometry[match(eventTableWithOccurrences$eventID, eventLocationsSF$eventID)]
  #eventTableWithOccurrences$acceptedScientificName <- speciesLegend$acceptedScientificName[match(eventTableWithOccurrences$species, speciesLegend$surveyedSpecies)]
  
  # Add final columns
  eventTableWithOccurrences$dataType <- "PA"
  eventTableWithOccurrences$taxonKey <-  genusLegend$taxonKey[match(eventTableWithOccurrences$order, genusLegend$surveyedSpecies)]
  eventTableWithOccurrences$taxa <- eventTableWithOccurrences$taxonKey #focalTaxon$taxa[match(eventTableWithOccurrences$taxonKey, focalTaxon$key)]
  eventTableWithOccurrences$taxonKeyProject <- eventTableWithOccurrences$taxonKey #focalTaxon$key[match(eventTableWithOccurrences$taxonKey, focalTaxon$key)]
  
  # We need to find out the total number of species in each order that we have summed over.
  totalSpeciesinOrder <- lapply(surveyedSpecies, function(x){
   ret <- focalData %>%
      dplyr::filter(order %in% x)
   
   ret <- length(unique(ret$acceptedScientificName))
   return(ret)
  })%>%
    do.call("rbind",. )%>%
    as.data.frame()%>%
    dplyr::mutate(surveyedSpecies = surveyedSpecies)
colnames(totalSpeciesinOrder) <- c("counts", "order")

eventTableWithOccurrences$counts <- totalSpeciesinOrder$counts[match(eventTableWithOccurrences$order, totalSpeciesinOrder$order)]
eventTableWithOccurrences$taxonKey <- surveySpecies$taxonKey[match(eventTableWithOccurrences$order, surveySpecies$surveyedSpecies)]
eventTableWithOccurrences$taxa <- surveySpecies$taxonKey[match(eventTableWithOccurrences$order, surveySpecies$surveyedSpecies)]
eventTableWithOccurrences$taxonKeyProject <- surveySpecies$taxonKey[match(eventTableWithOccurrences$order, surveySpecies$surveyedSpecies)]

  # eventTableWithOccurrences <- merge(eventTableWithOccurrences, as.data.frame(focalData[, c("acceptedScientificName", "genus", "order")])[, 1:3],
  #                                    all.x = TRUE, by = c("acceptedScientificName"))
  # 
  # New dataset is ready!
  newDataset <- st_as_sf(eventTableWithOccurrences,          
                         crs = newCrs)
  newDataset <- newDataset %>%
    dplyr::select(individualCount, geometry, dataType, taxa, year, taxonKeyProject, order, counts) %>%
    dplyr::group_by(individualCount, geometry, dataType, taxa, year, taxonKeyProject, order, counts) %>%
    distinct()%>%
    dplyr::ungroup() %>%
    filter(!is.na(order))%>%
    dplyr::mutate(Taxon = ifelse(order == "Lepidoptera", "Butterflies",
                                 ifelse(order == "Hymenoptera", "Bees",
                                        "Hoverflies")))
  return(newDataset)
}
