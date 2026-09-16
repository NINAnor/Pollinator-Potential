
#' @title \emph{processNationalInsectMonitoring}: Process data given by the National insect Monitoring survey from NINA.

#' @description The National Insect Monitoring Survey by NINA has a very specific data format (read more in its page on the GBIF website). As such, it needs to be processed ina  particular way.
# '
#' @param focalData The Insect monitoring dataset as imported by the GBIF download.
#' @param endpoint The endpoint for this dataset.
#' @param tempFolderName The temporary folder in which we save the data.
#' 
#' @return A new, processed dataset
#'
#' @import sf
#' @import stringr



processNationalInsectMonitoringWithSpecies <- function(focalData, endpoint, tempFolderName) {
  
  # Download and unzip file in temp folder
  
  # zippedDownload <- paste0(tempFolderName,"/NationaInsectMonitoring.zip")
  # download.file(endpoint, zippedDownload, mode = "wb")
  # unzip(paste0(tempFolderName,"/NationaInsectMonitoring.zip"), exdir = paste0(tempFolderName,"/NationaInsectMonitoring"))
  # 
  # # Read in occurrence and event data
  # event_raw <- read.delim(paste0(tempFolderName, "/NationaInsectMonitoring/event.txt"))
  #events <- events[events$coordinateUncertaintyInMeters <= 100,]
  
  if(!file.exists("GBIF_data")){
    options(timeout=100)
    
    dataset_id <- "19fe96b0-0cf3-4a2e-90a5-7c1c19ac94ee" ##From the webpage URL
    # Suggested citation: Take the citation as from downloaded from GBIF website, replace "via GBIF.org" by endpoint url. 
    tmp <- tempfile()
    download.file(paste0("http://api.gbif.org/v1/dataset/",dataset_id,"/document"),tmp) # get medatadata from gbif api
    meta <- read_xml(tmp) %>% as_list() # create list from xml schema
    gbif_citation <- meta$eml$additionalMetadata$metadata$gbif$citation[[1]] # extract citation
    
    
    # Fetching the dataset URL
    dataset_url <-  paste0("http://api.gbif.org/v1/dataset/",dataset_id,"/endpoint")
    dataset <- RJSONIO::fromJSON(dataset_url)
    endpoint_url <- dataset[[1]]$url # extracting URL from API call result
    
    # create directory and download the GBIF Data
    system("mkdir -p GBIF_data")
    download.file(endpoint_url, destfile = "GBIF_data/gbif_download.zip", mode="wb")
    
    
    unzip("GBIF_data/gbif_download.zip",
          exdir = "GBIF_data")
  }
  
  event_raw <- read_delim("GBIF_data/event.txt", 
                          delim = "\t",
                          locale = locale(encoding = "UTF-8"),
                          progress = FALSE,
                          show_col_types = FALSE)
  
  occurrence_raw <- read_delim("GBIF_data/occurrence.txt", 
                               delim = "\t",
                               locale = locale(encoding = "UTF-8"),
                               progress = FALSE,
                               show_col_types = FALSE
  )
  
  #pollinatorFolder <- "Data/pollinatorDataFolder"
  #dataList <- readRDS(paste0(pollinatorFolder, "/speciesDataImported.RDS"))
  
  
  
  
  identifications_raw <- event_raw %>% 
    filter(samplingProtocol == "Level_1")
  
  sampling_trap_raw <- event_raw %>% 
    filter(samplingProtocol == "Level_2")
  
  locality_sampling_raw <- event_raw %>% 
    filter(samplingProtocol == "Level_3")
  
  year_locality_raw <- event_raw %>% 
    filter(samplingProtocol == "Level_4")
  
  identifications <- identifications_raw %>% 
    #collect() %>% 
    select(-dynamicProperties)
  
  tempDynamic <- identifications_raw %>% 
    #head() %>% 
    select(dyn = dynamicProperties) %>% 
    #collect() %>% 
    unlist() %>% 
    spread_all()  %>% 
    as_tibble() %>% 
    select(-document.id)
  
  identifications <- identifications %>% 
    left_join(tempDynamic,
              by = c("id" = "id"))
  
  
  identifications <- identifications %>% 
    select(identification_id = eventID,
           sampling_trap_id = parentEventID,
           identification_name,
           identification_comment,
           read_abundance,
           ekstraksjonsdato,
           ekstraksjonskit,
           ekstraksjonskommentar)
  
  # Sampling trap
  sampling_trap <- sampling_trap_raw %>% 
    mutate(trap_name = stringr::str_split_i(locationRemarks, ": ", 2)) %>% 
    select(-dynamicProperties) 
  
  tempDynamic <- sampling_trap_raw %>% 
    select(dyn = dynamicProperties) %>% 
    unlist() %>% 
    spread_all()  %>% 
    as_tibble() %>% 
    select(-document.id)
  
  sampling_trap <- sampling_trap %>% 
    left_join(tempDynamic,
              by = c("id" = "id"))
  
  sampling_trap <- sampling_trap %>% 
    select(sampling_trap_id = eventID,
           locality_sampling_id = parentEventID,
           trap_name,
           sample_size_value = sampleSizeValue,
           sample_size_unit = sampleSizeUnit,
           trap_type,
           trap_model,
           liquid_name,
           ethanol_prc_at_lab,
           net_wet_weight,
           point_latitude = decimalLatitude,
           point_longitude = decimalLongitude,
           point_wkt = footprintWKT,
           coord_unc_in_meters = coordinateUncertaintyInMeters,
           coordinate_system = geodeticDatum) 
  
  
  # Location samplings
  
  locality_sampling <- locality_sampling_raw %>% 
    mutate(sampling_event = stringr::str_split_i(locationRemarks, ": ", 2)) %>% 
    select(-dynamicProperties) 
  
  tempDynamic <- locality_sampling_raw %>% 
    select(dyn = dynamicProperties) %>% 
    unlist() %>% 
    spread_all()  %>% 
    as_tibble() %>% 
    select(-document.id)
  
  locality_sampling <- locality_sampling %>% 
    left_join(tempDynamic,
              by = c("id" = "id"))
  
  locality_sampling <- locality_sampling %>% 
    separate(eventTime, 
             c("start_time", "end_time"), 
             sep = "/") %>% 
    mutate(start_time = as.POSIXct(start_time),
           end_time = as.POSIXct(end_time))
  
  locality_sampling <- locality_sampling %>% 
    select(locality_sampling_id = eventID,
           year_locality_id = parentEventID,
           sampling_event,
           sample_size_value = sampleSizeValue,
           sample_size_unit = sampleSizeUnit,
           start_time,
           end_time,
           sampling_min_temp,
           sampling_max_temp,
           sampling_avg_temp)
  
  
  # Year samplings
  year_locality <- year_locality_raw %>% 
    select(locality_vernacular_name = locality,
           everything()) %>% 
    select(-dynamicProperties)
  
  tempDynamic <- year_locality_raw %>% 
    select(dyn = dynamicProperties) %>% 
    unlist() %>% 
    spread_all()  %>% 
    as_tibble() %>% 
    select(-document.id)
  
  
  year_locality <- year_locality %>% 
    left_join(tempDynamic,
              by = c("id" = "id"))
  
  year_locality <- year_locality %>% 
    separate(eventTime, 
             c("start_time", "end_time"), 
             sep = "/") %>% 
    mutate(season_start_time = as.POSIXct(start_time),
           season_end_time = as.POSIXct(end_time)) %>%  
    mutate(ssbid = as.character(ssbid)) 
  
  year_locality <- year_locality %>% 
    select(year_locality_id = eventID,
           season_sample_size_value = sampleSizeValue,
           season_sample_size_unit = sampleSizeUnit,
           season_start_time,
           season_end_time,
           locality,
           year,
           habitat_type,
           region_name,
           ssbid,
           ano_flate_id,
           no_herb_spec,
           avg_prc_cov_herb_species,
           no_tree_spec,
           dom_tree_spec,
           avg_dom_tree_age,
           centroid_latitude = decimalLatitude,
           centroid_longitude = decimalLongitude,
           polygon_wkt = footprintWKT,
           coordinate_system = geodeticDatum
    ) 
  
  occurrence <- occurrence_raw %>% 
    select(occurrence_id = occurrenceID,
           identification_id = eventID,
           quantity = organismQuantity,
           quantity_type = organismQuantityType,
           scientificName = scientificName,
           vernacular_name = vernacularName,
           scientific_name_id = taxonID,
           class,
           order,
           family,
           genus,
           specific_epithet = specificEpithet)%>%
    dplyr::filter(order %in% c("Diptera", "Lepidoptera", "Hymenoptera"))
  
  
  events <- identifications %>%
    left_join(sampling_trap, by = c("sampling_trap_id" = "sampling_trap_id")) %>% 
    left_join(locality_sampling,
              by = c("locality_sampling_id" = "locality_sampling_id")) %>%
    left_join(year_locality,
              by = c("year_locality_id" = "year_locality_id")) %>% 
    filter(trap_type == "Malaise") %>% 
    mutate(julian_day = lubridate::yday(end_time),
           year = forcats::as_factor(year))%>%
    filter(habitat_type %in% "Semi-nat")
  
  
  
  #sum(occurrence$identification_id %in% identifications$identification_id)
  
  # There are four levels to this thing. Level 1 events are linked to level 2 by their parent ID and so forth.
  # Level 1 - Identification of an insect in a trap
  # Level 2 - Laying of the trap itself
  # Level 3 - Multiple traps within a sampling period
  # Level 4 - Multiple sampling periods within a sampling season
  
  # First, get a list of all species samples and create a legend for them
  surveyedSpecies <- unique(occurrence$scientificName)
  legendTable <- data.frame(scientificName = surveyedSpecies, genus = word(surveyedSpecies,1))
  
  # Narrow it down to the same genus we're looking for
  legendTable <- legendTable[legendTable$genus %in% word(unique(focalData$acceptedScientificName), 1),]
  legendTable$acceptedScientificName <- sapply(legendTable$scientificName, FUN = findGBIFName)
  
  # And now narrow it down to species
  legendTable <- legendTable[legendTable$acceptedScientificName %in% focalData$acceptedScientificName,]
  
  # Here I've taken level 1 events and linked them to level 2. Each location has a separate date for the trap sampling.
  occurrencesWithEvent <- occurrence[occurrence$identification_id %in% events$identification_id,]
  
  # Now get all key data columns from the events data
  occurrencesWithEvent$parentEventID <- events$sampling_trap_id[match(occurrencesWithEvent$identification_id , events$identification_id)]
  occurrencesWithEvent$locationID <- events$locality_sampling_id[match(occurrencesWithEvent$identification_id , events$identification_id)]
  occurrencesWithEvent$eventDate <- events$season_end_time[match(occurrencesWithEvent$identification_id , events$identification_id)]
  occurrencesWithEvent$decimalLongitude <- events$point_longitude[match(occurrencesWithEvent$identification_id , events$identification_id)]
  occurrencesWithEvent$decimalLatitude <- events$point_latitude[match(occurrencesWithEvent$identification_id , events$identification_id)]
  
  # Take only most recent event from a location
  occurrencesWithEvent$exactDate <- as.Date(substr(occurrencesWithEvent$eventDate,1,10))
  mostRecentSample <- occurrencesWithEvent %>%
    dplyr::select(parentEventID, locationID, exactDate) %>%
    group_by(locationID) %>%
    slice_max(exactDate, n = 1,na_rm = TRUE) %>% 
    distinct() %>%
    ungroup()
  mostRecentSample$year <- substr(mostRecentSample$exactDate, 1, 4)
  
  # Now add a line for every species
  ourDataset <- merge(mostRecentSample, legendTable["scientificName"], all = TRUE)
  ourDataset$acceptedScientificName <- legendTable$acceptedScientificName[match(ourDataset$scientificName, legendTable$scientificName)]
  ourDataset <- ourDataset[!is.na(ourDataset$acceptedScientificName),]
  
  # Now get abundances for those specie that have them
  ourDatasetAbundance <- merge(ourDataset, occurrencesWithEvent[,c("scientificName", "quantity", "parentEventID",
                                                                   "locationID", "genus")], all.x = TRUE, 
                               by = c("parentEventID", "locationID", "scientificName"))
  
  # Now create a locations table to get location data
  locations <- distinct(occurrencesWithEvent[,c("locationID", "decimalLatitude", "decimalLongitude")])
  ourDatasetLocated <- merge(ourDatasetAbundance, locations, all.x = TRUE, by = "locationID")
  
  # Convert all NAs to 0
  ourDatasetLocated$quantity[is.na(ourDatasetLocated$quantity)] <- 0
  
  # Convert this into an sf object
  newDataset <- st_as_sf(ourDatasetLocated, coords = c("decimalLongitude", "decimalLatitude"),
                         crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  
  # Cut down set
  newDataset <- newDataset[c("acceptedScientificName", "quantity", "year", "genus")] %>%
    rename(individualCount = quantity)%>%
    st_transform("+proj=longlat +ellps=WGS84")
  
  # Crop to relevant region
  #st_crs(newDataset) <- "+proj=longlat +ellps=WGS84"
  #newDataset <- st_intersection(newDataset, regionGeometry)
  newDataset <- st_transform(newDataset, crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 ")
  newDataset$dataType <- "Counts"
  
  taxaLegend <- distinct(st_drop_geometry(focalData[,c("taxa", "acceptedScientificName", "taxonKeyProject")]))
  
  newDataset$taxa <- taxaLegend$taxa[match(newDataset$acceptedScientificName, taxaLegend$acceptedScientificName)]
  newDataset$taxonKeyProject <- taxaLegend$taxonKeyProject[match(newDataset$acceptedScientificName, taxaLegend$acceptedScientificName)]
  return(newDataset)
  
  
}