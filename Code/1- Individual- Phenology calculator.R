## Luke Ozsanlav-Harris
##Created: 13/05/2020


## Script Aims:
## Calculating spring departure and arrival dates for the Iceland to Greenland migratory leg
## Will use the absolute displacement method  from Soriano-Redondo et al (2020) Ibis


##packages required
pacman::p_load(tidyverse, data.table, leaflet, 
               geosphere, lubridate, RColorBrewer)


##** FIXES **##
## 1. All of the checks and visualizations need updating in each section



##
#### 1. Importing Ecotone and Ornitela GPS data ####
##


## Read in ecotone data
Eco_data <- readRDS(file = "Tracking data/Eco_data_cleaned.RDS")

## Extract GPS data and a subset of the columns
Eco_gps <- Eco_data %>% filter(is.na(Acc_X) == TRUE) # extract just GPS data
setnames(Eco_gps, old = "Tag_ID", new = "device_id") # change name of this column
Eco_gps_short <- Eco_gps %>% select(Latitude, Longitude, UTC_datetime, device_id)
Eco_gps_short$UTC_datetime <- as.character(Eco_gps_short$UTC_datetime)
rm(Eco_data, Eco_gps) # remove large raw data set


## Read in Ornitela gps data
orn_gps <- readRDS("Tracking data/Ornitela_GPSClean.RDS")
orn_gps_short <- orn_gps %>% select(Latitude, Longitude, UTC_datetime, device_id)
orn_gps_short$UTC_datetime <- as.character(orn_gps_short$UTC_datetime)
rm(orn_gps)





##
#### 2. Combining GPS data sets ####
##

##bind all the GPS data together
All_GPS <- rbind(Eco_gps_short, orn_gps_short)
rm(Eco_gps_short, orn_gps_short)

## parse the datetime stamp
All_GPS$UTC_datetime <- as.POSIXct(All_GPS$UTC_datetime, format= "%Y-%m-%d %H:%M:%S")

##order by tag and then date
All_GPS <- All_GPS[order(All_GPS$device_id, All_GPS$UTC_datetime),]

## filter out this tag ID because it should not exist
## Need to go back to the raw data and fix it first
All_GPS <- filter(All_GPS, !ID == "17809")

##change column names so they match the code used from paper
All_GPS %>% setnames(old=c("Longitude", "Latitude",  "UTC_datetime", "device_id"), 
                     new=c("lon", "lat", "timestamp", "ID"))

##select a subset of tags to work with so code run quicker while trying, delete later
# All_GPS2 <- All_GPS %>%  filter(ID %in% c("WHIT17", "WHIT02", "UCOL40", 
#                                           "HAR27", "BLO09", "17762", "17813"))




##
#### 3. Absolute displacement (AD) fucntion for longitudinal migration ####
##

##Turning the loop in the paper into a function that you can reuse
##now have a function that does AD calculation for migration with large longitudinal component
AD_method_lon <- function(DayADf, IDs, col_heads, MinDDisp, Max_lon, Min_lon) {
  # loop for each bird 
  for(z in 1:length(IDs)) { 
    DayAD <- filter(DayADf, ID == IDs[z])
    
    # onset and end of migration 
    for(i in 2:(nrow(DayAD)-1)) {
      if(DayAD$Migration[i] != "Migration")
      {next} # find start of unmistakable migratory segment 
      if(DayAD$Migration[i-1] == "Migration" & DayAD$Migration[i+1] == "Migration") 
      {next} # skip intermidate locations of unmistakable migratory segment 
      # find onset of autumn migration or onset of spring migration 
      if(DayAD$Migration[i-1] != "Migration") {
        j=i-1 
        for (j in j:1){ 
          if(DayAD$ID[j] !=
             DayAD$ID[j+1]) {break} 
          if((DayAD$Start_lon[j] <= -Min_lon & DayAD$Start_lon[j] >= -Max_lon)
             | DayAD$DDisp[j] > MinDDisp | DayAD$BDISP[j] > MinDDisp)
          {DayAD$Migration[j] <- "Migration" 
          next}
          if((DayAD$Start_lon[j] > -Min_lon & DayAD$Start_lon[j] < -Max_lon) 
             | DayAD$DDisp[j] <= MinDDisp) {break}}}
      # find end of autumn migration or end of spring migration 
      if(DayAD$Migration[i+1] != "Migration") { 
        j=i+1
        for (j in j:nrow(DayAD)) {
          if(DayAD$ID[j] !=
             DayAD$ID[j-1]) {break}
          if((DayAD$Start_lon[j] <= -Min_lon & DayAD$Start_lon[j] >=-Max_lon) 
             | DayAD$DDisp[j] > MinDDisp | DayAD$BDISP[j] > MinDDisp) 
          {DayAD$Migration[j] <- "Migration" 
          next} 
          if((DayAD$Start_lon[j] > -Min_lon & DayAD$Start_lon[j] < -Max_lon) 
             | DayAD$DDisp[j] <= MinDDisp) {break}}}} 
    col_heads <- rbind(col_heads, DayAD)}
  return(col_heads)
  }
  




##
#### 4. Calculate phenology for Ice-->Green Spring migration ####
##


##
#### 4.1 Calculate daily travel distances ####
##

## change data frame for AD method
AD_Ice <- All_GPS 

## extract the date from the timestamp
AD_Ice$date <- as.Date(AD_Ice$timestamp) 
AD_Ice$yday <- yday(AD_Ice$timestamp) 

## label if a GPS fix is in the unmistakable migratory segment (Over the Atlantic between Iceland and Greenland)
AD_Ice$Migration <- ifelse(AD_Ice$lon <= -25 & AD_Ice$lon >= -45 & AD_Ice$lat >= -62.5 &  AD_Ice$yday < 185, 1, 0)

## Create data frame with first and last location for each tag day
DayAD1_Ice<- AD_Ice %>% 
              group_by(ID, date) %>% 
              summarise(Start_lat = first(lat), Start_lon = first(lon), 
                        End_lat = last(lat), End_lon = last(lon),
                        n=n(), Migration=sum(Migration)) 

## Change migration column to either 1 or zero depending if there were any fixes in the unmistakable migratory segment
DayAD1_Ice$Migration <- ifelse(DayAD1_Ice$Migration==0, 0, 1)

## Now calculate how far bird has gone been the start and end of each day (DDISP)
DayAD1_Ice$DDisp<- distGeo(cbind(DayAD1_Ice$Start_lon, DayAD1_Ice$Start_lat), 
                       cbind(DayAD1_Ice$End_lon, DayAD1_Ice$End_lat))/1000


## Now calculate how far bird has gone been start of one day and end of previous (BDISP)
## create a tag_year column
DayAD1_Ice <- DayAD1_Ice %>% mutate(year=year(ymd(date)), tag_year = paste(ID, year, sep="_"))

##function to calculate distance between start of one day and end of last
Diff_day_dist_calc <- function(x){
  
  distGeo(cbind(x$Start_lon, x$Start_lat), 
          cbind(lag(x$End_lon), lag(x$End_lat)))/1000
  
}

##now create column with this BDISP
DayAD1_Ice <- DayAD1_Ice %>% 
              group_by(tag_year) %>% 
              nest() %>% 
              mutate(BDISP = map(data, Diff_day_dist_calc)) %>% 
              unnest(cols = c(data, BDISP))

DayAD1_Ice$BDISP <- ifelse(is.na(DayAD1_Ice$BDISP) == TRUE, 1, DayAD1_Ice$BDISP)

##reset the BDISP colum to 1 if the number of days between rows is greater than 1
DayAD1_Ice$DayDiff <- ymd(DayAD1_Ice$date) - ymd(lag(DayAD1_Ice$date))
DayAD1_Ice$BDISP <- ifelse(DayAD1_Ice$DayDiff > 2, 1, DayAD1_Ice$BDISP)




##
#### 4.2 Classify migrations ####
##


## First let's remove any  birds that didn't make it to Greenland in a particular year
MaxDists <- DayAD1_Ice %>% 
            group_by(tag_year) %>% 
            summarise(Maxlon = min(End_lon))
## filter out the ones that did make it, -45 is half way across the Greenalnd Icecap
MaxDists <- MaxDists[MaxDists$Maxlon < -45,]
## now filter out the birds that made it to Greenland 
survivors <- unique(MaxDists$tag_year)
DayAD1_Ice <- filter(DayAD1_Ice, tag_year %in% survivors)


## Unmistakable migratory segment (e.g. over North Atlantic and Greenland Ice cap) was determined above
## Add another condition in case birds crossed this segment in one day
table(DayAD1_Ice$Migration)
DayAD1_Ice$Migration <- ifelse(DayAD1_Ice$Start_lon >= -25 & DayAD1_Ice$End_lon <= -45 &  yday(DayAD1_Ice$date) < 185, 1, DayAD1_Ice$Migration)
table(DayAD1_Ice$Migration)

## now relabel the migration column
DayAD1_Ice$Migration <- ifelse(DayAD1_Ice$Migration == 1, "Migration", "Non-migration")

##create list of all bird IDs
IDbird <- levels(as.factor(DayAD1_Ice$ID)) 

##create empty data frame of just column headings from DAYAD1 dataset
DayADALL_Ice = DayAD1_Ice[FALSE,]

# Minimum daily displacement value (in KMs) 
###long threshold: less then -24 and greater then -49
MinDDisp <- 50
Max_lon <- -24
Min_lon <- -49

##Run AD method with function
## DayADf is the DayAD1 data set created above
## IDs is a list of the Tag_IDs in the data set
## col_head is an empty data set with the col headings from DayAD1
DayADALL_Ice <- AD_method_lon(DayADf = DayAD1_Ice, IDs= IDbird, col_head = DayADALL_Ice, 
                              MinDDisp = MinDDisp, Max_lon = Max_lon , Min_lon = Min_lon)




##
#### 4.3 Extract migration start and end ####
##

## This extract the start and end dates for each migration 
## First it checks for multiple migrations for each individual and removes duplciatees that don't make it to the final destination.

## Filter out the days that have been classified as migration
find_Ice <- filter(DayADALL_Ice, Migration == "Migration")
## calculate number of days in between each row
find_Ice$date <- as.Date(find_Ice$date)
find_Ice$day_diff <- as.numeric(find_Ice$date -lag(find_Ice$date))
find_Ice$change <- ifelse(find_Ice$day_diff == 1, 0, 1)
find_Ice$change[1] <- 1

## use cumsum function to label each incubation
find_Ice$MigID <- cumsum(find_Ice$change)


##now use dplyr to create summary of each putative incubation
All_migs_Ice <- find_Ice %>% 
                group_by(tag_year, MigID) %>% 
                summarise(Ice_dep_spring = min(date), 
                          Green_arrive = max(date), 
                          time_taken = as.numeric(((Green_arrive-Ice_dep_spring)+1)),
                          min_lon = min(End_lon),
                          mean_fixes = mean(n),
                          ID = max(ID),
                          year = max(year)) 

## Identify any tag-years with multiple putative incubation
Dups21 <- duplicated(All_migs_Ice$tag_year)

## extract tag_years with multiple incubation attempts
Dup_tags_ice <- All_migs_Ice[Dups21 == T, ]

## separate out the tags with multiple incubation and those without multiple incubation
Dup_migs_ice <- All_migs_Ice %>% filter(tag_year %in% c(Dup_tags_ice$tag_year))
Non_dups_ice <- All_migs_Ice %>% filter(!tag_year %in% c(Dup_tags_ice$tag_year))

## get unique tag_years in the data set
unique_dups_ice <- unique(Dup_migs_ice$tag_year)

## run loop to select migration that makes it to Greenland
for(qq in 1:length(unique_dups_ice)) {
  
  ## extract the tag year
  dup_TY_ice <- filter(Dup_migs_ice, tag_year == unique_dups_ice[qq])
  
  ## remove any that don't make it to Greenland
  dup_TY_ice <- filter(dup_TY_ice, min_lon  < -49)
  
  ## select the one with the largest max end lat
  min <- min(dup_TY_ice$min_lon)
  dup_TY_ice <- filter(dup_TY_ice, min_lon == min)
  
  ## join together all of the top incubation
  if(qq == 1){dup_all_ice <- dup_TY_ice}
  else{dup_all_ice <- rbind(dup_all_ice, dup_TY_ice)}
  
}


## bind together the none duplicated tag years and the daat set with the duplciates removed
AllIce_deps <- rbind(Non_dups_ice, dup_all_ice)


## remove the migID column and then remove tag years where we don't want the data!
AllIce_deps$MigID <- NULL 
AllIce_deps <- AllIce_deps %>% filter(!(tag_year == "WHIT22_2018")) ##filter out this tag year as only went to E Greenland

##Extract spring dates for Icerland departure and Greenland arrival
Ice_spring_dates <- AllIce_deps




##
#### 4.4 Check phenology classification ####
##

### PLOTTING A MAP IN LEAFLET WITH COLOURED MARKERS ###

##filter rough time period for plotting that doesn't crash
DayADALL_Ice_leaflet <- DayADALL_Ice %>% filter(yday >= 120 & yday < 150)

## Create this variable and use it as a pop up when clicking on individual points on the map
DayADALL_Ice_leaflet$ID_year <- paste0(DayADALL_Ice_leaflet$ID, "_", DayADALL_Ice_leaflet$year, "_", DayADALL_Ice_leaflet$date)

## Create color pallet for the map
myColors <- brewer.pal(5,"Accent")
## assign each factor in the ringinglocations column to one of the colors
factpal <- colorFactor(myColors, DayADALL_Ice_leaflet$Migration)

##using leaflet to plot the lat/longs and colour= to colour the markers
leaflet(DayADALL_Ice_leaflet) %>% 
  addTiles()%>%
  addCircleMarkers(DayADALL_Ice_leaflet$End_lon, DayADALL_Ice_leaflet$End_lat, 
                   color = ~factpal(Migration), fillOpacity = 2, stroke = FALSE, radius = 5, 
                   popup = ~ID_year) %>% 
  addLegend("topleft", pal = factpal, values = ~Migration,
            title = "Ringing Location",
            labFormat = labelFormat(prefix = ""),
            opacity = 1)



### Checking if any tags that reached Greenland don't have a departure date ###

## filter tags that got to Greenland
DayADALL_Ice_check <- DayADALL_Ice %>% filter(End_lon < -49)

## extract the tag_years that got there
check <- as.data.frame(unique(DayADALL_Ice_check$tag_year))
colnames(check)[1] <- "tag_year"

## join to departures
missing_years <- anti_join(Ice_spring_dates, check, by = "tag_year")
stopifnot(nrow(missing_years)==0) ## should be zero, any tags in this data set got to Greenland but do not have an associated departure date






##
#### 5. Absolute displacement (AD) function for latitudinal migration ####
##

##Turning the loop in the paper into a function that you can resuse
##now have a function that does AD calculation for migration with large latitudinal component
AD_method_lat <- function(DayADf, IDs, col_heads, MinDDisp, Max_lat, Min_lat) {
  # loop for each bird 
  for(z in 1:length(IDs)) { 
    DayAD <- filter(DayADf, ID == IDs[z])
    
    # onset and end of migration 
    for(i in 2:(nrow(DayAD)-1)) {
      if(DayAD$Migration[i] != "Migration")
      {next} # find start of unmistakable migratory segment 
      if(DayAD$Migration[i-1] == "Migration" & DayAD$Migration[i+1] == "Migration") 
      {next} # skip intermidate locations of unmistakable migratory segment 
      # find onset of autumn migration or onset of spring migration 
      if(DayAD$Migration[i-1] != "Migration") {
        j=i-1 
        for (j in j:1){ 
          if(DayAD$ID[j] !=
             DayAD$ID[j+1]) {break} 
          if((DayAD$Start_lat[j] > Min_lat & DayAD$Start_lat[j] < Max_lat)
             | DayAD$DDisp[j] > MinDDisp | DayAD$BDISP[j] > MinDDisp)
          {DayAD$Migration[j] <- "Migration" 
          next}
          if((DayAD$Start_lat[j] < Min_lat | DayAD$Start_lat[j] > Max_lat) 
             | DayAD$DDisp[j] <= MinDDisp) {break}}}
      # find end of autumn migration or end of spring migration 
      if(DayAD$Migration[i+1] != "Migration") { 
        j=i+1
        for (j in j:nrow(DayAD)) {
          if(DayAD$ID[j] !=
             DayAD$ID[j-1]) {break}
          if((DayAD$Start_lat[j] > Min_lat & DayAD$Start_lat[j] < Max_lat) 
             | DayAD$DDisp[j] > MinDDisp | DayAD$BDISP[j] > MinDDisp) 
          {DayAD$Migration[j] <- "Migration" 
          next} 
          if((DayAD$Start_lat[j] < Min_lat | DayAD$Start_lat[j] > Max_lat) 
             | DayAD$DDisp[j] <= MinDDisp) {break}}}} 
    col_heads <- rbind(col_heads, DayAD)}
  return(col_heads)
}






##
#### 6. Calculate phenology for UK-->Ice Spring migration ####
##

##
#### 6.1 Calculate daily travel distances ####
##

## change data frame for AD method
AD_UK <- All_GPS

## extract the date from the timestamp
AD_UK$date <- as.Date(AD_UK$timestamp) 
AD_UK$yday <- yday(AD_UK$timestamp)

## label if a GPS fix is in the unmistakable migratory segment (Over the Atlantic between UK and Iceland)
AD_UK$Migration <- ifelse(AD_UK$lat <= 62.5 & AD_UK$lat >= 60 & AD_UK$yday < 185 & AD_UK$lon > -25, 1, 0)

## Create data frame with first and last location for each tag each day
DayAD_UK <- AD_UK %>% 
            group_by(ID, date) %>% 
            summarise(Start_lat = first(lat), Start_lon = first(lon), 
                      End_lat = last(lat), End_lon = last(lon),n=n(),
                      Migration=sum(Migration)) 

## Change migration column to either 1 or zero depending if there were any fixes in the unmistakable migratory segment
DayAD_UK$Migration <- ifelse(DayAD_UK$Migration==0, 0, 1)

## Now calculate how far bird has gone been start and end of each day
DayAD_UK$DDisp<- distGeo(cbind(DayAD_UK$Start_lon, DayAD_UK$Start_lat), 
                       cbind(DayAD_UK$End_lon, DayAD_UK$End_lat))/1000


## Now calculate how far bird has gone been start of one day and end of previous
##create a tag_year column
DayAD_UK <- DayAD_UK %>% mutate(year=year(ymd(date)), tag_year = paste(ID, year, sep="_"))

##now create column with this new distance in
DayAD_UK <- DayAD_UK %>% 
            group_by(tag_year) %>% 
            nest() %>% 
            mutate(BDISP = map(data, Diff_day_dist_calc)) %>% 
            unnest(cols = c(data, BDISP))

DayAD_UK$BDISP <- ifelse(is.na(DayAD_UK$BDISP) == TRUE, 1, DayAD_UK$BDISP)

##reset the BDISP colum to 1 if the number of days between rows is greater than 1
DayAD_UK$DayDiff <- ymd(DayAD_UK$date) - ymd(lag(DayAD_UK$date))
DayAD_UK$BDISP <- ifelse(DayAD_UK$DayDiff > 2, 1, DayAD_UK$BDISP)




##
#### 6.2 Classify migrations ####
##

## First let's remove any  birds that didn't make it to Greenland in a particular year
MaxDists2 <- DayAD_UK %>% 
             group_by(tag_year) %>% 
             summarise(Maxlat = max(End_lat))
## filter out the ones that did make it, -45 is half way across the Greenalnd Icecap
MaxDists2 <- MaxDists2[MaxDists2$Maxlat > 63,]
## now filter out the birds that made it to Greenland 
survivors2 <- unique(MaxDists2$tag_year)
DayAD_UK <- filter(DayAD_UK, tag_year %in% survivors2)


## Already determined unmistakable migratory segment (e.g. over North Atlantic between UK and Iceland) 
## Add another condition in case birds crossed this segment in one day
table(DayAD_UK$Migration)
DayAD_UK$Migration <- ifelse(DayAD_UK$Start_lat <= 59 & DayAD_UK$End_lat >= 63 &  yday(DayAD_UK$date) < 185, 1, DayAD_UK$Migration)
table(DayAD_UK$Migration)

## now relabel the migration column
DayAD_UK$Migration <- ifelse(DayAD_UK$Migration == 1, "Migration", "Non-migration")

##create list of all bird IDs
IDbird_UK <- levels(as.factor(DayAD_UK$ID)) 
##create empty data frame of just column headings from DAYAD1 dataset
DayADALL_UK = DayAD_UK[FALSE,]

# Minimum daily displacement value (in KMs) 
###long threshold: less then -24 and greater then -49
MinDDisp <- 50
Max_lat <- 63
Min_lat <- 58.7

##Run AD method with function
## DayADf is the DayAD_UK data set created above
## IDs is a list of the Tag_IDs in the data set
## col_head is an empty data set with the col headings from DayAD_UK
DayADALL_UK <- AD_method_lat(DayADf = DayAD_UK, IDs= IDbird_UK, col_head = DayADALL_UK, 
                           MinDDisp = MinDDisp, Max_lat = Max_lat , Min_lat = Min_lat)





##
#### 6.3 Extract migration start and end ####
##

## This extract the start and end dates for each migration 
## First it checks for multiple migrations for each individual and removes duplciatees that don't make it to the final destination.


## Filter out the days that have been classified as migration
find_UK <- filter(DayADALL_UK, Migration == "Migration")
## calculate number of days in between each row
find_UK$date <- as.Date(find_UK$date)
find_UK$day_diff <- as.numeric(find_UK$date -lag(find_UK$date))
find_UK$change <- ifelse(find_UK$day_diff == 1, 0, 1)
find_UK$change[1] <- 1

## use cumsum function to label each incubation
find_UK$MigID <- cumsum(find_UK$change)

##now use dplyr to create summary of each putative incubation
All_migs <- find_UK %>% 
                group_by(tag_year, MigID) %>% 
                summarise(UK_dep = min(date), 
                          Ice_arrive_spring = max(date), 
                          time_taken = as.numeric(((Ice_arrive_spring-UK_dep)+1)),
                          max_lat = max(End_lat),
                          mean_fixes = mean(n),
                          ID = max(ID),
                          year = max(year)) 

## Identify any tag-years with multiple putative incubation
Dups22 <- duplicated(All_migs$tag_year)

## extract tag_years with multiple incubation attempts
Dup_tags <- All_migs[Dups22 == T, ]

## separate out the tags with multiple incubation and those without multiple incubation
Dup_migs <- All_migs %>% filter(tag_year %in% c(Dup_tags$tag_year))
Non_dups <- All_migs %>% filter(!tag_year %in% c(Dup_tags$tag_year))

## get unique tag_years in the data set
unique_dups <- unique(Dup_migs$tag_year)

## run loop to select migraiton that makes it to Iceland
for(qq in 1:length(unique_dups)) {
  
  ## extract the tag year
  dup_TY <- filter(Dup_migs, tag_year == unique_dups[qq])
  
  ## remove any that don't make it to Iceland
  dup_TY <- filter(dup_TY, max_lat  > 63)
  
  ## select the one with the largest max end lat
  max <- max(dup_TY$max_lat)
  dup_TY <- filter(dup_TY, max_lat == max)
  
  ## join together all of the top incubation
  if(qq == 1){dup_all <- dup_TY}
  else{dup_all <- rbind(dup_all, dup_TY)}
  
}

## bind together the none duplicated tag years and the daat set with the duplciates removed
AllUK_deps <- rbind(Non_dups, dup_all)


## remove the migID column and then remove tag years where we don't want the data!
AllUK_deps$MigID <- NULL 
AllUK_deps <- AllUK_deps %>% filter(!(tag_year == "WHIT22_2018")) ##filter out this tag year as only went to E Greenland

##Extract spring dates for Icerland departure and Greenland arrival
UK_spring_dates <- AllUK_deps




##
#### 6.4 Check phenology classification ####
##

###PLOTTING A MAP IN LEAFLET WITH COLOURED MARKERS###

##filter for faster plotting
DayADALL_UK_leaflet <- DayADALL_UK %>% filter(yday > 20 & yday <= 125)
DayADALL_UK_leaflet$ID_year <- paste(DayADALL_UK_leaflet$ID, DayADALL_UK_leaflet$year, sep="_")

##first create colour pallet with 10 colours from the paired pallet
myColors <- brewer.pal(5,"Accent")
##Next assign each factor in the ringinglocations column to one of the colors
factpal <- colorFactor(myColors, DayADALL_UK_leaflet$Migration)

##using leaflet to plot the lat/longs and colour= to colour the markers
leaflet(DayADALL_UK_leaflet) %>% 
  addTiles()%>%
  addCircleMarkers(DayADALL_UK_leaflet$End_lon, DayADALL_UK_leaflet$End_lat, 
                   color = ~factpal(Migration), fillOpacity = 2, stroke = FALSE, radius = 5, 
                   popup = ~ID_year) %>% 
  addLegend("topleft", pal = factpal, values = ~Migration,
            title = "Ringing Location",
            labFormat = labelFormat(prefix = ""),
            opacity = 1)


### Checking if any tags that reached Greenland don't have a departure date ###
## filter tags that got to Greenland
DayADALL_UK_check <- DayADALL_UK %>% filter(End_lat > 63)

## extract the tag_years that got there
check2 <- as.data.frame(unique(DayADALL_UK_check$tag_year))
colnames(check2)[1] <- "tag_year"

## join to departures
missing_years2 <- anti_join(UK_spring_dates, check2, by = "tag_year")
stopifnot(nrow(missing_years2)==0) ## should be zero, any tags in this data set got to Greenland but do not have an associated departure date












##
#### 7. Calculate phenology for Green-->Green Spring migration ####
##


##
#### 7.1 Calculate daily travel distances ####
##

## change data frame for AD method
AD_Green <- All_GPS 

## extract the date from the timestamp
AD_Green$date <- as.Date(AD_Green$timestamp) 
AD_Green$yday <- yday(AD_Green$timestamp) 

## label if a GPS fix is in the unmistakable migratory segment (Over the Atlantic between Iceland and Greenland)
## This time has to be after year day = 185 to get the spring instead of the autumn migraiton
AD_Green$Migration <- ifelse(AD_Green$lon <= -25 & AD_Green$lon >= -45 & AD_Green$lat >= -62.5 &  AD_Green$yday > 185, 1, 0)

## Create data frame with first and last location for each tag day
DayAD1_Green<- AD_Green %>% 
              group_by(ID, date) %>% 
              summarise(Start_lat = first(lat), Start_lon = first(lon), 
                        End_lat = last(lat), End_lon = last(lon),
                        n=n(), Migration=sum(Migration)) 

## Change migration column to either 1 or zero depending if there were any fixes in the unmistakable migratory segment
DayAD1_Green$Migration <- ifelse(DayAD1_Green$Migration==0, 0, 1)

## Now calculate how far bird has gone been the start and end of each day (DDISP)
DayAD1_Green$DDisp<- distGeo(cbind(DayAD1_Green$Start_lon, DayAD1_Green$Start_lat), 
                           cbind(DayAD1_Green$End_lon, DayAD1_Green$End_lat))/1000


## Now calculate how far bird has gone been start of one day and end of previous (BDISP)
## create a tag_year column
DayAD1_Green <- DayAD1_Green %>% mutate(year=year(ymd(date)), tag_year = paste(ID, year, sep="_"))

##function to calculate distance between start of one day and end of last
Diff_day_dist_calc <- function(x){
  
  distGeo(cbind(x$Start_lon, x$Start_lat), 
          cbind(lag(x$End_lon), lag(x$End_lat)))/1000
  
}

##now create column with this BDISP
DayAD1_Green <- DayAD1_Green %>% 
                group_by(tag_year) %>% 
                nest() %>% 
                mutate(BDISP = map(data, Diff_day_dist_calc)) %>% 
                unnest(cols = c(data, BDISP))

DayAD1_Green$BDISP <- ifelse(is.na(DayAD1_Green$BDISP) == TRUE, 1, DayAD1_Green$BDISP)

##reset the BDISP colum to 1 if the number of days between rows is greater than 1
DayAD1_Green$DayDiff <- ymd(DayAD1_Green$date) - ymd(lag(DayAD1_Green$date))
DayAD1_Green$BDISP <- ifelse(DayAD1_Green$DayDiff > 2, 1, DayAD1_Green$BDISP)




##
#### 7.2 Classify migrations ####
##


## First let's remove any  birds that didn't make it to Iceland in a particular year (after year day 185)
MaxDists <- DayAD1_Green %>% 
            mutate(yday = yday(date)) %>% 
            filter(yday > 185) %>% 
            group_by(tag_year) %>% 
            summarise(Maxlon = max(End_lon))

## filter out the ones that did make it back to Iceland after year day = 185
MaxDists <- MaxDists[MaxDists$Maxlon > -24,]
## now filter out the birds that made it to Greenland 
survivors <- unique(MaxDists$tag_year)
DayAD1_Green <- filter(DayAD1_Green, tag_year %in% survivors)


## Unmistakable migratory segment (e.g. over North Atlantic and Greenland Ice cap) was determined above
## Add another condition in case birds crossed this segment in one day
table(DayAD1_Green$Migration)
DayAD1_Green$Migration <- ifelse(DayAD1_Green$End_lon >= -25 & DayAD1_Green$Start_lon <= -45 &  yday(DayAD1_Green$date) > 185, 1, DayAD1_Green$Migration)
table(DayAD1_Green$Migration)

## now relabel the migration column
DayAD1_Green$Migration <- ifelse(DayAD1_Green$Migration == 1, "Migration", "Non-migration")

##create list of all bird IDs
IDbird <- levels(as.factor(DayAD1_Green$ID)) 

##create empty data frame of just column headings from DAYAD1 dataset
DayADALL_Green = DayAD1_Green[FALSE,]

# Minimum daily displacement value (in KMs) 
###long threshold: less then -24 and greater then -49
MinDDisp <- 50
Max_lon <- -24
Min_lon <- -49

##Run AD method with function
## DayADf is the DayAD1 data set created above
## IDs is a list of the Tag_IDs in the data set
## col_head is an empty data set with the col headings from DayAD1
DayADALL_Green <- AD_method_lon(DayADf = DayAD1_Green, IDs= IDbird, col_head = DayADALL_Green, 
                                MinDDisp = MinDDisp, Max_lon = Max_lon , Min_lon = Min_lon)




##
#### 7.3 Extract migration start and end ####
##

## This extract the start and end dates for each migration 
## First it checks for multiple migrations for each individual and removes duplciatees that don't make it to the final destination.

## Filter out the days that have been classified as migration
find_Green <- filter(DayADALL_Green, Migration == "Migration")
## calculate number of days in between each row
find_Green$date <- as.Date(find_Green$date)
find_Green$day_diff <- as.numeric(find_Green$date -lag(find_Green$date))
find_Green$change <- ifelse(find_Green$day_diff >= 1, 0, 1)
find_Green$change[1] <- 1

## use cumsum function to label each migration
find_Green$MigID <- cumsum(find_Green$change)


##now use dplyr to create summary of each putative migration
All_migs_Green <- find_Green %>% 
                  group_by(tag_year, MigID) %>% 
                  summarise(Green_dep_aut = min(date), 
                            Ice_arrive_aut = max(date), 
                            time_taken = as.numeric(((Ice_arrive_aut-Green_dep_aut)+1)),
                            max_lon = max(End_lon),
                            min_lon = min(Start_lon),
                            mean_fixes = mean(n),
                            ID = max(ID),
                            year = max(year)) 

## Identify any tag-years with multiple putative incubation
Dups21 <- duplicated(All_migs_Green$tag_year)

## extract tag_years with multiple incubation attempts
Dup_tags_Green <- All_migs_Green[Dups21 == T, ]

## separate out the tags with multiple incubation and those without multiple incubation
Dup_migs_Green <- All_migs_Green %>% filter(tag_year %in% c(Dup_tags_Green$tag_year))
Non_dups_Green <- All_migs_Green %>% filter(!tag_year %in% c(Dup_tags_Green$tag_year))
nrow(Dup_migs_Green);nrow(Non_dups_Green)

## get unique tag_years in the data set
unique_dups_Green <- unique(Dup_migs_Green$tag_year)

if(nrow(Dup_migs_Green)==0){print("No duplcicates")
                            AllGreen_deps <- Non_dups_Green}else{
                              
  ## run loop to select migration that makes it to Iceland
  for(qq in 1:length(unique_dups_Green)) {
    
    ## extract the tag year
    dup_TY_Green <- filter(Dup_migs_Green, tag_year == unique_dups_Green[qq])
    
    ## remove any that don't make it to Iceland and any that didnt at least start in Greenland
    dup_TY_Green <- filter(dup_TY_Green, max_lon  > -24)
    dup_TY_Green <- filter(dup_TY_Green, min_lon  < -45)
    
    ## select the one with the largest max end lon
    max <- max(dup_TY_Green$max_lon)
    dup_TY_Green <- filter(dup_TY_Green, max_lon == max)
    
    ## join together all of the top incubation
    if(qq == 1){dup_all_Green <- dup_TY_Green}
    else{dup_all_Green <- rbind(dup_all_Green, dup_TY_Green)}
    
  }
  
  
  ## bind together the none duplicated tag years and the daat set with the duplciates removed
  AllGreen_deps <- rbind(Non_dups_Green, dup_all)
  
}



## remove the migID column and then remove tag years where we don't want the data!
AllGreen_deps$MigID <- NULL 
AllGreen_deps <- AllGreen_deps %>% filter(!(tag_year == "WHIT22_2018")) ##filter out this tag year as only went to E Greenland

##Extract spring dates for Icerland departure and Greenland arrival
Green_aut_dates <- AllGreen_deps




##
#### 7.4 Check phenology classification ####
##

### PLOTTING A MAP IN LEAFLET WITH COLOURED MARKERS ###

##filter rough time period for plotting that doesn't crash
DayADALL_Green$yday <- yday(DayADALL_Green$date)
DayADALL_Green_leaflet <- DayADALL_Green %>% filter(yday > 185 & yday < 300)

## Create this variable and use it as a pop up when clicking on individual points on the map
DayADALL_Green_leaflet$ID_year <- paste0(DayADALL_Green_leaflet$ID, "_", DayADALL_Green_leaflet$year, "_", DayADALL_Green_leaflet$date)

## Create color pallet for the map
myColors <- brewer.pal(5,"Accent")
## assign each factor in the ringinglocations column to one of the colors
factpal <- colorFactor(myColors, DayADALL_Green_leaflet$Migration)

##using leaflet to plot the lat/longs and colour= to colour the markers
leaflet(DayADALL_Green_leaflet) %>% 
  addTiles()%>%
  addCircleMarkers(DayADALL_Green_leaflet$End_lon, DayADALL_Green_leaflet$End_lat, 
                   color = ~factpal(Migration), fillOpacity = 2, stroke = FALSE, radius = 5, 
                   popup = ~ID_year) %>% 
  addLegend("topleft", pal = factpal, values = ~Migration,
            title = "Ringing Location",
            labFormat = labelFormat(prefix = ""),
            opacity = 1)



### Checking if any tags that reached Greenland don't have a departure date ###

## filter tags that got to Greenland
DayADALL_Green_check <- DayADALL_Green %>% filter(End_lon > -24)

## extract the tag_years that got there
check3 <- as.data.frame(unique(DayADALL_Green_check$tag_year))
colnames(check3)[1] <- "tag_year"

## join to departures
missing_years3 <- anti_join(Green_aut_dates, check3, by = "tag_year")
stopifnot(nrow(missing_years3)==0) ## should be zero, any tags in this data set got to Greenland but do not have an associated departure date








##
#### 8. Combine all phonological dates ####
##

UK_spring <- UK_spring_dates %>% select(UK_dep, Ice_arrive_spring, year, ID)
Ice_spring <- Ice_spring_dates %>% select(Ice_dep_spring, Green_arrive, ID, year)
Green_autumn <- Green_aut_dates %>% select(Green_dep_aut, Ice_arrive_aut, ID, year)

Phenology <- full_join(UK_spring, Ice_spring, by = c("ID", "year", "tag_year")) %>% full_join(Green_autumn, by = c("ID", "year", "tag_year"))

## write out the file
write_csv(Phenology, file = "Outputs/Full_phenology_upto2021.csv")
