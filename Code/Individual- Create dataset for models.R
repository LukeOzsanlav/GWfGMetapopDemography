## Luke Ozsanlav-Harris

## Created:06/08/2020
## Updates 01/10/2020 to add in calculations for year centered phenology dates

## Create data set for all incubating females that can be used in later scripts for modeling
## Explanaotories should include:
##    Important departure/arrival and incubation dates
##    Binary columns for attempt/defer and success/fail
##    Pre breeding and incubation lengths
##    Ordinal breeding outcomes
##    Ringing location
##    Previous breeding success
##    Breeding lat/long
##    Migration length in time
##    Year centered phenological dates


## Packages required
pacman::p_load(tidyverse, data.table)




#---------------------------------------------#
#### 1. Read in Incubation data and format ####
#---------------------------------------------#

## read in Incubation attempts (GPS + Acc Tags)
#setwd("~/PhD Documents/1_Tracking Data Chapters/Incuabtion Lengths/Breeding attempts")
Inc <- read.csv("Outputs/Incubation_attempt_lengths.csv")
Inc$Acc <- "Y"
## set dates to as.date format
Inc$Ice_arrive_spring <- as.Date(Inc$Ice_arrive_spring)
Inc$Ice_dep_spring <- as.Date(Inc$Ice_dep_spring)

## read in Incubation attempts (GPS only Tags)
#setwd("~/PhD Documents/1_Tracking Data Chapters/Incuabtion Lengths/Breeding attempts")
Inc_gps <- read.csv("Outputs/Incubation_attempt_lengths_GPSonly.csv")
Inc_gps$Acc <- "N"
Inc_gps$X <- NULL
## set dates to as.date format
Inc_gps$Ice_arrive_spring <- as.Date(Inc_gps$Ice_arrive_spring)
Inc_gps$Ice_dep_spring <- as.Date(Inc_gps$Ice_dep_spring)


## bind the two togther
Inc_all <- rbind(Inc, Inc_gps) # bind so have dataset of all tags


## create binary column for success/fail and attempt/defer
Inc_all$attempt <- ifelse(Inc_all$length == 0, 0, 1)
Inc_all$success24 <- ifelse(Inc_all$length >=24, 1, 
                            ifelse(Inc_all$length == 0, NA, 0))

Inc_all$successall <- ifelse(Inc_all$length >=24, 1, 
                             ifelse(Inc_all$length == 0, 0, 0))

## create ordinal column with different breeding categories
Inc_all$ordinal <-  ifelse(Inc_all$length == 0, "defer", 
                           ifelse(Inc_all$length > 0 & Inc_all$length <= 11, "early.fail",
                                  ifelse(Inc_all$length > 11 & Inc_all$length <= 22, "late.fail", 
                                         ifelse(Inc_all$length > 22, "success", "error"))))

## create column for staging length
Inc_all$staging_length <- as.numeric(Inc_all$Ice_dep_spring - Inc_all$Ice_arrive_spring)




#------------------------------------------------------------------#
#### 2. Read in additional tag data and bind to incubation data ####
#------------------------------------------------------------------#

## importing dataset with additonal information on individuals 
#setwd("~/PhD Documents/Luke's Data")
tag_info <- read.csv("Data/Metadata/Tagged bird summary data new.csv")

## creating ID column to bind everthing by
tag_info$Bird.ID <- as.character(tag_info$Bird.ID)
tag_info$S.N <- as.character(tag_info$S.N)
tag_info$S.N <- ifelse(is.na(tag_info$S.N) == TRUE, tag_info$Bird.ID, tag_info$S.N)

## selecting important columns and sorting column names
setnames(tag_info, old = c("S.N"), new = c("ID"))
tag_info2 <- subset(tag_info, select = c("Age", "Mass", "Sex", "Ringing.location", "ID", "X2013.Gos", "X2014.Gos", 
                                         "X2016.Gos", "X2017.Gos","X2018.Gos", "X2019.Gos", "X2020.Gos"))

## changing HVAN to WEXF as all HVAN caugh birds wintered in Wexford
tag_info2$Ringing.location <- ifelse(tag_info2$Ringing.location == "HVAN", "WEXF", paste(tag_info2$Ringing.location))

###joining together Incubation attempt data, env data and tag info 
Inc_all$ID <- as.character(Inc_all$ID)
Inc_ext <- inner_join(Inc_all, tag_info2, by = "ID")
stopifnot(length(unique(Inc_all$Tag_year)) == length(unique(Inc_ext$Tag_year))) ## two lengths should be the same or have lost tag_years

## Create column for breeding in year - 1
Inc_ext <- Inc_ext %>% 
           mutate(across(c(X2013.Gos, X2014.Gos, X2016.Gos, X2017.Gos, X2018.Gos, X2019.Gos), as.numeric))
prev_success <- ifelse(Inc_ext$year == 2014 & Inc_ext$X2013.Gos > 0 |
                       Inc_ext$year == 2015 & Inc_ext$X2014.Gos > 0 |
                       Inc_ext$year == 2017 & Inc_ext$X2016.Gos > 0 |
                       Inc_ext$year == 2018 & Inc_ext$X2017.Gos > 0 |
                       Inc_ext$year == 2019 & Inc_ext$X2018.Gos > 0 |
                       Inc_ext$year == 2020 & Inc_ext$X2019.Gos > 0, 1, 0)

Inc_ext <- Inc_ext %>% select(-c(X2013.Gos, X2014.Gos, X2016.Gos, X2017.Gos, X2018.Gos, X2019.Gos, X2020.Gos)) # remove the breeding success by year columns
Inc_ext$prev_success <- prev_success # add to data set




#---------------------------------#
#### 3. Join in nest locations ####
#---------------------------------#

## read in nest locations from a previous script
nest_locs <- fread("Outputs/Nest_site_locations.csv")

## bind to the incubation dataset
nest_locs$Tag_year <- as.character(nest_locs$Tag_year)
Inc_ext2 <- full_join(Inc_ext, nest_locs, by = "Tag_year")




#-------------------------------------#
#### 4. Join on migratory duration ####
#-------------------------------------#

## read in migration ODBA data
#setwd("~/PhD Documents/1_Tracking Data Chapters/Migration Flight/Ice-Green Migration/Script Outputs")
ODBA <- fread("Outputs/Ice_Green_migration_duration_ODBA.csv")

## bind to the incubation dataset
ODBA <- subset(ODBA, select = c("tag_year", "duration", "total_ODBA"))
setnames(ODBA, old = c("duration", "tag_year"), new = c("mig_duration", "Tag_year"))
Inc_ext3 <- left_join(Inc_ext2, ODBA, by = "Tag_year")



## **** FIXED UP TO HERE ****##


########
## 5. ##
######## Add in Year centred dates

#### Add year centred Greenland arrival date ####

## Create Greenland arival data columns for modelling
## Change date into year day
Inc_ext3$Green_yday <- as.numeric(yday(as.character(Inc_ext3$Green_arrive)))

## center arrival date by each year
## calculate median arrival date for each year
median_dates <- Inc_ext3 %>% 
                dplyr::group_by(year) %>% 
                dplyr::summarise(median_arrival = floor(median(Green_yday, na.rm = T)))

## create list of unique year in the data set
Years <- unique(Inc_ext3$year)

## run loop to extract each year from incubation and median dates data sets and caculate
## a year centred Greenland arrival data
for(i in 1:length(Years)){
  
  birds = filter(Inc_ext3, year == Years[i])
  birds$Green_centre <- NA
  median = filter(median_dates, year == Years[i])
  
  for(j in 1:nrow(birds)){
    birds$Green_centre[j] = birds$Green_yday[j] - median$median_arrival
  }
  
  if(i == 1){Inc_ext4 <- birds} 
  else{Inc_ext4 <- rbind(Inc_ext4, birds)}
  
}




#### Add year centred laying date ####

## Creeat attempt start explanatories for models
## Change date into year day
Inc_ext4$attempt_yday <- as.numeric(yday(as.character(Inc_ext4$attempt_start)))

## centre laying date by each year
## calcualte median laying date for each year
median_dates <- Inc_ext4 %>% 
  group_by(year) %>% 
  summarise(median_laying = floor(median(attempt_yday, na.rm = T)))

## create list of unique year in the data set
Years <- unique(Inc_ext4$year)

## run loop to extract each year from incubation and median dates data sets and caculate
## a year centred Greenland arrival data
for(i in 1:length(Years)){
  
  birds = filter(Inc_ext4, year == Years[i])
  birds$laying_centre <- NA
  median = filter(median_dates, year == Years[i])
  
  for(j in 1:nrow(birds)){
    birds$laying_centre[j] = birds$attempt_yday[j] - median$median_laying
  }
  
  if(i == 1){Inc_ext5 <- birds} 
  else{Inc_ext5 <- rbind(Inc_ext5, birds)}
  
}




#### Add year centred failure date ####

## Creeat attempt start explanatories for models
## Change date into year day
Inc_ext5$end_yday <- as.numeric(yday(as.character(Inc_ext5$attempt_end)))

## centre laying date by each year
## calcualte median laying date for each year
median_dates <- Inc_ext5 %>% 
  group_by(year) %>% 
  summarise(median_end = floor(median(end_yday, na.rm = T)))

## create list of unique years in the data set
Years <- unique(Inc_ext5$year)

## run loop to extract each year from incubation and median dates data sets and caculate
## a year centred Greenland arrival data
for(i in 1:length(Years)){
  
  birds = filter(Inc_ext5, year == Years[i])
  birds$failure_centre <- NA
  median = filter(median_dates, year == Years[i])
  
  for(j in 1:nrow(birds)){
    birds$failure_centre[j] = birds$end_yday[j] - median$median_end
  }
  
  if(i == 1){Inc_ext6 <- birds} 
  else{Inc_ext6 <- rbind(Inc_ext6, birds)}
  
}

## Now remove the year day columns created
Inc_ext6$Green_yday <- NULL
Inc_ext6$attempt_yday <- NULL
Inc_ext6$end_yday <- NULL

stopifnot(nrow(Inc_all) == nrow(Inc_ext6)) # stop if we have lost birds from the original data set

########
## 6. ##
######## Read out the final incubation data set

## read out this final incubation data set
write.csv(Inc_ext6, file = "Outputs/Incubation_attempts_extra_update.csv", row.names = F)

