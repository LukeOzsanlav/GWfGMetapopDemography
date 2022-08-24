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
pacman::p_load(tidyverse, data.table, lubridate)




#---------------------------------------------#
#### 1. Read in Incubation data and format ####
#---------------------------------------------#

## read in Incubation attempts (GPS + Acc Tags)
Inc <- read.csv("Outputs/Incubation_attempt_lengths_new.csv")
Inc$UK_dep <- NULL

## set dates to as.date format
Inc$Ice_arrive_spring <- as.Date(Inc$Ice_arrive_spring)
Inc$Ice_dep_spring <- as.Date(Inc$Ice_dep_spring)


## change name
Inc_all <- Inc 


## create binary column for success/fail and attempt/defer
Inc_all$attempt <- ifelse(Inc_all$length == 0, 0, 1)
Inc_all$success24 <- ifelse(Inc_all$length >=24, 1, 
                            ifelse(Inc_all$length == 0, NA, 0))

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
#### 4. Add in Year centered dates ####
#-------------------------------------#

## Add year centered Greenland arrival date, laying data and failure data
Inc_ext4 <- Inc_ext2 %>% 
            drop_na(Green_arrive) %>% 
            mutate(Green_yday = yday(ymd(Green_arrive)),
                   attempt_yday =yday(ymd(attempt_start)),
                   end_yday =yday(ymd(attempt_end))) %>% 
            group_by(year) %>% 
            mutate(Green_centre = Green_yday-median(Green_yday),
                   laying_centre = attempt_yday-median(attempt_yday, na.rm = T),
                   failure_centre = end_yday-median(end_yday, na.rm = T)) %>% 
            ungroup() %>% 
            select(-c(Green_yday, attempt_yday, end_yday))



#-------------------------------------------------#
#### 5. Read out the final incubation data set ####
#-------------------------------------------------#

## read out this final incubation data set
write_csv(Inc_ext4, file = "Outputs/Incubation_attempts_for_models.csv")

