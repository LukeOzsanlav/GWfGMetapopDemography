## Luke Ozsanlav-Harris

## Calculate the nest locations for all breeding females
## These locations will be combined with the rest of the data for each female breeding season

## Created: 06/08/2020

## packages required
pacman::p_load(tidyverse, data.table, lubridate, move)



#---------------------------------------------#
#### 1. Read in Incubation data and format ####
#---------------------------------------------#

## read in Incubation attempts (GPS + Acc Tags)
Inc <- read.csv("Outputs/Incubation_attempt_lengths_new.csv")
Inc$Acc <- "Y"

## set date format
Inc <- Inc %>% mutate(across(c(attempt_end, Green_arrive, attempt_start, Ice_arrive_spring), ymd))

## **DELETE**
# ## read in Incubation attempts (GPS only Tags)
# Inc_gps <- read.csv("Outputs/Incubation_attempt_lengths_GPSonly.csv")
# Inc_gps$Acc <- "N"
# Inc_gps$X <- NULL
# 
# ## bind the two togther
# Inc_all <- rbind(Inc, Inc_gps) # bind so have dataset of all tags
## **DELETE**


##filter the birds that were known to incubate
Inc_br <- Inc %>% filter(is.na(attempt_start) == F)

## set dates to as.date format
Inc_br$attempt_start <- ymd_hms(paste0(Inc_br$attempt_start, "_", "12:00:00"))
Inc_br$attempt_end <- ymd_hms(paste0(Inc_br$attempt_end, "_", "12:00:00"))




#----------------------------------------------#
#### 2. Read in tracking data from Movebank ####
#----------------------------------------------#

## load in GPS tracking data 
GWF.df <- readRDS("Tracking data/Eco_Orn_All_MBformat.RDS")

## use luibridate to parse timestamp
GWF.df$UTC_datetime <- ymd_hms(GWF.df$UTC_datetime)

## create tag_year column
GWF.df$Tag_year <- paste0(GWF.df$device_id, "_", year(GWF.df$UTC_datetime))




#-----------------------------------#
#### 3. Calculate the nest sites ####
#-----------------------------------#

#### filter the GPS data to extract periods when birds were known to be incubating
## Join the two data sets
GWF.df <- inner_join(GWF.df, Inc_br, by = "Tag_year")

## filter the incubation windows
GWF.inc <- GWF.df %>% filter(UTC_datetime > attempt_start & UTC_datetime < attempt_end)

## calculate median lat and median long for each tag year
length(unique(GWF.df$Tag_year))
med_nest_locs <- GWF.inc %>% 
                 group_by(Tag_year) %>% 
                 summarise(breeding_lat = median(Latitude), 
                           breeding_long = median(Longitude)) 



#----------------------------------#
#### 4. Read out the nest sites ####
#----------------------------------#

## read out the nest locations
write.csv(med_nest_locs, file = "Outputs/Nest_site_locations.csv", row.names = F)



