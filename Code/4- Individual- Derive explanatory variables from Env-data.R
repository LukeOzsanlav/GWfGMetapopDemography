## Luke Ozsanlav-Harris

## Derive explanatory variables for models from the annotation of tacking data by the Movebank Env-data system. 
## Environmental vairbales annotated were Snow cover, NDVI, Temp and Precip. 
## Products used: MODIS Snow 500m Daily Terra NDSI Snow Cover, MODIS Land Vegetation Indices 500m 16d Terra NDVI,
##                ECMWF Interim Full Daily SFC Temperature (2 m above Ground) & ECMWF Interim Full Daily SFC-FC Total Precipitation

## Created: 28/05/2020
## Updated: 10/07/2020 added in temp and precip data
## Updated: 05/10/2020 resembled tracking data so that there is an equal number of fixes per day for all birds
## Updates: 23/08/2022 cleaned and added to repo

##packages required
pacman::p_load(tidyverse, data.table, lubridate, zoo)




## 1. Read in tracking data sets annotated with Env data
## 2. Combine and format Env data sets
## 3. Import migratory phenology data then filter greenland period
## 4. Resample the tracking data
## 5. Calculate rolling averages of Env variables over the breeding season for each day
## 6. Import breeding attempt data and create custom windows for each bird 
## 7. loop to calculate average Env variables over custom window lengths


## NEED to add in the new env data as for some reason my old Env data does not have all of the tag_years that i need!!!




#-------------------------------------------------------------#
#### 1. Read in tracking data sets annotated with Env data ####
#-------------------------------------------------------------#

## MODIS 500m daily Terra NDSI dataset
NDVI <- readRDS("Env data/NDVI.RDS")

## MODIS 500m 16day Terra NDVI dataset
Snow <- readRDS("Env data/Snow_cover.RDS")

## NCEP precip and temp annotation, will have to be handled seperately
GWF_clim <- fread("Env data/NCEP_precip_temp_anno_EXTRABITS.csv")

## Rename env data columns
setnames(NDVI, old=c("MODIS Land Vegetation Indices 500m 16d Terra NDVI"), new=c("NDVI")) # NDVI
setnames(Snow, old = c("MODIS Snow 500m Daily Terra NDSI Snow Cover"), new = c("snow")) # snow cover






#-------------------------------------------#
#### 2. Combine and format Env data sets ####
#-------------------------------------------#

## Format NDVI and Snow cover data set

## Combine the NDVI and snow cover data sets
NDVI$Snow_cover <- Snow$snow
Env <- NDVI

## set timestamp as date time object
Env$timestamp <- ymd_hms(Env$timestamp)

## create tag_year column
Env$Tag_year <- paste0(Env$`tag-local-identifier`, "_", year(Env$timestamp))

##changing NaN to NA in Env data coloumns
Env$Snow_cover <- ifelse(Env$Snow_cover == "NaN", NA, Env$Snow_cover)
Env$NDVI <- ifelse(Env$NDVI == "NaN", NA, Env$NDVI)

## set env data columns to numeric
Env$Snow_cover <- as.numeric(Env$Snow_cover); Env$NDVI <- as.numeric(Env$NDVI)



## Format temp and precip data set
GWF_clim$timestamp <- as.character(ymd_hms(GWF_clim$timestamp))

## convert temperatures in Kelvin to degrees centigrade
GWF_clim$temp_2mC <- GWF_clim$temp_2m - 273.15

## convert precip from m to mm
GWF_clim$precip_ratemm <- GWF_clim$precip_rate*1000

## Create a absolute value of precip colum as currently it is a rate
## first create a time dif column, number of minutes between current fix and the next
GWF_clim$timedif <- as.numeric(as.character(difftime(lead(GWF_clim$timestamp), GWF_clim$timestamp, unit = "secs")))

GWF_clim$timedif <- ifelse(!GWF_clim$Tag_year == lead(GWF_clim$Tag_year) | !GWF_clim$Tag_year == lag(GWF_clim$Tag_year),
                      NA, GWF_clim$timedif)
GWF_clim$precip_tot <- GWF_clim$precip_ratemm*GWF_clim$timedif
## units of precip_tot are now g/m^2

unique(GWF_clim$Tag_year)



#------------------------------------------------------------------#
#### 3. Import migratory phenology data then filter Greenland period
#------------------------------------------------------------------#

## read in phenology data
Phenology <- fread("Outputs/Full_phenology_upto2021.csv")

## set a tag_year column
Phenology$Tag_year <- paste(Phenology$ID, Phenology$year, sep = "_")

## set phenology columns to lubridate timestamps and give 1 day buffer before departure and after arrival date
Phenology$Green_dep2 <- (ymd_hms(paste0(Phenology$Green_dep, " 00:00:00"))) - 86400
Phenology$Green_arrive2 <- (ymd_hms(paste0(Phenology$Green_arrive, " 00:00:00"))) + 86400

## now join the two data sets
Env2 <- inner_join(Env, Phenology, by = "Tag_year")

## filter out the period in Greenland
GWF_Gr <- Env2 %>%  filter(timestamp > Green_arrive2 & timestamp < Green_dep2) 

## check env-data columns
summary(GWF_Gr$Snow_cover)
summary(GWF_Gr$NDVI)


#### filter env data to just tag_years used in analysis ####

## Read in the incubation data set
Inc <- read.csv("Outputs/Incubation_attempt_lengths_new.csv") # GPS+Acc tags

# Inc_gps <- read.csv("Outputs/Incubation_attempt_lengths_GPSonly.csv") # GPS only tags
# Inc_gps$X <- NULL
# ## bind the two together
# Inc <- rbind(Inc, Inc_gps)

## Join to env data sets
Inc <- subset(Inc, select = c(Tag_year, length))
GWF_GR <- inner_join(GWF_Gr, Inc, by = "Tag_year")
stopifnot( length(table(GWF_GR$Tag_year)) == length(table(Inc$Tag_year)) ) # check join worked properly
length(unique(GWF_GR$Tag_year))
## lost WHIT01_2018 but had very little data anyway so may as well get rid of





########
## 4. ##
######## Resample the tracking data

## Resample  tracking data so each Tag_year has a similar number of daily fixes to caculate average Env conditions
trackSubSampyear <- function(TD, dt=1, unit='days', resamp_times = resamp_times){
  
  TD <- TD[order(TD$timestamp),] 
  years = unique(TD$Tag_year)
  
  for(j in 1:length(years)){
    
    print(j)
    
    TD_sub = filter(TD, Tag_year == years[j]) # extract one Tag_year at a time
    
    # breakdown to datasets per bird
    unid = unique(TD_sub$Tag_year) 
    nrid = length(unid)
    TDall = list(nrid)  
    TDsubred = list(nrid)
    timestep = paste(dt,unit) # 1 day timestep to create list of consequetive days
    # create time sequence 
    # Start and end date for each Tag_year
    dati_start = min(TD_sub$date, na.rm= T)
    dati_end = max(TD_sub$date, na.rm= T)
    #create sequence of dates for each Tag_year
    datiseq = as.data.frame(seq(from=dati_start, to=dati_end, by=timestep)); datiseq$Join <- "Join"; colnames(datiseq)[1] <- "dates"
    # Create list of times that want to sub sample to
    times = as.data.frame(resamp_times); times$Join <- "Join"; colnames(times)[1] <- "times"
    # create list of date times
    times_data = full_join(datiseq, times)
    times_data$seq = as.POSIXct(paste(times_data$dates, times_data$times, sep = " "), format= "%Y-%m-%d %H:%M:%OS", tz= "GMT")
    datiseq2 = times_data$seq
    
    for (i in 1:nrid) { # slightly defunct now but code doesnt need to be any quicker
      Dtemp = TD_sub[TD_sub$Tag_year == unid[i],]
      idx = sapply(datiseq2, function(x) which.min( abs( difftime( Dtemp$timestamp, x, units='mins')))) # finds closest time in data to your created time series
      TDall[[i]] = Dtemp
      TDsubred[[i]] = unique(Dtemp[idx,]) # the function unique makes sure that the rows in Dtemp[idx,] are unique - so no duplicate points
    }
    
    TDsubred2 <- do.call("rbind", TDsubred)
    if(j == 1){Resamp_data <- TDsubred2}
    else{Resamp_data <- rbind(Resamp_data, TDsubred2)}
    
  }
  
  return(Resamp_data)
  
}



## preapre data for function, timestamp needs to be POSIX object and add year column
GWF_GR$timestamp <- as.POSIXct(GWF_GR$timestamp, format= "%Y-%m-%d %H:%M:%OS", tz= "GMT")
GWF_GR$year <- year(GWF_GR$timestamp)
GWF_GR$date <- as.Date(GWF_GR$timestamp)
GWF_GR$tag_date <- paste0(GWF_GR$`individual-local-identifier`, "_", GWF_GR$date)


## Filter out rows past the end of August, remove September and October rows and those not contiaing precip or temp
GWF_GR$month <- month(GWF_GR$timestamp)
GWF_GR <- filter(GWF_GR, month < 9)
table(GWF_GR$month) # shouldn't have month values greater than 9

## Remove rows with missing NDVI values
GWF_GRndvi <- filter(GWF_GR, is.na(NDVI) == F)

## remove rows with missing snow cover values
GWF_GRsnow <- filter(GWF_GR, is.na(Snow_cover) == F)

## apply function to full data set
Mig_thindvi <- trackSubSampyear(GWF_GRndvi, resamp_times = c("08:00:00", "12:00:00", "16:00:00"))
Mig_thinsnow <- trackSubSampyear(GWF_GRsnow, resamp_times = c("12:00:00"))

## NDVI
summary(Mig_thindvi$NDVI)
length(unique(GWF_GR$tag_date)); length(unique(Mig_thindvi$tag_date))

## Snow cover
summary(Mig_thinsnow$Snow_cover)
length(unique(GWF_GR$tag_date)); length(unique(Mig_thinsnow$tag_date))


## NOTES on the resampling process
## 1. All fixes were annotated with weather variables using the Movebank env data system
## 2. NA values were removed from each varables before resampling
## 3. for NDVI, temp and precip the closet fix to 8:00, 12:00 and 16:00 was used for each day
## 4. For snow cover the closest fix to 12:00 was used for each day due to the high proportion of NAs in the data set






########
## 5. ##
######## Calculate rolling averages of temp and precip over the breeding season

## Calculate Rolling avergaes for Temperature and precipitation throughout the breeding season
## Not all days have all three fixes but will just create an average value for each day first

## IMPROVEMENTS:
## Rolling averges are from daily averages so may want to create them with the RAW data instead using a loop, like in section 6

## set timestamp as lubridate object
GWF_clim$timestamp <- ymd_hms(GWF_clim$timestamp)

## summarise daily weather per tag day
Daily_weather <- GWF_clim %>% group_by(tag_date) %>% summarise(avg_temp = mean(x= temp_2mC, na.rm= T), 
                                                               sum_precip = sum(precip_tot, na.rm= T),
                                                               ID= head(ID, n= 1L))

##changing NaN to NA in Env data coloumns
Daily_weather$avg_temp <- ifelse(Daily_weather$avg_temp == "NaN", NA, Daily_weather$avg_temp)
Daily_weather$sum_precip <- ifelse(Daily_weather$sum_precip == "NaN", NA, Daily_weather$sum_precip)


## summarise rolling weather averages per tag day
## fucntion to use in nested data frame for rolling mean and rolling sum
roll_mean_temp <- function(x, width) {rollapply(x$avg_temp, width = width, FUN = mean, align = "right", fill = NA, na.rm= T)}
roll_sum_precip <- function(x, width) {rollapply(x$sum_precip, width = width, FUN = sum, align = "right", fill = NA, na.rm= T)}

## applying rolling mean/sum on nested data set with varying window widths
Rolling_weather <- Daily_weather %>% 
                   group_by(ID) %>%  
                   nest() %>% 
                   mutate(avg_temp2 = purrr::map(data, roll_mean_temp, width = 2), sum_precip2 = purrr::map(data, roll_sum_precip, width = 2),
                          avg_temp3 = purrr::map(data, roll_mean_temp, width = 3), sum_precip3 = purrr::map(data, roll_sum_precip, width = 3),
                          avg_temp4 = purrr::map(data, roll_mean_temp, width = 4), sum_precip4 = purrr::map(data, roll_sum_precip, width = 4),
                          avg_temp5 = purrr::map(data, roll_mean_temp, width = 5), sum_precip5 = purrr::map(data, roll_sum_precip, width = 5),
                          avg_temp10 = purrr::map(data, roll_mean_temp, width = 10), sum_precip10 = purrr::map(data, roll_sum_precip, width = 10)) %>% 
                   unnest(cols = c(data, avg_temp2, sum_precip2, avg_temp3, sum_precip3, avg_temp4, sum_precip4, avg_temp5, sum_precip5, avg_temp10, sum_precip10))

## write out this file for use in other scripts
write_csv(Rolling_weather, file = "Outputs/Rolling_Env_data_per_day_breeding_season_evenfixes.csv")



#-------------------------------------------------------------------------------#
#### 6. Import breeding attempt data and create custom windows for each bird ####
#-------------------------------------------------------------------------------#

## read in Incubation attempts for GPS only and GPS + Acc tags
Inc <- read.csv("Outputs/Incubation_attempt_lengths_new.csv") # GPS+Acc tags
## set as date objects
Inc$attempt_end <- as.Date(Inc$attempt_end, format = "%Y-%m-%d")
Inc$Green_arrive <- as.Date(Inc$Green_arrive, format = "%Y-%m-%d")
Inc$attempt_start <- as.Date(Inc$attempt_start, format = "%Y-%m-%d")

## **DELETE**
# Inc_gps <- read.csv("Outputs/Incubation_attempt_lengths_GPSonly.csv") # GPS only tags
# Inc_gps$X <- NULL
# 
# ## bind the two together
# Inc <- rbind(Inc, Inc_gps)
## **DELETE**

## create some extra columns for periods to calculate averge env conditions over
Inc$Green_arrive10 <- ifelse(is.na(Inc$Green_arrive) == F, as.Date(Inc$Green_arrive) + 10, NA)
Inc$Green_arrive10 <- as.Date(Inc$Green_arrive10, origin = "1970-01-01")

Inc$Green_arrive20 <- ifelse(is.na(Inc$Green_arrive) == F, as.Date(Inc$Green_arrive) + 20, NA)
Inc$Green_arrive20 <- as.Date(Inc$Green_arrive20, origin = "1970-01-01")

Inc$Green_arrive30 <- ifelse(is.na(Inc$Green_arrive) == F, as.Date(Inc$Green_arrive) + 30, NA)
Inc$Green_arrive30 <- as.Date(Inc$Green_arrive30, origin = "1970-01-01")






#---------------------------------------------------------------------------------------------------------------#
#### 7. loop to calculate average Env variables over custom window lengths (first 10 & 20 days post arrival) ####
#---------------------------------------------------------------------------------------------------------------#

##---- Average TEMP and PRECIP ----##

## Create weights column so summing of temp isnt effected by days with only 1 or two fixes
daily_fixestemp <- as.data.frame(table(GWF_clim$tag_date))
colnames(daily_fixestemp)[1] <- "tag_date"
GWF_clim <- full_join(GWF_clim, daily_fixestemp, by = "tag_date")
GWF_clim$weight <- 3/GWF_clim$Freq

## parse timestamp with lubridate
GWF_clim$timestamp <- ymd_hms(GWF_clim$timestamp)
#GWF_Gr$timestamp <- as.POSIXct(GWF_Gr$timestamp, format= "%Y-%m-%d %H:%M:%S", tz= "GMT")

## create empty data set
Env_clim <- subset(Inc, select= c("Tag_year"))
Env_clim$Gr10temp <- NA
Env_clim$Gr20temp <- NA
Env_clim$Gr10precip <- NA
Env_clim$Gr20precip <- NA

## functions for use within the loop
mean_temp <- function(x){weighted.mean(x= x$temp_2mC, w= x$weight, na.rm= T)}
sum_precip <- function(x){sum(x$precip_tot, na.rm = T)}

## create list of unique tag years in Inc data set 
tagyears <- unique(Env_clim$Tag_year)
stopifnot(length(tagyears) == nrow(Env_clim))

## loop to filter out tag_years in the Inc data set then calculate ave conditions over various time windows
for(i in 1:length(tagyears)) {
  
  ## filter out the data and phenological data for the current tag_year
  sub_data = filter(GWF_clim, Tag_year == paste(tagyears[i]))
  Dates = filter(Inc, Tag_year == paste(tagyears[i]))
  
  message(i, " out of ", length(tagyears), " tags")
  
     ## calculate mean env varibales over various window lengths in relation to Greenland arrival
  
     Gr10 = sub_data %>% filter(timestamp <= ymd(Dates$Green_arrive10) & timestamp >= ymd(Dates$Green_arrive))
     Env_clim$Gr10temp[i] <- mean_temp(Gr10); Env_clim$Gr10precip[i] <- sum_precip(Gr10)
      
     Gr20 = sub_data %>% filter(timestamp <= Dates$Green_arrive20 & timestamp >= Dates$Green_arrive)
     Env_clim$Gr20temp[i] <- mean_temp(Gr20); Env_clim$Gr20precip[i] <- sum_precip(Gr20)
  
}

## remove the variables created in the loop from memory
rm(Dates, sub_data, Gr10, Gr20)





##---- Average NDVI ----##

## Create weights column so summing of NDVI isnt effected by days with only 1 or two fixes
daily_fixesndvi <- as.data.frame(table(Mig_thindvi$tag_date))
colnames(daily_fixesndvi)[1] <- "tag_date"
Mig_thindvi <- full_join(Mig_thindvi, daily_fixesndvi, by = "tag_date")
Mig_thindvi$weight <- 3/Mig_thindvi$Freq

## parse timestamp with lubridate
Mig_thindvi$timestamp <- ymd_hms(Mig_thindvi$timestamp)
#GWF_Gr$timestamp <- as.POSIXct(GWF_Gr$timestamp, format= "%Y-%m-%d %H:%M:%S", tz= "GMT")

## create empty data set
Env_ndvi <- subset(Inc, select= c("Tag_year"))
Env_ndvi$Gr10NDVI <- NA
Env_ndvi$Gr20NDVI <- NA

## functions for use within the loop
mean_NDVI <- function(x){weighted.mean(x= x$NDVI, w= x$weight, na.rm= T)}

## create list of unique tag years in Inc data set 
tagyears <- unique(Env_ndvi$Tag_year)
stopifnot(length(tagyears) == nrow(Env_ndvi))

## loop to filter out tag_years in the Inc data set then calculate ave conditions over various time windows
for(i in 1:length(tagyears)) {
  
  ## filter out the data and phenological data for the current tag_year
  sub_data = filter(Mig_thindvi, Tag_year == paste(tagyears[i]))
  Dates = filter(Inc, Tag_year == paste(tagyears[i]))
  
  message(i, " out of ", length(tagyears), " tags")
  
  ## calculate mean env varibales over various window lengths in relation to Greenland arrival
  
  Gr10 = sub_data %>% filter(timestamp <= ymd(Dates$Green_arrive10) & timestamp >= ymd(Dates$Green_arrive))
  Env_ndvi$Gr10NDVI[i] <- mean_NDVI(Gr10)
  
  Gr20 = sub_data %>% filter(timestamp <= Dates$Green_arrive20 & timestamp >= Dates$Green_arrive)
  Env_ndvi$Gr20NDVI[i] <- mean_NDVI(Gr20)
  
}

## remove the variables created in the loop from memory
rm(Dates, sub_data, Gr10, Gr20)





##---- Average SNOW COVER ----##

## Create weights column so summing of precip isnt effected by days with only 1 or two fixes
daily_fixessnow <- as.data.frame(table(Mig_thinsnow$tag_date))
colnames(daily_fixessnow)[1] <- "tag_date"
Mig_thinsnow <- full_join(Mig_thinsnow, daily_fixessnow, by = "tag_date")
Mig_thinsnow$weight <- 3/Mig_thinsnow$Freq

## parse timestamp with lubridate
Mig_thinsnow$timestamp <- ymd_hms(Mig_thinsnow$timestamp)
#GWF_Gr$timestamp <- as.POSIXct(GWF_Gr$timestamp, format= "%Y-%m-%d %H:%M:%S", tz= "GMT")

## create empty data set
Env_snow <- subset(Inc, select= c("Tag_year"))
Env_snow$Gr10Sn <- NA
Env_snow$Gr20Sn <- NA

## functions for use within the loop
mean_snow <- function(x){mean(x= x$Snow_cover, na.rm= T)}

## create list of unique tag years in Inc data set 
tagyears <- unique(Env_snow$Tag_year)
stopifnot(length(tagyears) == nrow(Env_snow))

## loop to filter out tag_years in the Inc data set then calculate ave conditions over various time windows
for(i in 1:length(tagyears)) {
  
  ## filter out the data and phenological data for the current tag_year
  sub_data = filter(Mig_thinsnow, Tag_year == paste(tagyears[i]))
  Dates = filter(Inc, Tag_year == paste(tagyears[i]))
  
  message(i, " out of ", length(tagyears), " tags")
  
  ## calculate mean env varibales over various window lengths in relation to Greenland arrival
  
  Gr10 = sub_data %>% filter(timestamp <= ymd(Dates$Green_arrive10) & timestamp >= ymd(Dates$Green_arrive))
  Env_snow$Gr10Sn[i] <- mean_snow(Gr10)
  
  Gr20 = sub_data %>% filter(timestamp <= Dates$Green_arrive20 & timestamp >= Dates$Green_arrive)
  Env_snow$Gr20Sn[i] <- mean_snow(Gr20)
  
}

## remove the variables created in the loop from memory
rm(Dates, sub_data, Gr10, Gr20)




##---- Join together all variable ----##

Env_all <- full_join(Env_clim, Env_ndvi, by = "Tag_year") %>% full_join(Env_snow, by = "Tag_year")
Env_all <- filter(Env_all, !Tag_year == "WHIT01_2018")

## read out this data set for use in other scripts
write_csv(Env_all, file = "Outputs/Average_Env_data_from_arrival_windows_equalfixes.csv")



