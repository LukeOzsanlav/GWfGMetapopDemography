## Luke Ozsanlav-Harris

## Create a time series for each of the unique Islay and Wexford nesting sites
## This time series will be from the start of May till the end of August for 1983 to 2019

## These data will then be annotated with temp and precip values from ECMWF on Movebank
## The following products needed to be annotated with Movebank:
    # 1. ECMWF Interim Full Daily SFC Temperature (2 m above Ground)
    # 2. ECMWF Interim Full Daily SFC-FC Total Precipitation


## Packages required
pacman::p_load(tidyverse, data.table, rnaturalearth, lubridate, svMisc)




#-------------------------------#
#### 1. Read in nesting data ####
#-------------------------------#

## Read in the Nesting data
Nests <- fread("Outputs/Incubation_attempts_extra_update.csv")

## filter out the Islay and Wexford birds
Nests <- Nests %>% dplyr::filter(is.na(attempt_start) == F & Ringing.location %in% c("ISLA", "WEXF"))

## filter out repeat individuals
index <- duplicated(Nests$ID)
Nest_solo <- dplyr::filter(Nests, index == F & Acc == "Y")




#------------------------------------------------------------------------------#
#### 2. Extract sites to sample climate at for Wexford and Islay population ####
#------------------------------------------------------------------------------#

## Extract all the nest sites for Wexford birds with repeats of individual removed
Wexf_clust <- dplyr::filter(Nest_solo, Ringing.location == "WEXF")

## Extract all the nest sites for Islay birds with repeats of individual removed
Islay_clust <- dplyr::filter(Nest_solo, Ringing.location == "ISLA")

## Join together the Islay and Wexford point and then plot them
Sample_points <- rbind(Wexf_clust, Islay_clust)





#--------------------------------------------#
#### 3. Plot the nest sites on a base map ####
#--------------------------------------------#

## create a base map that can be used for the plotting of locations
## Get a basemap of the world
world <- ne_countries(scale = "medium", returnclass = "sf")

## set bounds of the basemap
Base_map <- ggplot() +
            geom_sf(data = world, fill="lightblue",col="black") +
            coord_sf(xlim = c(min(Nest_solo$breeding_long)-1, max(Nest_solo$breeding_long)+1),
                     ylim = c(min(Nest_solo$breeding_lat)-1, max(Nest_solo$breeding_lat)+1), expand = FALSE) + 
            theme_bw()

## plot the locations on a basemap
Base_map + geom_point(data = Sample_points, aes(x= breeding_long, y =breeding_lat, colour = Ringing.location), size=3, alpha=1)





#-------------------------------------------------------------------------#
### 4. Create the list of dates for every breeding period back to 1983 ####
#-------------------------------------------------------------------------#

## create  lubridate dates for the start and end times
start_date <- ymd_hms("2020-04-01 00:00:01")
end_date <- ymd_hms("2020-08-31 00:00:01")

#calculate how many days in this time interval
n_starts <- interval(start_date,end_date)/hours(6)

## create sequence of dates
dates_2019 <- as.data.frame(start_date + hours(seq(from = 0, to = (6*n_starts), by = 6)))
colnames(dates_2019)[1] <- "Date"
dates_2019$Date <- ymd_hms(dates_2019$Date)

## Create the sequence of dates back to 1983 now
## Number of seconds in one year and one day
n_years <- 2020-1983 # number of years to go back for
p <- 0 # counter of leap years
d <- 60*60*24

## loop though the data set and minus off one year at a time while taking into account leap years
for(j in 1:n_years){
  
  svMisc::progress(j, max.value = n_years)
  
  dates_2019.2 <- dates_2019
  test <- leap_year(ymd_hms(dates_2019.2$Date[1])- (365*d*(j-1)))
  
  ## minus of number of years while correcting for leap years
  if(test == F){ dates_2019.2$Date <- ymd_hms(dates_2019.2$Date) - ((365*d*j) + (d*p)) }
  else{ p <- p + 1
  dates_2019.2$Date <- ymd_hms(dates_2019.2$Date) - ((365*d*j)+(d*p)) 
  }
  
  if(j == 1){Dates_all <- dates_2019.2}
  else{Dates_all <- rbind(Dates_all, dates_2019.2)}
  
}





#----------------------------------------------------------------------#
#### 5. Join together the sampling locations with the list of dates ####
#----------------------------------------------------------------------#

## Now add each of the locations to all of the dates
for(j in 1:nrow(Sample_points)){
  
  svMisc::progress(j, max.value = nrow(Sample_points))
  
  Dates_all2 <- Dates_all
  Dates_all2$Lat <- Sample_points$breeding_lat[j]
  Dates_all2$Long <- Sample_points$breeding_long[j]
  Dates_all2$Source_ID <- Sample_points$Tag_year[j]
  Dates_all2$Ringing_loc <- Sample_points$Ringing.location[j]
  
  if(j ==1){Date_locs <- Dates_all2}
  else{Date_locs <- rbind(Date_locs,Dates_all2)}
  
  
}




#------------------------------------------------#
#### 6. Prepare data for import into Movebank ####
#------------------------------------------------#

## rename columns so movebank recognises them
setnames(Date_locs, old = c("Date", "Lat", "Long"), new = c("timestamp", "location-lat", "location-long"))

## write out data set
write_csv(Date_locs, file = "Outputs/All_sampling_sites.csv")




