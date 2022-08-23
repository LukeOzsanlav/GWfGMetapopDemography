##Luke Ozsanlav-Harris

## Run incubation classifier when user has GPS and accelerometer data
## This approach requires some way of creating a training data set from known breeders
## Alternatively users can use the values from GWfG if working on a similar species. 


##loading Packages
pacman::p_load(tidyverse, lubridate, data.table, zoo, geosphere, ggpubr, DescTools)


#--------------------#
#### 1. Data Prep ####
#--------------------#


#-----------------------------------------#
#### 1.1 Importing/preparing ODBA data ####
#-----------------------------------------#

## read in ODBA data
ODBA <- readRDS("Tracking data/All_orni_ODBAQuantile.RDS")

## Format datetime, add tag_year column and select just the columns required for analysis
ODBA <- ODBA %>% 
        as.data.frame() %>% 
        mutate(device_id = as.character(device_id),
               UTC_datetime= ymd_hms(UTC_datetime),
               tag_year = paste0(device_id, "_", year(UTC_datetime))) %>% 
        select(device_id, UTC_datetime, burst_no_ID, ODBA_quant, tag_year) %>% 
        rename(ODBA = ODBA_quant)

## read in Ecotone ODBA
ODBAEco <- readRDS("Tracking data/Eco_ODBA_per_burst.RDS")

## Format datetime, add tag_year column and select just the columns required for analysis
ODBAEco <- ODBAEco %>% 
            as.data.frame() %>% 
            mutate(Tag_ID = as.character(Tag_ID),
                   UTC_datetime= ymd_hms(UTC_datetime),
                   tag_year = paste0(Tag_ID, "_", year(UTC_datetime))) %>% 
            select(Tag_ID, UTC_datetime, burst_no_ID, ODBA, tag_year) %>% 
            rename(device_id = Tag_ID)

## Now bind the two together
ODBA <- rbind(ODBA, ODBAEco)
rm(ODBAEco)


#----------------------------------------#
#### 1.2 Importing/preparing GPS data ####
#----------------------------------------#

## Read in Ornitela GPS data
orn_gps <- readRDS("Tracking data/Eco_Orn_All_MBformat.RDS")

## Format datetime, add tag_year column and select just the columns required for analysis
orn_gps <- orn_gps %>% 
           mutate(device_id = as.character(device_id),
                  UTC_datetime= ymd_hms(UTC_datetime),
                  tag_year = paste0(device_id, "_", year(UTC_datetime))) %>% 
           select(device_id, UTC_datetime, Latitude, Longitude, tag_year) %>% 
           filter(tag_year %in% unique(ODBA$tag_year))

## chekc that they have pretty much the same number of years
length(unique(orn_gps$tag_year))
length(unique(ODBA$tag_year))


#----------------------------------------------#
#### 1.3 Importing/preparing Phenology data ####
#----------------------------------------------#

## read in phenology data
Phenology <- fread("Outputs/Full_phenology_upto2021.csv")
glimpse(Phenology)

##set a tag_year column, parse greenland arrival/departure and add a day, select column needed
Phen <- Phenology %>% 
        drop_na(Green_dep_aut, Green_arrive) %>% 
        mutate(tag_year = paste0(ID, "_", year),
               Green_dep_aut = (ymd(Green_dep_aut) + days(1)),
               Green_arrive = (ymd(Green_arrive) + days(1))) %>% 
        select(tag_year, Green_dep_aut, Green_arrive)
  



#---------------------------------------------#
#### 1.4 Importing/preparing breeding data ####
#---------------------------------------------#

## importing dataset with breeding histories and adding info Inca data
tag_info <- read.csv("Data/Metadata/Tagged bird summary data new.csv")

##creating ID column to bind everything by
tag_info <- tag_info %>% mutate(across(c(Bird.ID, S.N), as.character))
tag_info$S.N <- ifelse(is.na(tag_info$S.N) == TRUE, tag_info$Bird.ID, tag_info$S.N)

##selecting important columns
tag_info <- tag_info %>% 
             rename(device_id = S.N) %>% 
             select(Age, Sex, Ringing.location, device_id, X2017.Gos, X2018.Gos, X2019.Gos, X2020.Gos)




#---------------------------------------------------#
#### 1.5 Trim biologging data to breeding period ####
#---------------------------------------------------#

## Trim the biologging data to just the period birds are on the breeding ground
## This will speed up computation as the data sets will be much smaller
## NOTE: this may not be possible in non-migrants but there are periods for many species when breeding does not occur


## Trim ODBA data ##

## remove tag_years without phenology dates first
ODBA <- ODBA %>% filter(tag_year %in% unique(Phen$tag_year))

## Join on the phenology data set and filter for Greenland breeding period
ODBA <- ODBA %>% 
        left_join(Phen, by = "tag_year") %>% 
        mutate(Date = as.Date(UTC_datetime)) %>% 
        filter(Date >= Green_arrive & Date <= Green_dep_aut) %>% 
        select(!c(Green_arrive, Green_dep_aut, Date))


## Trim GPS data ##

## remove tag_years without phenology dates first
orn_gps <- orn_gps %>% filter(tag_year %in% unique(Phen$tag_year))

## Join on the phenology data set and filter for Greenland breeding period
orn_gps <- orn_gps %>% 
           left_join(Phen, by = "tag_year") %>% 
           mutate(Date = as.Date(UTC_datetime)) %>% 
           filter(Date >= Green_arrive & Date <= Green_dep_aut) %>% 
           select(!c(Green_arrive, Green_dep_aut, Date))





#--------------------------#
#### 2 Data re-sampling ####
#--------------------------#


#----------------------------#
#### 2.1 OBDA re-sampling ####
#----------------------------#

## **This section is optional**
## Because we derive average ODBA the number of ODBA burst in each day does'nt effect the results too much
## I have left it in here in case users wanted to re-sample their ODBA so it matched one of the sampling rates I provided with training values


## TD = the tracking data set
## dt = the new sampling rate required in minutes

## Define function
trackSubSampyear <- function(TD = TD, dt=dt, unit=unit){
  
  TD <- TD[order(TD$UTC_datetime),] 
  years = unique(TD$year)
  
  for(j in 1:length(years)){
    
    message(j, " out of ", length(years), "years")
    
    TD_sub = filter(TD, year == years[j]) # extract one year at a time
    
    # breakdown to datasets per bird
    unid = unique(TD_sub$tag_year) 
    nrid = length(unid)
    TDall = list(nrid)  
    TDsubred = list(nrid)
    timestep = paste(dt,unit)
    # create time sequence from min to max time with step sizes that are defined at the start of the function
    dati_start = min(TD_sub$UTC_datetime, na.rm= T)
    dati_end = max(TD_sub$UTC_datetime, na.rm= T)
    datiseq = seq(from=dati_start, to=dati_end, by=timestep)
    
    for (i in 1:nrid) {
      Dtemp = TD_sub[TD_sub$tag_year == unid[i],]
      idx = sapply(datiseq, function(x) which.min( abs( difftime( Dtemp$UTC_datetime, x, units='mins')))) # finds closest time in data to your created time series
      TDall[[i]] = Dtemp
      TDsubred[[i]] = unique(Dtemp[idx,]) # the function unique makes sure that the rows in Dtemp[idx,] are unique - so no duplicate points
    }
    
    TDsubred2 <- do.call("rbind", TDsubred)
    if(j == 1){Resamp <- TDsubred2}
    else{Resamp <- rbind(Resamp, TDsubred2)}
    
  }
  
  return(Resamp)
  
}


## re-sample at desired sampling interval
ODBA <- ODBA %>% mutate(year = year(UTC_datetime))
#ODBA <- trackSubSampyear(TD = ODBA, dt = 6, unit = 'mins') # I won't run this now as it won't save any time




#---------------------------#
#### 2.2 GPS re-sampling ####
#---------------------------#

## **This section is not as optional**
## The classification likely works better with a more even sampling rate but it is not required
## I have left it in here in case users wanted to re-sample their GPS so it matched one of the sampling rates I provided with training values



## TD = the tracking data set
## dt = the new sampling rate required in minutes

## re-sample at desired sampling interval
orn_gps <- orn_gps %>% mutate(year = year(UTC_datetime))
#orn_gps <- trackSubSampyear(TD = orn_gps, dt = 15, unit = 'mins')





#-----------------------------------#
#### 3 Summarize Biologging Data ####
#-----------------------------------#


#------------------------------------------#
#### 3.1 Summarizing average daily ODBA ####
#------------------------------------------#

## Average daily ODBA for each tag 
ODBASum <- ODBA %>% 
            mutate(Date = as.Date(UTC_datetime)) %>% 
            group_by(device_id, Date, tag_year) %>% 
            summarise(avg_odba = mean(ODBA, na.rm=T))
  
## filtering out tags that don't have data for the full breeding season 
min_days <- 30 # set minimum number of days
ODBASum <- as.data.frame(table(ODBASum$tag_year)) %>% 
           filter(Freq > min_days) %>% 
           rename(tag_year = Var1) %>% 
           inner_join(ODBASum, by = "tag_year") %>% 
           select(-Freq)

  

#-----------------------------------------------#
#### 3.2 Calculate Movement metrics from GPS ####
#-----------------------------------------------#

## Going to make three movement metrics 
## 1. Daily Median NSD
## 2. Distance between median daily locations
## 3. Average daily distance between each location and each the median daily locations


## add on a daily nsd column
orn_gps <- orn_gps %>% 
           mutate(Date = as.Date(UTC_datetime)) %>% 
           group_by(device_id, Date) %>% 
           mutate(nsd_ = (distGeo(cbind(Longitude[1], Latitude[1]), 
                                  cbind(Longitude, Latitude))/1000)) %>% 
           ungroup()
  

## summarize by each tag_year
GPSSum <- orn_gps %>% 
          group_by(device_id, tag_year, Date) %>% 
          summarise(med_long = median(Longitude),
                    med_lat = median(Longitude),
                    nsd_ = median(nsd_))


## calculate distance between daily median locations remove unnecessary columns
GPSSum <- GPSSum %>% 
          group_by(tag_year) %>% 
          mutate(ddist = (distGeo(cbind(med_long, med_lat), 
                                  cbind(lag(med_long), lag(med_lat)))/1000)) %>% 
          ungroup() %>% 
          select(!c(med_long, med_lat))


# ##Added: 02/06/2020
# ## Created my own metric of movement
# ## For each day it calculates a median lat and a long
# ## then within each tag day it calculates the distance between a fix and
# ## the median location. These distances are then averaged over each tag_day
# 
# ## turn track into data frame
# track <- as.data.frame(trk)
# 
# ##first work out median lat/long for each tag day
# track <- track %>% 
#   group_by(day_ID) %>% 
#   summarise(med_long = median(x_), med_lat = median(y_)) %>% 
#   full_join(track, by = "day_ID") ##join back to original data set, each tag day has own med lat/long
# 
# ## calculate distance between fix and med lat/long for that tag day
# track$distance <- distGeo(cbind(track$x_, track$y_), 
#                           cbind(track$med_long, track$med_lat))/1000
# 
# ## calculate the average distance to the daily median lat/long
# ave_dist_med_loc <- track %>%
#   group_by(id, date) %>%
#   summarise(ave_dist = mean(distance), n_fixes = length(x_)) %>% 
#   mutate(year= year(ymd(date)), Tag_year = paste(id, year, sep = "_"))



#----------------------------------------#
#### 3.3 Combine biologging summaries ####
#----------------------------------------#

## check both summary data sets quickly
glimpse(ODBASum)
glimpse(GPSSum)

## join both biologging summaries together
Inc <- full_join(ODBASum, GPSSum, by = c("device_id", "tag_year", "Date"))
  
  



#--------------------------------------#
#### 4. Create Training data Values ####
#--------------------------------------#


#----------------------------------------------------#
#### 4.1 Separate into training/unknown tag years ####
#----------------------------------------------------#

## Joining together biologging summary and tag info
## define how many tag years I have before the join
NoIncs <- length(unique(Inc$tag_year))

Inc <- Inc %>% mutate(device_id = as.character(device_id))
Inc <- inner_join(Inc, tag_info, by = "device_id") %>% drop_na(avg_odba)

## check i didn't lose any tag years during the join
stopifnot(NoIncs== length(unique(Inc$tag_year)))



## Filter out birds that were observed to have breed successfully
Breeders <- Inc %>%
            mutate(year = year(ymd(Date))) %>% 
            filter((X2017.Gos %in% c(1:8) & year == 2017) | 
                   (X2018.Gos %in% c(1:8) & year == 2018) |
                   (X2019.Gos %in% c(1:8) & year == 2019) |
                   (X2020.Gos %in% c(1:8) & year == 2020) &
                    Sex == "F" & Age == "A")


## now remove the tag_years that are in the breeders training set
## so we have a data set of tag_years with unknown breeding outcome
Females <- filter(Inc, !tag_year %in% unique(Breeders$tag_year))
Females <- filter(Females, Sex == "F")


## Function that can be used form plotting for the reminder of the script
ODBA_plot <- function(h){
  ggplot(h, aes(x = Date, y = avg_odba, col = tag_year)) + geom_point() + facet_wrap(~tag_year, scales = "free_x")
}

nsd_plot <- function(h){
  ggplot(h, aes(x = Date, y = log(nsd_), col = tag_year)) + geom_point() + facet_wrap(~tag_year, scales = "free_x")
}

ODBA_plot2 <- function(h){
  ggplot(h, aes(x = Date, y = avg_odba, col = Incubating.cert)) + geom_point() +
    facet_wrap(~tag_year, scales = "free_x") + theme(legend.position = "none")
}

nsd_plot2 <- function(h){
  ggplot(h, aes(x = Date, y = log(nsd_), col = Incubating.cert)) + geom_point() +
    facet_wrap(~tag_year, scales = "free_x") + theme(legend.position = "none")
}

ODBA_plot3 <- function(h){
  ggplot(h, aes(x = Date, y = avg_odba, col = Incubating.cert2)) + geom_point() +
    facet_wrap(~tag_year, scales = "free_x") + theme(legend.position = "none")
}

nsd_plot3 <- function(h){
  ggplot(h, aes(x = Date, y = log(nsd_), col = Incubating.cert2)) + geom_point() +
    facet_wrap(~tag_year, scales = "free_x") + theme(legend.position = "none")
}




#------------------------------------#
#### 4.2 Label training tag years ####
#------------------------------------#

## create 24 day rolling average of daily ODBA that is left centered
Breed_train <- Breeders %>% 
                group_by(tag_year) %>% 
                mutate(Roll24_ODBA = rollmean(avg_odba, k = 24, align = "left", fill = NA)) 


## Now Identify the lowest value 24 day period using a loop
## label that as incubating, the 2 day either side as buffers and everything else not incubating

## set inputs for lo0p
Tag_years <- unique(Breed_train$tag_year)
Inc_days <- 24 ## the length of minimum breeding period
Buff <- 3 ## the buffer period

## loop through each tag_year
for(j in 1:length(Tag_years)){
  
  message(j, " out of ", length(Tag_years), " Breeding periods")
  
  ## filter out one tag at a time
  Tag_sub <- filter(Breed_train, tag_year == Tag_years[j])
  
  ## find the start of the 24 day period with the lowest ODBA
  low <- which.min(Tag_sub$Roll24_ODBA)
  
  Tag_sub$status <- NA
  
  ## label incubating period
  Tag_sub$status[low:(low+Inc_days-1)] <- "Incubating"
  
  ## label buffer before
  Tag_sub$status[(low-Buff):(low-1)] <- "Buffer"
  
  ## label buffer after
  Tag_sub$status[(low+Inc_days):(low+Inc_days+Buff-1)] <- "Buffer"
  
  ## label everything else as not-incubating
  Tag_sub$status <- ifelse(is.na(Tag_sub$status) == T, "N_incubating", Tag_sub$status)
  
  if(j ==1){Breed_train2 <- Tag_sub}
  else{Breed_train2 <- rbind(Breed_train2, Tag_sub)}
  
}




#-------------------------------------------------#
#### 4.3 Calculate threshold from training set ####
#-------------------------------------------------#

## From the training data calculate quantities for ODBA, nsd and ddist of incubating and non incubating days

# for incubating days calculate 95th quantile
inc_days <- Breed_train2 %>% filter(status == "Incubating")
qinc_acc <- as.numeric(quantile(inc_days$avg_odba, probs = 0.975)) 
qinc_nsd <- as.numeric(quantile(inc_days$nsd_, probs = 0.975))
qinc_ddist <- as.numeric(quantile(inc_days$ddist, probs = 0.975))

# for non-incubating days calculate 5th quantile
no_inc_days <- Breed_train2 %>% filter(status == "N_incubating")
qnoinc_acc <- as.numeric(quantile(no_inc_days$avg_odba, probs = 0.025, na.rm = TRUE))
qnoinc_nsd <- as.numeric(quantile(no_inc_days$nsd_, probs = 0.025, na.rm = TRUE))
qnoinc_ddist <- as.numeric(quantile(no_inc_days$ddist, probs = 0.025, na.rm = TRUE))


## Save the 95% CIs from the training data
error <- qt(0.975,df=length(inc_days$avg_odba)-1)*sd(inc_days$avg_odba)/sqrt(length(inc_days$avg_odba))
mean(inc_days$avg_odba)

  


  
#-----------------------------------#
#### 5. Classify breeding events ####
#-----------------------------------#


#-----------------------------------------------------#
#### 5.1 Use training set to label all other birds ####
#-----------------------------------------------------#

##Use quantiles from nsd and ddist to label days; 2 = below qinc, 1 = in between qinc and qnoinc and 0 = above qnoinc

## label using ddist
if(qinc_ddist < qnoinc_ddist){
  Females$Incubating1 <- ifelse(Females$ddist < qinc_ddist, 2, 0)
  Females$Incubating1 <- ifelse(Females$ddist > qinc_ddist & Females$ddist < qnoinc_ddist, 1, Females$Incubating1)
}else{
  
  Females$Incubating1 <- ifelse(Females$ddist < qnoinc_ddist, 2, 0)
  Females$Incubating1 <- ifelse(Females$ddist < qinc_ddist & Females$ddist > qnoinc_ddist, 1, Females$Incubating1)
}

## End up with NAs in this column which causes problems later
Females$Incubating1 <- ifelse(is.na(Females$Incubating1) == T, 0, Females$Incubating1)



## label using nsd
if(qinc_nsd < qnoinc_nsd){
  Females$Incubating2 <- ifelse(Females$nsd_ < qinc_nsd, 2, 0)
  Females$Incubating2 <- ifelse(Females$nsd_ > qinc_nsd & Females$nsd_ < qnoinc_nsd, 1, Females$Incubating2)
}else{
  
  Females$Incubating2 <- ifelse(Females$nsd_ < qinc_nsd, 2, 0)
  Females$Incubating2 <- ifelse(Females$nsd_ < qinc_nsd & Females$nsd_ > qnoinc_nsd, 1, Females$Incubating2)
}



## label using acc
if(qinc_acc < qnoinc_acc){
  Females$Incubating.acc <- ifelse(Females$avg_odba < qinc_acc, 2, 0)
  Females$Incubating.acc <- ifelse(Females$avg_odba > qinc_acc & Females$avg_odba < qnoinc_acc, 1, Females$Incubating.acc)
}else{
  
  Females$Incubating.acc <- ifelse(Females$avg_odba < qnoinc_acc, 2, 0)
  Females$Incubating.acc <- ifelse(Females$avg_odba < qinc_acc & Females$avg_odba > qnoinc_acc, 1, Females$Incubating.acc)
}


##combine nsd, ddist and ODBA labeling by adding
Females <- Females %>% 
           mutate(Incubating.comb = (Incubating1 + Incubating2 + Incubating.acc),
                  Incubating.cert = ifelse(Incubating.comb >= 6, 1, 0))


## Data quality check ##
## find birds which had missing data for too long a period during breedin
NACheckFails <- Females %>% 
                group_by(tag_year) %>% 
                summarise(NAs = sum(ifelse(is.na(Incubating.comb)==T, 1, 0))) %>% 
                filter(NAs > 30)
NACheckFails$tag_year # lose two tag_years here

## remove the tag_years failing the check
Females <- Females %>% filter(!tag_year %in% NACheckFails$tag_year)

## Now change any NAs in the incubation columns to zero
Females <- Females %>% 
          mutate(Incubating.comb = ifelse(is.na(Incubating.comb) == T, 0, Incubating.comb),
                 Incubating.cert = ifelse(is.na(Incubating.cert) == T, 0, Incubating.cert))

## plots to check labelling at this stage if needed
#ODBA_plot2(Females)
#nsd_plot2(Females)




#---------------------------------------------------#
#### 5.2 Linearly interpolate incubation periods ####
#---------------------------------------------------#

## loop to extend incubation periods based of days were I have labeled as certainly incubating based off previous steps

## Order females by tag_year and data
Females <- Females[order(Females$tag_year, Females$Date),]

## create new column to put incubation status in. This it just extending interviews forever
Females$Incubating.cert2 <- Females$Incubating.cert

## create a list of all the tag years
tag_years <- (unique(Females$tag_year))


## This loop does the linear interpolation, I added a load of next statements that prevent it from overwriting 1s and speed it up 
for (i in 1:length(tag_years)) {
  
  #message("loop", i)
  single_tag <- filter(Females, tag_year == tag_years[i]) ## subset individual tag years for the loop
  
  for(j in 4:(nrow(single_tag)-3)) {
    
    #message("row", j)
    if(single_tag$Incubating.cert2[j] == 1){next} # need this to skip the rest of the loop
    
    for(z in 1:3) { 
      
      #message("forward", z)
      if(single_tag$Incubating.cert2[j] == 1){next} 
      ## if we are below acc qinc threshold search up to 3 days ahead for certainly incubating day
      single_tag$Incubating.cert2[j] <-  ifelse(single_tag$Incubating.comb[j] >= 2 & single_tag$Incubating.cert[j+z] == 1, 
                                                1, single_tag$Incubating.cert2[j])
    }
    
    for(z in 1:3) { 
      
      #message("backward", z)
      if(single_tag$Incubating.cert2[j] == 1){next}
      ## if we are below acc qinc threshold search up to 3 day behind for certainly incubating day
      single_tag$Incubating.cert2[j] <-  ifelse(single_tag$Incubating.comb[j] >= 2 & single_tag$Incubating.cert[j-z] == 1,
                                                1, single_tag$Incubating.cert2[j])
    }
    
    #message("missing")
    if(single_tag$Incubating.cert2[j] == 1){next}
    ## if acc runs out while incubating then use nsd and ddist to label the rest of incubation
    single_tag$Incubating.cert2[j] <-  ifelse(is.na(single_tag$Incubating.cert[j]) == TRUE & single_tag$Incubating.cert[j-1] == 1 & 
                                                single_tag$Incubating.comb[j] >= 2, 
                                              1, single_tag$Incubating.cert2[j])
  }
  
  #message("joining")
  ## re join all the individual tag years into a new data frame
  if(i == 1){all_tags <- single_tag
  } else {all_tags <- rbind(all_tags, single_tag)} 
  
}


## make sure there are no NAs in the Incubating.cert column
all_tags$Incubating.cert2 <- ifelse(is.na(all_tags$Incubating.cert2) == T, 0, all_tags$Incubating.cert2)

## plot the outputs to check the classification in the loop 
ODBA_plot3(all_tags)
#nsd_plot2(all_tags)



## Here is a verbal run through of the rules used to classify incubation
##  0. For each day in Greenland calculate average ODBA, net squared displacement (NSD) and distance between median locations of successive days (ddist)
##  1. Filter out adults females from the original data set (> 3rd breeding season as adult (17812 will be in 3rd BS in 2020))
##  2. Label known breeders based off low periods of ODBA. Label middle 24-26 days of low ODBA period as incubating and place buffers 1-2 days either side
##  3. Remove buffers from training data set sp just have the definitely incubating or not incubating days
##  4. For Incubating days calculate 95th quantile (qinc) for ODBA, nsd and ddist. For non Incubating days calculate 5th quantile (qnoinc) for ODBA, nsd and ddist.
##  5. Use quantiles from nsd and ddist to label days as; 2 = below qinc, 1 = in between qinc and qnoinc and 0 = above qnoinc
##  6. Next label days as certainly incubating if they are below the qinc threshold for ODBA, nsd and ddist
##  7. Label all days that have ODBA below max value from training data. 
##  8. Use a loop to extend the incubation periods based off days labeled as certainly incubating to those days that fall below the ODBA threshold
##  9. If acc runs out while incubating then use nsd and ddist to label the rest of incubation, must have a value of 2 in either the quantile labeling from nsd or ddist to be extended





#-----------------------------------------------#
#### 5.3 Extract Attempt start and end dates ####
#-----------------------------------------------# 

## label each incubation with it's own code ##

## Filter out the days that have been classified as incubation then label the change tag years with a 1
all_tags_inc <- all_tags %>% 
                filter(Incubating.cert2 == 1) %>% 
                mutate(day_diff = as.numeric(Date -lag(Date)),
                       change = ifelse(day_diff == 1, 0, 1))

all_tags_inc$change[1] <- 1

## use cumsum function to label each incubation
all_tags_inc$Inc_ID <- cumsum(all_tags_inc$change)

## extract all the tag years used so i don't miss the birds that didn't have any incubation
unique_tag_year <- subset(all_tags, select = c("tag_year")) %>% unique()

##now use dplyr to create summary of each putative incubation
All_attempts <- all_tags_inc %>% 
                group_by(tag_year, Inc_ID) %>% 
                summarise(attempt_start = min(Date), 
                          attempt_end = max(Date), 
                          length = as.numeric(((attempt_end-attempt_start)+1)),
                          score_tot = sum(Incubating.comb, na.rm = T),
                          max_score = max(Incubating.comb, na.rm = T),
                          ave_score = score_tot/length) %>% 
                full_join(unique_tag_year, by = "tag_year")

##change NAs in the length column to 0's
All_attempts <- All_attempts %>% 
                mutate(length = ifelse(is.na(length) == T, 0, length),
                       yday = yday(attempt_start))

## Identify any tag-years with multiple putative incubation
Dups22 <- duplicated(All_attempts$tag_year)

## extract tag_years with multiple incubation attempts
Dup_tags <- All_attempts[Dups22 == T, ]


## separate out the tags with multiple incubation and those without multiple incubation
Dup_incs <- All_attempts %>% filter(tag_year %in% c(Dup_tags$tag_year))
Non_dups <- All_attempts %>% filter(!tag_year %in% c(Dup_tags$tag_year))

## If there are no duplicates then skip the whole next section as it doesn't need doing
if(nrow(Dup_incs) == 0){Attempt_length <- Non_dups}else{
  ## Rules for picking the best putative incubation
  ## 1) Remove incubation that don't have at least one day that was classified as a 6
  ## 2) Remove incubation during the molting period if all in the molting period then just pick the one with the best average daily score 
  ## 3) Pick the incubation with the best average score
  ## (Single days with a score of 6 might mess this up)
  
  ## work out the earliest start date and maximum incubation length for each tag_year in the duplicated data set
  Dup_sum <- Dup_incs %>% 
              group_by(tag_year) %>% 
              summarise(earlist_start = min(yday, na.rm = T),
                        max_length = max(length, na.rm = T)) %>% 
              full_join(Dup_incs, by = "tag_year")
  
  ## first step remove any rows without a max score of 6
  ## this removes those 1/2 days of incubation that aren't connected to the main incubation
  ## Need to change this score to 4 or 2 in the GPS only and Acc only incubation
  Dup_sum <- filter(Dup_sum, max_score == 6)
  
  ## loop through each tag_year
  ## get unique tag_years in the data set
  unique_dups <- unique(Dup_sum$tag_year)
  
  ## run loop to select top incubation
  for(qq in 1:length(unique_dups)) {
    
    ## extract the tag year
    dup_TY <- filter(Dup_sum, tag_year == unique_dups[qq])
    
    ## if only one incubation left then just skip to the end
    if(nrow(dup_TY) == 1){dup_top <- dup_TY}
    else{
      ## take different paths if all putative incubation are in July
      ## If all in July then pick the one with the highest ave_score
      ## If ones before July then filter out for July
      if(dup_TY$earlist_start[1] <= 183){dup_TY <- filter(dup_TY, yday <= 183)}
      else{dup_TY <- filter(dup_TY, ave_score == max(dup_TY$ave_score))}
      
      ## remove 1 day incubation if there are longer ones
      if(max(dup_TY$length) == 1){dup_TY <- dup_TY}
      else{dup_TY <- filter(dup_TY, length > 1)}
      
      ## If any duplicates left then just pick the max score
      if(nrow(dup_TY) == 1){dup_top <- dup_TY}
      else{dup_top <- filter(dup_TY, ave_score == max(dup_TY$ave_score))
      dup_top <- filter(dup_top, length == max(dup_top$length))
      dup_top <- filter(dup_top, yday == min(dup_top$yday))}
    }
    
    ## join together all of the top incubation
    if(qq == 1){dup_all <- dup_top}
    else{dup_all <- rbind(dup_all, dup_top)}
    
  }
  
  ## Now join together the duplicate set and the non duplicate set
  dup_all <- select(dup_all, c("tag_year", "Inc_ID", "attempt_start", "attempt_end",   
                               "length", "score_tot", "max_score", "ave_score", "yday"))
  Attempt_length <- rbind(as.data.frame(Non_dups), as.data.frame(dup_all))
  
  
}

## Make sure there are no more duplicates
stopifnot(duplicated(Attempt_length$tag_year) == F)

## Had a check and this output is basically the same as the Incubation lengths from the methods chapter






#----------------------------------#
#### 6. LOOCV of training birds ####
#----------------------------------#


#------------------------------#
#### 6.1  Start of LOO loop ####
#------------------------------#
B=1
## Extract all the unique tag_years in the data set
Test_birds <- unique(Breeders$tag_year)

for(B in 1:length(Test_birds)){
  
  message("Training bird ", B, " out of ", length(Test_birds))

  ## remove one successful breeder at a time
  Breed_train <- Breeders %>% filter(!tag_year == Test_birds[B])
  Breed_test <- Breeders %>% filter(tag_year == Test_birds[B])
  

  #-------------------------------------------------#
  # The REST of the code is a repeat of 4.2 --> 5.3 #
  #-------------------------------------------------#
  
  ## create 24 day rolling average of daily ODBA that is left centered
  LOOBreed_train <- Breed_train %>% 
                     group_by(tag_year) %>% 
                     mutate(Roll24_ODBA = rollmean(avg_odba, k = 24, align = "left", fill = NA)) 
  
  
  ## Now Identify the lowest value 24 day period using a loop
  ## label that as incubating, the 2 day either side as buffers and everything else not incubating
  
  ## set inputs for lo0p
  Tag_years <- unique(LOOBreed_train$tag_year)
  Inc_days <- 24 ## the length of minimum breeding period
  Buff <- 3 ## the buffer period
  
  ## loop through each tag_year
  for(j in 1:length(Tag_years)){
    
    message(j, " labelling training birds")
    
    ## filter out one tag at a time
    Tag_sub <- filter(LOOBreed_train, tag_year == Tag_years[j])
    
    ## find the start of the 24 day period with the lowest ODBA
    low <- which.min(Tag_sub$Roll24_ODBA)
    
    Tag_sub$status <- NA
    
    ## label incubating period
    Tag_sub$status[low:(low+Inc_days-1)] <- "Incubating"
    
    ## label buffer before
    Tag_sub$status[(low-Buff):(low-1)] <- "Buffer"
    
    ## label buffer after
    Tag_sub$status[(low+Inc_days):(low+Inc_days+Buff-1)] <- "Buffer"
    
    ## label everything else as not-incubating
    Tag_sub$status <- ifelse(is.na(Tag_sub$status) == T, "N_incubating", Tag_sub$status)
    
    if(j ==1){LOOBreed_train2 <- Tag_sub}
    else{LOOBreed_train2 <- rbind(LOOBreed_train2, Tag_sub)}
    
  }
  

  #-------------------------------------------------#
  
  
  ## From the training data calculate quantities for ODBA, nsd and ddist of incubating and non incubating days
  
  # for incubating days calculate 95th quantile
  inc_days <- LOOBreed_train2 %>% filter(status == "Incubating")
  qinc_acc <- as.numeric(quantile(inc_days$avg_odba, probs = 0.975)) 
  qinc_nsd <- as.numeric(quantile(inc_days$nsd_, probs = 0.975))
  qinc_ddist <- as.numeric(quantile(inc_days$ddist, probs = 0.975))
  
  # for non-incubating days calculate 5th quantile
  no_inc_days <- LOOBreed_train2 %>% filter(status == "N_incubating")
  qnoinc_acc <- as.numeric(quantile(no_inc_days$avg_odba, probs = 0.025, na.rm = TRUE))
  qnoinc_nsd <- as.numeric(quantile(no_inc_days$nsd_, probs = 0.025, na.rm = TRUE))
  qnoinc_ddist <- as.numeric(quantile(no_inc_days$ddist, probs = 0.025, na.rm = TRUE))
  
  
  ## Save the 95% CIs from the training data
  error <- qt(0.975,df=length(inc_days$avg_odba)-1)*sd(inc_days$avg_odba)/sqrt(length(inc_days$avg_odba))
  mean(inc_days$avg_odba)
  

  #-----------------------------------------------------#
  

  ##Use quantiles from nsd and ddist to label days; 2 = below qinc, 1 = in between qinc and qnoinc and 0 = above qnoinc
  
  ## label using ddist
  if(qinc_ddist < qnoinc_ddist){
    Breed_test$Incubating1 <- ifelse(Breed_test$ddist < qinc_ddist, 2, 0)
    Breed_test$Incubating1 <- ifelse(Breed_test$ddist > qinc_ddist & Breed_test$ddist < qnoinc_ddist, 1, Breed_test$Incubating1)
  }else{
    
    Breed_test$Incubating1 <- ifelse(Breed_test$ddist < qnoinc_ddist, 2, 0)
    Breed_test$Incubating1 <- ifelse(Breed_test$ddist < qinc_ddist & Breed_test$ddist > qnoinc_ddist, 1, Breed_test$Incubating1)
  }
  
  ## End up with NAs in this column which causes problems later
  Breed_test$Incubating1 <- ifelse(is.na(Breed_test$Incubating1) == T, 0, Breed_test$Incubating1)
  
  
  
  ## label using nsd
  if(qinc_nsd < qnoinc_nsd){
    Breed_test$Incubating2 <- ifelse(Breed_test$nsd_ < qinc_nsd, 2, 0)
    Breed_test$Incubating2 <- ifelse(Breed_test$nsd_ > qinc_nsd & Breed_test$nsd_ < qnoinc_nsd, 1, Breed_test$Incubating2)
  }else{
    
    Breed_test$Incubating2 <- ifelse(Breed_test$nsd_ < qinc_nsd, 2, 0)
    Breed_test$Incubating2 <- ifelse(Breed_test$nsd_ < qinc_nsd & Breed_test$nsd_ > qnoinc_nsd, 1, Breed_test$Incubating2)
  }
  
  
  
  ## label using acc
  if(qinc_acc < qnoinc_acc){
    Breed_test$Incubating.acc <- ifelse(Breed_test$avg_odba < qinc_acc, 2, 0)
    Breed_test$Incubating.acc <- ifelse(Breed_test$avg_odba > qinc_acc & Breed_test$avg_odba < qnoinc_acc, 1, Breed_test$Incubating.acc)
  }else{
    
    Breed_test$Incubating.acc <- ifelse(Breed_test$avg_odba < qnoinc_acc, 2, 0)
    Breed_test$Incubating.acc <- ifelse(Breed_test$avg_odba < qinc_acc & Breed_test$avg_odba > qnoinc_acc, 1, Breed_test$Incubating.acc)
  }
  
  
  ##combine nsd, ddist and ODBA labeling by adding
  Breed_test <- Breed_test %>% 
    mutate(Incubating.comb = (Incubating1 + Incubating2 + Incubating.acc),
           Incubating.cert = ifelse(Incubating.comb >= 6, 1, 0))
  
  
  ## Data quality check ##
  ## find birds which had missing data for too long a period during breedin
  NACheckFails <- Breed_test %>% 
                  group_by(tag_year) %>% 
                  summarise(NAs = sum(ifelse(is.na(Incubating.comb)==T, 1, 0))) %>% 
                  filter(NAs > 30)
  NACheckFails$tag_year # lose two tag_years here
  
  ## remove the tag_years failing the check
  Breed_test <- Breed_test %>% filter(!tag_year %in% NACheckFails$tag_year)
  
  ## Now change any NAs in the incubation columns to zero
  Breed_test <- Breed_test %>% 
                mutate(Incubating.comb = ifelse(is.na(Incubating.comb) == T, 0, Incubating.comb),
                       Incubating.cert = ifelse(is.na(Incubating.cert) == T, 0, Incubating.cert))
  
  ## plots to check labelling at this stage if needed
  #ODBA_plot2(Breed_test)
  #nsd_plot2(Breed_test)
  

  #---------------------------------------------------#
  
  
  ## loop to extend incubation periods based of days were I have labeled as certainly incubating based off previous steps
  
  ## Order females by tag_year and data
  Breed_test <- Breed_test[order(Breed_test$tag_year, Breed_test$Date),]
  
  ## create new column to put incubation status in. This it just extending interviews forever
  Breed_test$Incubating.cert2 <- Breed_test$Incubating.cert
  
  ## create a list of all the tag years
  tag_years <- (unique(Breed_test$tag_year))
  
  
  ## This loop does the linear interpolation, I added a load of next statements that prevent it from overwriting 1s and speed it up 
  for (i in 1:length(tag_years)) {
    
    #message("loop", i)
    single_tag <- filter(Breed_test, tag_year == tag_years[i]) ## subset individual tag years for the loop
    
    for(j in 4:(nrow(single_tag)-3)) {
      
      #message("row", j)
      if(single_tag$Incubating.cert2[j] == 1){next} # need this to skip the rest of the loop
      
      for(z in 1:3) { 
        
        #message("forward", z)
        if(single_tag$Incubating.cert2[j] == 1){next} 
        ## if we are below acc qinc threshold search up to 3 days ahead for certainly incubating day
        single_tag$Incubating.cert2[j] <-  ifelse(single_tag$Incubating.comb[j] >= 2 & single_tag$Incubating.cert[j+z] == 1, 
                                                  1, single_tag$Incubating.cert2[j])
      }
      
      for(z in 1:3) { 
        
        #message("backward", z)
        if(single_tag$Incubating.cert2[j] == 1){next}
        ## if we are below acc qinc threshold search up to 3 day behind for certainly incubating day
        single_tag$Incubating.cert2[j] <-  ifelse(single_tag$Incubating.comb[j] >= 2 & single_tag$Incubating.cert[j-z] == 1,
                                                  1, single_tag$Incubating.cert2[j])
      }
      
      #message("missing")
      if(single_tag$Incubating.cert2[j] == 1){next}
      ## if acc runs out while incubating then use nsd and ddist to label the rest of incubation
      single_tag$Incubating.cert2[j] <-  ifelse(is.na(single_tag$Incubating.cert[j]) == TRUE & single_tag$Incubating.cert[j-1] == 1 & 
                                                  single_tag$Incubating.comb[j] >= 2, 
                                                1, single_tag$Incubating.cert2[j])
    }
    
    #message("joining")
    ## re join all the individual tag years into a new data frame
    if(i == 1){LOOall_tags <- single_tag
    } else {LOOall_tags <- rbind(LOOall_tags, single_tag)} 
    
  }
  
  
  ## make sure there are no NAs in the Incubating.cert column
  LOOall_tags$Incubating.cert2 <- ifelse(is.na(LOOall_tags$Incubating.cert2) == T, 0, LOOall_tags$Incubating.cert2)
  
  ## plot the outputs to check the classification in the loop 
  #ODBA_plot3(LOOall_tags)
  #nsd_plot2(LOOall_tags)
  
  
  
  ## Here is a verbal run through of the rules used to classify incubation
  ##  0. For each day in Greenland calculate average ODBA, net squared displacement (NSD) and distance between median locations of successive days (ddist)
  ##  1. Filter out adults females from the original data set (> 3rd breeding season as adult (17812 will be in 3rd BS in 2020))
  ##  2. Label known breeders based off low periods of ODBA. Label middle 24-26 days of low ODBA period as incubating and place buffers 1-2 days either side
  ##  3. Remove buffers from training data set sp just have the definitely incubating or not incubating days
  ##  4. For Incubating days calculate 95th quantile (qinc) for ODBA, nsd and ddist. For non Incubating days calculate 5th quantile (qnoinc) for ODBA, nsd and ddist.
  ##  5. Use quantiles from nsd and ddist to label days as; 2 = below qinc, 1 = in between qinc and qnoinc and 0 = above qnoinc
  ##  6. Next label days as certainly incubating if they are below the qinc threshold for ODBA, nsd and ddist
  ##  7. Label all days that have ODBA below max value from training data. 
  ##  8. Use a loop to extend the incubation periods based off days labeled as certainly incubating to those days that fall below the ODBA threshold
  ##  9. If acc runs out while incubating then use nsd and ddist to label the rest of incubation, must have a value of 2 in either the quantile labeling from nsd or ddist to be extended
  
  
  #-----------------------------------------------# 
  
  
  ## label each incubation with it's own code ##
  
  ## Filter out the days that have been classified as incubation then label the change tag years with a 1
  LOOall_tags_inc <- LOOall_tags %>% 
                      filter(Incubating.cert2 == 1) %>% 
                      mutate(day_diff = as.numeric(Date -lag(Date)),
                             change = ifelse(day_diff == 1, 0, 1))
  
  LOOall_tags_inc$change[1] <- 1
  
  ## use cumsum function to label each incubation
  LOOall_tags_inc$Inc_ID <- cumsum(LOOall_tags_inc$change)
  
  ## extract all the tag years used so i don't miss the birds that didn't have any incubation
  unique_tag_year <- subset(LOOall_tags, select = c("tag_year")) %>% unique()
  
  ##now use dplyr to create summary of each putative incubation
  LOOAll_attempts <- LOOall_tags_inc %>% 
                      group_by(tag_year, Inc_ID) %>% 
                      summarise(attempt_start = min(Date), 
                                attempt_end = max(Date), 
                                length = as.numeric(((attempt_end-attempt_start)+1)),
                                score_tot = sum(Incubating.comb, na.rm = T),
                                max_score = max(Incubating.comb, na.rm = T),
                                ave_score = score_tot/length) %>% 
                      full_join(unique_tag_year, by = "tag_year")
  
  ##change NAs in the length column to 0's
  LOOAll_attempts <- LOOAll_attempts %>% 
                      mutate(length = ifelse(is.na(length) == T, 0, length),
                             yday = yday(attempt_start))
  
  ## Identify any tag-years with multiple putative incubation
  Dups22 <- duplicated(LOOAll_attempts$tag_year)
  
  ## extract tag_years with multiple incubation attempts
  Dup_tags <- LOOAll_attempts[Dups22 == T, ]
  
  
  ## separate out the tags with multiple incubation and those without multiple incubation
  Dup_incs <- LOOAll_attempts %>% filter(tag_year %in% c(Dup_tags$tag_year))
  Non_dups <- LOOAll_attempts %>% filter(!tag_year %in% c(Dup_tags$tag_year))
  
  ## If there are no duplicates then skip the whole next section as it doesn't need doing
  if(nrow(Dup_incs) == 0){LOOAttempt_length <- Non_dups}else{
    ## Rules for picking the best putative incubation
    ## 1) Remove incubation that don't have at least one day that was classified as a 6
    ## 2) Remove incubation during the molting period if all in the molting period then just pick the one with the best average daily score 
    ## 3) Pick the incubation with the best average score
    ## (Single days with a score of 6 might mess this up)
    
    ## work out the earliest start date and maximum incubation length for each tag_year in the duplicated data set
    Dup_sum <- Dup_incs %>% 
      group_by(tag_year) %>% 
      summarise(earlist_start = min(yday, na.rm = T),
                max_length = max(length, na.rm = T)) %>% 
      full_join(Dup_incs, by = "tag_year")
    
    ## first step remove any rows without a max score of 6
    ## this removes those 1/2 days of incubation that aren't connected to the main incubation
    ## Need to change this score to 4 or 2 in the GPS only and Acc only incubation
    Dup_sum <- filter(Dup_sum, max_score == 6)
    
    ## loop through each tag_year
    ## get unique tag_years in the data set
    unique_dups <- unique(Dup_sum$tag_year)
    
    ## run loop to select top incubation
    for(qq in 1:length(unique_dups)) {
      
      ## extract the tag year
      dup_TY <- filter(Dup_sum, tag_year == unique_dups[qq])
      
      ## if only one incubation left then just skip to the end
      if(nrow(dup_TY) == 1){dup_top <- dup_TY}
      else{
        ## take different paths if all putative incubation are in July
        ## If all in July then pick the one with the highest ave_score
        ## If ones before July then filter out for July
        if(dup_TY$earlist_start[1] <= 183){dup_TY <- filter(dup_TY, yday <= 183)}
        else{dup_TY <- filter(dup_TY, ave_score == max(dup_TY$ave_score))}
        
        ## remove 1 day incubation if there are longer ones
        if(max(dup_TY$length) == 1){dup_TY <- dup_TY}
        else{dup_TY <- filter(dup_TY, length > 1)}
        
        ## If any duplicates left then just pick the max score
        if(nrow(dup_TY) == 1){dup_top <- dup_TY}
        else{dup_top <- filter(dup_TY, ave_score == max(dup_TY$ave_score))
        dup_top <- filter(dup_top, length == max(dup_top$length))
        dup_top <- filter(dup_top, yday == min(dup_top$yday))}
      }
      
      ## join together all of the top incubation
      if(qq == 1){dup_all <- dup_top}
      else{dup_all <- rbind(dup_all, dup_top)}
      
    }
    
    ## Now join together the duplicate set and the non duplicate set
    dup_all <- select(dup_all, c("tag_year", "Inc_ID", "attempt_start", "attempt_end",   
                                 "length", "score_tot", "max_score", "ave_score", "yday"))
    LOOAttempt_length <- rbind(as.data.frame(Non_dups), as.data.frame(dup_all))
    
    
  }
  
  ## Make sure there are no more duplicates
  stopifnot(duplicated(LOOAttempt_length$tag_year) == F)
  
  
  LOOAttempt_length <- LOOAttempt_length %>% ungroup() %>% as.data.frame()
  
  
  if(B == 1){LOOAll_samples <- LOOAttempt_length}else{LOOAll_samples <- rbind(LOOAll_samples, LOOAttempt_length)}
  
}
  


#----------------------------#
#### 6.2  End of LOO loop ####
#----------------------------#



#----------------------------------------------#
#### 7. Create final breeding event summary ####
#----------------------------------------------#

## Combine the classification of all unknown birds with LOO classification of training birds
All_Incubations <- rbind(Attempt_length, LOOAll_samples)
nrow(All_Incubations) 

# I have 3 tag years less, other scripts in this project have 90 tag years
# I lost two tag years because of no GPS data; "17791_2020" and "17791_2021"
# The other lost tag year is.....


## Add on the arrival dates and calculate pre-breeding length
All_Incubations2 <- Phenology %>% 
                    drop_na(Green_arrive) %>% 
                    mutate(tag_year = paste0(ID, "_", year),
                           Green_arrive = ymd(Green_arrive)) %>%
                    select(-Ice_arrive_aut) %>% 
                    right_join(All_Incubations, by = "tag_year") %>% 
                    mutate(pre_breed_length = as.numeric(ymd(attempt_start)-ymd(Green_arrive))) %>% 
                    select(!c(yday, Inc_ID, ave_score, max_score, score_tot)) %>% 
                    rename(Tag_year = tag_year, Green_dep = Green_dep_aut)

glimpse(All_Incubations2)

## write out this incubation file
write_csv(All_Incubations2, file = "Outputs/Incubation_attempt_lengths_new.csv")
