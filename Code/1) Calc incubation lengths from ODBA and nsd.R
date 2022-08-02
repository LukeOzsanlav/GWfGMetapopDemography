##Luke Ozsanlav-Harris
## combine the ODBA data from Ectone tags and quantile mapped ODBA data fro ornitela tags into one data set
## add in additional information on sex, age, ringing location then visulaise using ggplot
## finally come up with ways to class days as incubating (1) or not incubating (0)

## 1. Importing and preparing Ornitela ODBA data
## 2. Importing and preparing Ecotone ODBA data
## 3. Combing eco and orn ODBA then summarising avergae daily ODBA
## 4. Importing and organising daily movement data
## 5. Adding additional tag info
## 6. Plotting Average Daily ODBA for breeders and different regions
## 7. Prepparing labelled training data set
## 8. Classifying incubations
## 9. Extract Attempt start and end dates

## Created: 20/05/2020
## Updated: 05/05/2020. Updates was to include nsd in the classification process

##loading Packages
library(tidyverse)
library(lubridate)
library(data.table)




########
## 1. ##
######## Importing and preparing Ornitela ODBA data

##reading in and binding ODBA values for ll bursts
setwd("~/Acc/ALL_ORNI_DATA/ODBA per burst/Quantile mapped ODBA")
ODBA_orn <- fread("All_Orn_ODBA_quantile_map.csv")

##selecting required columns  
ODBA_orn <- subset(ODBA_orn, select = c("UTC_datetime", "device_id", "ODBA_quant", "burst_no_ID",
                                        "datatype"))
ODBA_orn %>% setnames(old=c("device_id", "ODBA_quant", "datatype"), 
                      new=c("Tag_ID", "ODBA", "Data_Type"))

##filtering out rows with only NAs in
ODBA_orn <- ODBA_orn %>% select_if(~sum(!is.na(.)) > 0)

##setting timestamp object
ODBA_orn$UTC_datetime <- ymd_hms(ODBA_orn$UTC_datetime)
  
##adding year month and hour columns
ODBA_orn <- ODBA_orn %>% mutate(year =  year(UTC_datetime),
                                        month = month(UTC_datetime), 
                                        hour = hour(UTC_datetime),
                                        tag_date = paste(Tag_ID, date(UTC_datetime), sep = "_"),
                                        tag_year =paste(Tag_ID, year(UTC_datetime), sep = "_"),
                                        date= date(UTC_datetime))  

####Now can check for duplicated observations####
check4 <- subset(ODBA_orn, select = c("UTC_datetime", "Tag_ID", "ODBA")) %>% duplicated
sum(check4) ##how many rows are duplicates
ODBA_orn <-  ODBA_orn %>% filter(!check4)





########
## 2. ##
######## Importing and preparing Ecotone ODBA data

##reading in and binding ODBA values for ll bursts
setwd("~/Ecotone Data/ODBA_per_burst_ecotone")
ODBA_eco <- fread("Eco_ODBA_per_burst.csv")

##selecting required columns  
ODBA_eco <- subset(ODBA_eco, select = c("UTC_datetime", "Tag_ID", "ODBA", "burst_no_ID", "Data_Type"))

##filtering out rows with only NAs in
ODBA_eco <- ODBA_eco %>% select_if(~sum(!is.na(.)) > 0)

##setting timestamp object
ODBA_eco$UTC_datetime <- ymd_hms(ODBA_eco$UTC_datetime)

##adding year month and hour columns
ODBA_eco <- ODBA_eco %>% mutate(year =  year(UTC_datetime),
                                        month = month(UTC_datetime), 
                                        hour = hour(UTC_datetime),
                                        tag_date = paste(Tag_ID, date(UTC_datetime), sep = "_"),
                                        tag_year =paste(Tag_ID, year(UTC_datetime), sep = "_"),
                                        date= date(UTC_datetime))  

####Now can check for duplicated observations####
check4 <- subset(ODBA_eco, select = c("UTC_datetime", "Tag_ID", "ODBA")) %>% duplicated
sum(check4) ##how many rows are duplicates
ODBA_eco <-  ODBA_eco %>% filter(!check4)





########
## 3. ##
######## Combing eco and orn ODBA then summarizing avergae daily ODBA

##check two data sets have same column names
colnames(ODBA_eco)
colnames(ODBA_orn)

##bind the two together
ODBA_data <- rbind(ODBA_eco, ODBA_orn)

## write out the combined ODBA data set if required
#write_csv(ODBA_data, path = "~/Incubation Lengths/All ODBA data/ODBA_data_combined.csv")

##average daily ODBA for each tag then put back in original data set so retain all additional columns
ODBA_wAvg <- ODBA_data %>% group_by(tag_date) %>% summarise(
  avg_odba = mean(ODBA)) %>% full_join(ODBA_data, by = "tag_date") 
  
##removing lines so just have one value per tag day
ODBA_Avg_pday <- distinct(ODBA_wAvg, tag_date, .keep_all= TRUE)


#### trim ODBA data to just Greenland period ####
##read in phenology data
setwd("~/Migration/Phenology")
Phenology <- fread("Full_phenology_fixes.csv")
Phenology$V1 <- NULL

##set a tag_year column
Phenology$tag_year <- paste(Phenology$ID, Phenology$year, sep = "_")

##set phenology columns to lubridate timestamps and give 1 day buffer before departure and after arrival date
Phenology$Green_dep <- (ymd_hms(paste0(Phenology$Green_dep, " 00:00:00"))) - 86400
Phenology$Green_arrive <- (ymd_hms(paste0(Phenology$Green_arrive, " 00:00:00"))) + 86400

##now join the two data sets
ODBA_Avg_pday <- inner_join(ODBA_Avg_pday, Phenology, by = "tag_year")
unique(ODBA_Avg_pday$tag_year) ## check haven't lost lots of tag_years

## filter out the period in Greenland
ODBA_Avg_pday_Gr <- ODBA_Avg_pday %>%  filter(UTC_datetime > Green_arrive & UTC_datetime < Green_dep)


#### filtering out tags that don't have data for the full breeding season or were dead ####
## create data set with tag years with little data removed
no_of_days <- as.data.frame(table(ODBA_Avg_pday_Gr$tag_year)) %>% filter(Freq > 30)
colnames(no_of_days)[1] <- "tag_year"

##joining the above tag list with ODBA data so only have tag years that have data for most of a breeding sesaon
ODBA_Avg_pday_Gr <- inner_join(ODBA_Avg_pday_Gr, no_of_days, by = "tag_year")
ODBA_Avg_pday_Gr$Freq <- NULL
length(unique(ODBA_Avg_pday_Gr$tag_year)) ##should be the length as no_of_days data set




########
## 4. ##
######## Importing and organising daily movement data

##reading in movement data
setwd("~/Incubation Lengths/Daily movement characteristics")
Eco_move <- fread("Eco_movement_metrics.csv")
Orn_move <- fread("Orni_movement_metrics.csv")

##bind Eco and Orni data together
move <- rbind(Eco_move, Orn_move)
move$V1 <- NULL

##now bind to the ODBA data
##select the columns needed from ODBA data and rename so they match the movement data
ODBA_Avg_pday_Gr <- subset(ODBA_Avg_pday_Gr, select = c("tag_date","avg_odba", "Tag_ID", "year.x", "tag_year", "date"))
ODBA_Avg_pday_Gr %>% setnames(old=c("Tag_ID", "year.x", "tag_year"), 
                      new=c("id", "year", "Tag_year"))
ODBA_Avg_pday_Gr$date <- as.character(ODBA_Avg_pday_Gr$date)

##create tag_date column in move data
move$tag_date <- paste(move$id, move$date, sep = "_")
move$Tag_year <- paste(move$id, move$year, sep = "_")

##now full join the two together
Inc <- full_join(move, ODBA_Avg_pday_Gr)

##filter out the tags that had accelerometers; remove GPS only tags
acc <- Inc %>% group_by(Tag_year) %>% summarise(
  avg = mean(avg_odba, na.rm= TRUE)) %>% filter(avg > 0)
Inca <- inner_join(Inc, acc[,1], by = "Tag_year")







########
## 5. ##
######## Adding additional tag info

##importing dataset with breeding histories and adding info Inca data
setwd("~/Additional data files")
tag_info <- read.csv("Tagged bird summary data.csv")

##creating ID column to bind everthing by
tag_info$Bird.ID <- as.character(tag_info$Bird.ID)
tag_info$S.N <- as.character(tag_info$S.N)
tag_info$S.N <- ifelse(is.na(tag_info$S.N) == TRUE, tag_info$Bird.ID, tag_info$S.N)

##selecting important columns
tag_info2 <- subset(tag_info, select = c("Age", "Mass", "Sex", "Ringing.location", 
                                         "S.N", "X2017.Gos","X2018.Gos", "X2019.Gos", "X2020.Gos"))
colnames(tag_info2)[5] <- "id"

###joining together odba data and tag info
Incex <- inner_join(Inca, tag_info2, by = "id")
length(unique(Inca$Tag_year))
length(unique(Incex$Tag_year)) ## two lengths should be the same or have lost tag_years





########
## 6. ##
######## Plotting Average Daily ODBA and nsd for breeders and different locations

##filering out the required data for each plot
Breeders_2017 <- Incex %>% filter(year == 2017 & X2017.Gos %in% c(1:8) & Sex == "F")

Breeders_2018 <- Incex %>% filter(year == 2018 & X2018.Gos %in% c(1:8) & Sex == "F")

Breeders_2019 <- Incex %>% filter(year == 2019 & X2019.Gos %in% c(1:8) & Sex == "F")

Breeders_2020 <- Incex %>% filter(year == 2020 & X2020.Gos %in% c(1:8) & Sex == "F")

Breeders <- rbind(Breeders_2020, Breeders_2019, Breeders_2018, Breeders_2017)

Irish <- Incex %>% filter(Ringing.location %in% c("LIRO", "SESK", "WEXF", "HVAN", "LOSW") & Sex == "F")

Scot <- Incex %>% filter(!Ringing.location %in% c("LIRO", "SESK", "WEXF", "HVAN", "LOSW") & Sex == "F")

Males <- Incex %>% filter(Sex == "M")

Females <- Incex %>% filter(Sex == "F" & Age == "A") ## needs adapting after next season so keep birds ringed as juvs but not in 2nd/3rd breeding season

##plotting ODBA through breeding season
##functions for plotting ODBA
ODBA_plot <- function(h){
  ggplot(h, aes(x = date, y = avg_odba, col = Tag_year)) + geom_point() +
    facet_wrap(~Tag_year, scales = "free_x")
}

nsd_plot <- function(h){
  ggplot(h, aes(x = date, y = log(nsd_), col = Tag_year)) + geom_point() +
    facet_wrap(~Tag_year, scales = "free_x")
}
  
###Plotting the various groups
ODBA_plot(Breeders)
#ODBA_plot(Irish)
#ODBA_plot(Scot)
#ODBA_plot(Males)

nsd_plot(Breeders)
#nsd_plot(Irish)
#nsd_plot(Scot)
#nsd_plot(Males)
#nsd_plot(Females)




########
## 7. ##
######## Preparing labeled training data set

##reading in data set I created of labeled incubation days for known breeders
setwd("~/Incubation Lengths")
labelled_breeders <- fread("Confrimed breeders labelled ODBA.csv")

##creating column for Not incubating
labelled_breeders$N_incubating <- ifelse(labelled_breeders$Incubating == 0 & labelled_breeders$Buffer == 0, 1, 0)

##join the labeled breeders and Phenology data sets
labelled_breeders <- inner_join(labelled_breeders, Phenology, by = "tag_year")

## filter out the period in Greenland
labelled_breeders$UTC_datetime <- as.POSIXct(labelled_breeders$UTC_datetime, format = "%d/%m/%Y %H:%M")
labelled_breeders <- labelled_breeders %>%  filter(UTC_datetime > Green_arrive & UTC_datetime < Green_dep)

## join the labeled days with the ODBA, nsd and ddist data
labbe <- subset(labelled_breeders, select = c("Incubating",  "Buffer", "tag_date", "N_incubating"))
Label <- full_join(Breeders, labbe)

##Now sort out the rest of the incubating columns (these were dates near the end of the season that i did not label)
Label$Incubating <- ifelse(is.na(Label$Incubating) == TRUE, 0, Label$Incubating)
Label$Buffer <- ifelse(is.na(Label$Buffer) == TRUE, 0, Label$Buffer)
Label$N_incubating <- ifelse(is.na(Label$N_incubating) == TRUE, 1, Label$N_incubating)

##removed: 21/05/2020; went with a threshold method instead of logistic regression
##running logistic regression of Incuabating~ avg_ODBA
#labelled_breeders_logistic <- labelled_breeders %>% filter(Buffer == 0) ##filter out buffer values
#logistic_model <- glm(Incubating ~ avg_odba, family=binomial(link='logit'), data=labelled_breeders) ##run the model
#summary(logistic_model)
##plot the logisitc regression model
#ggplot(labelled_breeders, aes(x=avg_odba, y=Incubating)) + geom_point() + 
  #stat_smooth(method="glm", method.args=list(family="binomial"), se=FALSE)





########
## 8. ##
######## Classifying incubations

##From the traning data calculate quantiles for ODBA, nsd and ddist of incubating and non incubating days
# for incubating days calculate 95th quanile
incubating_days <- Label %>% filter(Incubating == 1)
qinc_acc <- as.numeric(quantile(incubating_days$avg_odba, probs = 0.95)) 
qinc_nsd <- as.numeric(quantile(incubating_days$nsd_, probs = 0.95))
qinc_adist <- as.numeric(quantile(incubating_days$ave_dist, probs = 0.95))
qinc_ddist <- as.numeric(quantile(incubating_days$ddist, probs = 0.95))

# for non-incubating days calculate 5th quanile
non_incubating_days <- Label %>% filter(Incubating == 0)
qnoinc_acc <- as.numeric(quantile(non_incubating_days$avg_odba, probs = 0.05, na.rm = TRUE))
qnoinc_nsd <- as.numeric(quantile(non_incubating_days$nsd_, probs = 0.05, na.rm = TRUE))
qnoinc_adist <- as.numeric(quantile(non_incubating_days$ave_dist, probs = 0.05, na.rm = TRUE))
qnoinc_ddist <- as.numeric(quantile(non_incubating_days$ddist, probs = 0.05, na.rm = TRUE))

##Use quanitles from nsd and ddist to label days; 2 = below qinc, 1 = in between qinc and qnoinc and 0 = above qnoinc
##label using ddist
Females$Incubating1 <- ifelse(Females$ddist < qinc_ddist, 2, 0)
Females$Incubating1 <- ifelse(Females$ddist > qinc_ddist & Females$ddist < qnoinc_ddist, 1, Females$Incubating1)

##lable using nsd
Females$Incubating2 <- ifelse(Females$nsd_ < qinc_nsd, 2, 0)
Females$Incubating2 <- ifelse(Females$nsd_ > qinc_nsd & Females$nsd_ < qnoinc_nsd, 1, Females$Incubating2)

##combine nsd and ddist labelling by adding
Females$Incubating.comb <- Females$Incubating1 + Females$Incubating2

##label as definately incubating if is labelled with 2 for nsd and ddist and below ODBA threshold
Females$Incubating.cert <- ifelse(Females$Incubating.comb == 4 & Females$avg_odba < qinc_acc, 1, 0)

##now label points below acc threshold (qinc) to calculate the full incubation length
Females$Incubating.acc <- ifelse(Females$avg_odba < max(incubating_days$avg_odba), 1, 0)


## loop to extend incubation periods based of days were I have labbelled as certainly incubating based off previoous steps
##create a list of all the tag years
tag_years <- (unique(Females$Tag_year))


for (i in 1:length(tag_years)) {
  
  single_tag <- filter(Females, Tag_year == tag_years[i]) ## subset individual tag years for the loop
  
    for(j in 2:nrow(single_tag)) {
        
      for(z in 1:3) { 
            
            ## if we are below acc qinc threshold search up to 3 days ahead for certainly incubating day
            single_tag$Incubating.cert[j] <-  ifelse(single_tag$Incubating.acc[j] == 1 & single_tag$Incubating.cert[j+z] == 1, 
                                                 1, single_tag$Incubating.cert[j])
      }
        ## if we are below acc qinc threshold search 1 day behind for certainly incubating day
        single_tag$Incubating.cert[j] <-  ifelse(single_tag$Incubating.acc[j] == 1 & single_tag$Incubating.cert[j-1] == 1,
                                                 1, single_tag$Incubating.cert[j])
        ## if acc runs out while incubaitng then use nsd and ddist to lable the rest of incubation
        single_tag$Incubating.cert[j] <-  ifelse(is.na(single_tag$Incubating.cert[j]) == TRUE & single_tag$Incubating.cert[j-1] == 1 & 
                                                   single_tag$Incubating.comb[j] >= 2, 
                                                 1, single_tag$Incubating.cert[j])
   }
  ## re join all the individual tag years into a new data frame
  if(i == 1){all_tags <- single_tag
  } else {all_tags <- rbind(all_tags, single_tag)} 
  
}

## plot the outputs to check loop hasn't gone awol
ODBA_plot(all_tags)
nsd_plot(all_tags)


## Here is a verbal run through of the rules used to classify incubation
##  0. For each day in Greenland calculate average ODBA, net squared displacement (NSD) and distance between median locations of successive days (ddist)
##  1. Filter out adults females from the original data set (> 3rd breeding season as adult (17812 will be in 3rd BS in 2020))
##  2. Label known breeders based off low periods of ODBA. Label middle 24-26 days of low ODBA period as incubating and place buffers 1-2 days either side
##  3. Remove buffers from training data set sp just have the definitely incubating or not incubating days
##  4. For Incubating days calculate 95th quantile (qinc) for ODBA, nsd and ddist. For non Incubating days calculate 5th quantile (qnoinc) for ODBA, nsd and ddist.
##  5. Use quanitles from nsd and ddist to label days as; 2 = below qinc, 1 = in between qinc and qnoinc and 0 = above qnoinc
##  6. Next label days as certainly incubating if they are below the qinc threshold for ODBA, nsd and ddist
##  7. Label all days that have ODBA below max value from training data. 
##  8. Use a loop to extend the incubation periods based off days labeled as certainly incubating to those days that fall below the ODBA threshold
##  9. If acc runs out while incubating then use nsd and ddist to label the rest of incubation, must have a value of 2 in either the quantile labeling from nsd or ddist to be extended




########
## 9. ##
######## Extract Attempt start and end dates

## Changed: 04/06/2020
## Changed method to include movement metrics so this bit of the code changed as well

##Set date as a date object
all_tags$date <- as.Date(all_tags$date)

##create data frame with attmpt lengths
unique_tag_year <- subset(all_tags, select = c("Tag_year")) %>% unique()

##now use dplyr to create summary
Attempt_length <- all_tags %>% 
                  filter(Incubating.cert == 1) %>% 
                  group_by(Tag_year) %>% 
                  summarise(attempt_start = min(date), attempt_end = max(date), 
                            length =((attempt_end-attempt_start)+1), tag_year = min(Tag_year)) %>% 
                  full_join(unique_tag_year)

##change NAs in the length column to 0's
Attempt_length$length <- ifelse(is.na(Attempt_length$length) == TRUE, 0, Attempt_length$length)

##historgam of attmept length
hist(Attempt_length$length, breaks = 8)


##calculate length between breeding attempt start and W Greenland arrival
##join with Phenology data
Phenology %>% setnames(old=c("tag_year"), new=c("Tag_year"))
Attempt_length <- inner_join(Attempt_length, Phenology, by = "Tag_year")

##Set dates as date object
Attempt_length$attempt_start <- as.Date(Attempt_length$attempt_start)
Attempt_length$Green_arrive <- as.Date(Attempt_length$Green_arrive)

##calculate pre breeding length
Attempt_length$pre_breed_length <- ifelse(is.na(Attempt_length$attempt_start) == FALSE, 
                                    Attempt_length$attempt_start - Attempt_length$Green_arrive, NA)

##boxplot of pre breeding length
Attempt_length$year <- as.factor(Attempt_length$year)
ggplot(Attempt_length, aes(x= year, y= pre_breed_length, fill= year)) + geom_boxplot()

##write out attempts lengths file
Attempt_length <- subset(Attempt_length, select = c("Tag_year", "attempt_start", "attempt_end", "length", "Green_dep", "UK_dep",
                                                    "Ice_arrive_spring", "year", "ID", "Green_arrive", "pre_breed_length", "Ice_dep_spring"))
Attempt_length <- filter(Attempt_length, !Tag_year == "17791_2020") # remove this bird as the GPS failed and took no fixes before or during incubation
setwd("~/Incubation Lengths/Incubation attempts output")
write.csv(Attempt_length, file = "Incubation_attempt_lengths_new.csv", row.names = F)

