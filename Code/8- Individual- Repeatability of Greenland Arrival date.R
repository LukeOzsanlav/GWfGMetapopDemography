## Luke Ozsanlav-Harris

## Are Greenland Arrivla dates repeatable??

## Luke Ozsanlav-Harris

## Repeatability of phenolgical dates within individuals
## Specifically to test whether birds that arrive early into Greenland for that year are consistently early arrivers

## Created: 08/01/2020

## packages required
pacman::p_load(tidyverse, data.table, rptR) # link to rptR vignette: https://cran.r-project.org/web/packages/rptR/vignettes/rptR.html




#----------------------------#
#### 1. Read in data sets ####
#----------------------------#

## read in the phenology data sets
Phenol <- fread("Outputs/Full_phenology_upto2021.csv")

## read in tag info
Info <- fread("Data/Metadata/Tagged bird summary data new.csv")



#-----------------------------------#
#### 2. Format and join datasets ####
#-----------------------------------#


## Format Phenolgy data set
Phen <- Phenol %>% 
        group_by(year) %>% 
        mutate(Green_centre = yday(Green_arrive)-mean(yday(Green_arrive), na.rm = T)) %>% 
        ungroup() %>% 
        drop_na(Green_centre) %>% 
        select(tag_year, year, ID, Green_centre)

## Format tag info
tag_info <- Info %>% 
            mutate(across(c(Bird.ID, S.N), as.character)) %>% 
            mutate(S.N = ifelse(is.na(S.N) == TRUE, Bird.ID, S.N)) %>% 
            rename(ID = S.N) %>% 
            select(Age, Sex, Ringing.location, ID, Accelerometer.)

## join the two data sets together
Rep <- left_join(Phen, tag_info, by = "ID")


## Remove birds that are not adult females with accelerometer tags
Rep <- Rep %>% filter(Age == "A" & Sex == "F" & Accelerometer. == "Y")



#-----------------------------------#
#### 3. Summarise data collected ####
#-----------------------------------#

## How many breeding seaosns worth of data did I collect
DataSum <- as.data.frame(table(Rep$ID))
table(DataSum$Freq)


#------------------------------------#
#### 4. Format data for modeling ####
#------------------------------------#

## remove tags that dpo not have repeats 
Dups <- Rep %>% select(ID) %>% duplicated() %>% filter(Rep, .) %>% select(ID) %>% unique()
Reps <- Rep %>% filter(ID %in% Dups$ID)



#-------------------------------------#
#### 5. Run repeatability analysis ####
#-------------------------------------#

## NOTE: I tried running this analysis seperately for Scotland/Irish birds but the Irish model 
##       had lots of problems with singular fits so have to run them together

## use the rptRpackage to test the repeatability of Greenland arrival dates
rep1 <- rpt(Green_centre ~ (1 | ID), grname = "ID", data = Reps, 
            datatype = "Gaussian", nboot = 1000, npermut = 0)
summary(rep1)
plot(rep1)

## Model Output
# Repeatability estimation using the lmm method
# 
# Call = rpt(formula = Green_centre ~ (1 | ID), grname = "ID", data = Reps, datatype = "Gaussian", nboot = 1000, npermut = 0)
# 
# Data: 101 observations
# ----------------------------------------
#         
#         ID (38 groups)
# 
# Repeatability estimation overview: 
#    R     SE   2.5%  97.5% P_permut  LRT_P
# 0.44  0.107  0.212  0.627       NA      0
# 
# Bootstrapping and Permutation test: 
#             N   Mean Median   2.5%  97.5%
# boot     1000  0.431  0.437  0.212  0.627
# permut      1     NA     NA     NA     NA
# 
# Likelihood ratio test: 
# logLik full model = -221.522
# logLik red. model = -228.0073
# D  = 13, df = 1, P = 0.000158
# 
# ----------------------------------------
