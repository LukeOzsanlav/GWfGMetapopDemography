## LukOzsanlav-Harris

## Created 27/03/2023

## Compare morphologicla measurement of two sub-populations

## Load packages required
pacman::p_load(tidyverse, data.table, smatr, performance, ggpubr)



#---------------------------------------------#
#### 1. Read in morphology data and format ####
#---------------------------------------------#

## importing dataset with additonal information on individuals 
Morph <- read.csv("Data/Metadata/Tagged bird summary data new.csv")

## Extract the birds used in this analysis
Morph <- filter(Morph, is.na(S.N) == F & Age == "A") %>% rename(ID = S.N)

## selecting important columns 
Mor <- subset(Morph, select = c("ID", "Sex", "Mass", "Wing", "H.B", "Ringing.location"))

## Add sub-population column
Mor <- mutate(Mor, sub_pop = ifelse(Ringing.location == "WEXF" | Ringing.location == "SESK" | Ringing.location == "HVAN", "Wexford", "Islay"))

## drop rows if we don't have the morpho data, only 1 bird that did not have weight
Mor <- Mor %>% drop_na(Mass, Wing)

## Make sure certain columns are  numeric
Mor <- Mor %>% mutate(across(c(Mass, Wing, H.B), as.character)) %>%  mutate(across(c(Mass, Wing, H.B), as.numeric))



#--------------------------------------#
#### 2. Calculate scaled mass index ####
#--------------------------------------#

## This 2013 paper suggests using head and bill as the size metric for scaled mass index as this most closely correlated with fat stores
## https://onlinelibrary.wiley.com/doi/full/10.1111/jofo.12019

## I am using the scaled mass index in Peig and Green 2009
## https://onlinelibrary.wiley.com/doi/full/10.1111/j.1600-0706.2009.17643.x


## plot the data first
plot(Mor$Mass, Mor$H.B) # one outlier, probable transcription error as it down as this in the RAW data sheet
plot(Mor$Mass, Mor$Wing)

## Remove the outlier from head and bill
Mor <- filter(Mor, H.B < 120)

## standardised major axis (SMA) regression on ln-transformed data
Mod <- sma(Mass ~ H.B, data = Mor, log = 'xy', method= "SMA")
summary(Mod) # slope estimate is 3.212287
plot(Mod)

## calcualte scaled mass index
slope <- as.numeric(coef(Mod)[2]) # first get the slope estimate
Mor <- Mor %>% mutate(SMI = Mass*((mean(H.B)/H.B)^slope) )

## plot of the data
ggplot() + geom_point(data = Mor, aes(x = H.B, y = Mass, colour = SMI)) + theme_classic()



#------------------------------------------#
#### 3. Model comparing sub-populations ####
#------------------------------------------#

## set sub population as a factor
Mor$sub_pop <- as.factor(Mor$sub_pop)

## Model SMI vs sub_pop
SMI_Mod <- lm(SMI ~ sub_pop, data = Mor)
summary(SMI_Mod)
confint(SMI_Mod)
check_model(SMI_Mod)

## Model Mass vs sub_pop
Mass_Mod <- lm(Mass ~ sub_pop, data = Mor)
summary(Mass_Mod)
confint(Mass_Mod)
check_model(Mass_Mod)

## Model Wing vs sub_pop
Wing_Mod <- lm(Wing ~ sub_pop, data = Mor)
summary(Wing_Mod)
confint(Wing_Mod)
check_model(Wing_Mod)




#-----------------------------------#
#### 4. Box plots of comparisons ####
#-----------------------------------#


## create the box plot of SMI
SMIbp <- ggplot(data = Mor, aes(y = SMI, x = sub_pop, fill = sub_pop)) + 
  geom_jitter(stroke = 1, alpha = 0.5, colour = "darkgrey", width = 0.1) +
  geom_boxplot(alpha = 0.5, width = 0.4, outlier.shape = NA) +
  theme_light() +
  #scale_y_continuous(breaks = seq(-140, -20, by = 20)) +
  scale_fill_manual(values=c("#0072B2", "#D55E00")) +
  ylab("Scaled Mass Index") + xlab("Sub-population") + labs(fill = "Sub-population") +
  theme(axis.text=element_text(size=14), panel.grid.minor = element_blank(), legend.title = element_text(size=18),
        axis.text.x = element_text(size=14), axis.title=element_text(size=18), panel.grid.major.x = element_blank(),
        legend.text = element_text(size=14),  axis.title.x = element_blank(), legend.spacing.x = unit(0.8, 'cm'))

## create the box plot of SMI
Wingbp <- ggplot(data = Mor, aes(y = Wing, x = sub_pop, fill = sub_pop)) + 
  geom_jitter(stroke = 1, alpha = 0.5, colour = "darkgrey", width = 0.1) +
  geom_boxplot(alpha = 0.5, width = 0.4, outlier.shape = NA) +
  theme_light() +
  #scale_y_continuous(breaks = seq(-140, -20, by = 20)) +
  scale_fill_manual(values=c("#0072B2", "#D55E00")) +
  ylab("Wing length (mm)") + xlab("Sub-population") + labs(fill = "Sub-population") +
  theme(axis.text=element_text(size=14), panel.grid.minor = element_blank(), legend.title = element_text(size=18),
        axis.text.x = element_text(size=14), axis.title=element_text(size=18), panel.grid.major.x = element_blank(),
        legend.text = element_text(size=14),  axis.title.x = element_blank(), legend.spacing.x = unit(0.8, 'cm'))

## create the box plot of SMI
Massbp <- ggplot(data = Mor, aes(y = Mass, x = sub_pop, fill = sub_pop)) + 
  geom_jitter(stroke = 1, alpha = 0.5, colour = "darkgrey", width = 0.1) +
  geom_boxplot(alpha = 0.5, width = 0.4, outlier.shape = NA) +
  theme_light() +
  #scale_y_continuous(breaks = seq(-140, -20, by = 20)) +
  scale_fill_manual(values=c("#0072B2", "#D55E00")) +
  ylab("Mass (g)") + xlab("Sub-population") + labs(fill = "Sub-population") +
  theme(axis.text=element_text(size=14), panel.grid.minor = element_blank(), legend.title = element_text(size=18),
        axis.text.x = element_text(size=14), axis.title=element_text(size=18), panel.grid.major.x = element_blank(),
        legend.text = element_text(size=14), axis.title.x = element_blank(), legend.spacing.x = unit(0.8, 'cm'))

## Combine the plot
Combo <- ggarrange(SMIbp, Wingbp, Massbp, common.legend = T, legend = "bottom", nrow= 1)


## save the plot
ggsave(plot = Combo, filename =  "Paper Plots/Supp Fig 3- Morphology comparison.png",
       width = 30, height = 22, units = "cm")


