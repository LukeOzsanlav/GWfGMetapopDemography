## Luke Ozsanlav-Harris

## Create Graphs of population trends 

## packages required
library(tidyverse)
library(data.table)
library(RColorBrewer)
library(ggpubr)
library(zoo)

##Read in the data
Pop_counts <- read.csv("Data/GWfG Pop counts.csv")
Age_counts <- read.csv("Data/GWfG Islay Wexford young counts.csv")
young <- fread("Data/GWfG Islay Wexford young counts long format.csv")


##sort out pop_counts
Pop_counts$X <- NULL
Pop_counts <- filter(Pop_counts, is.na(Islay) == F)

## manipulate data frame for plotting
A <- data.frame(Pop_counts$spring, Pop_counts$Population.estimate) 
colnames(A)[2] <- "Count"; A$Population <- "Global"
B <- data.frame(Pop_counts$spring, Pop_counts$Wexford)
colnames(B)[2] <- "Count"; B$Population <- "Wexford"
C <- data.frame(Pop_counts$spring, Pop_counts$Islay)
colnames(C)[2] <- "Count"; C$Population <- "Islay"
Pop_counts2 <- rbind(A,B) %>% rbind(C)

## Plot the population count data
Pop_counts2 <- dplyr::filter(Pop_counts2, !Pop_counts.spring == 2020)
Pop_plot <- ggplot(Pop_counts2,aes(x=Pop_counts.spring, y = Count, group = Population, colour = Population)) + 
            geom_point(size = 2.5, colour = "grey") +
            geom_path(size = 1.5) +
            xlab("Year") + 
            ylab("Winter population count") + 
            theme_bw() + 
            theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), 
                  legend.position = c(0.90,0.88), axis.text=element_text(size=13), axis.title=element_text(size=17), 
                  plot.title = element_text(size=16, face="bold"), panel.grid.major.x = element_blank(), 
                  panel.grid.minor.x = element_blank(), legend.text = element_text(size=14), legend.title = element_text(size=15)) + 
            ylim(0, 36000) +  
            scale_colour_manual(values=c("#009E73", "#0072B2", "#D55E00")) +
            geom_vline(xintercept = 1999, linetype="dashed", color = "grey", size=1.25) +
            annotate(geom="text", x=1982, y=36000, label="a)",color="black", size =8)
            



## add running average smoother
roller <- function(x){rollmean(x$Perc_young, k=3, fill = NA, align = "center")}
young2 <- young %>% 
          group_by(Ringing.location) %>% 
          nest() %>% 
          mutate(run_mean = map(data, roller)) %>% 
          unnest(cols = c(data, run_mean))



Young_plot <- ggplot(young2, aes(x=data_year, group = Ringing.location, colour = Ringing.location)) + 
  geom_path(aes(y = Perc_young), size = 1, alpha = 0.3) +
  geom_line(aes(y = run_mean), size = 2) +
  xlab("Year") + 
  ylab("Percentage of juveniles") + 
  theme_bw() + 
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), 
        legend.position = "none", axis.text=element_text(size=13), axis.title=element_text(size=17), 
        plot.title = element_text(size=14, face="bold"), panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank()) + 
  ylim(0, 40) +  
  scale_colour_manual(values=c("#0072B2", "#D55E00")) +
  geom_vline(xintercept = 1999, linetype="dashed", color = "grey", size=1.25) +
  annotate(geom="text", x=1982, y=40, label="b)",color="black", size =8)


## combine the 2 plots and save as a PNG in the correct place
ggarrange(Pop_plot, Young_plot, ncol=1, nrow=2)
ggsave("Paper Plots/Figure 2- pop trends.png", 
       width = 20, height = 32, units = "cm")
