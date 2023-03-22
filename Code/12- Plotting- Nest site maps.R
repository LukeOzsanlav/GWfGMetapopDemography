## Luke Ozsanlav-Harris

## Create a plot that shows the nest sites for all individuals colour by sub-popualtion
## Then create a plot for individuals with nest sites in mutiple years

## Load packages
pacman::p_load(tidyverse, data.table, lubridate, move, rnaturalearth, sf, ggspatial, MetBrewer)

## Workflow
# Read in tracking data
# Crop to the nesting period
# Read in basemsap of Greenland
# Plot map of all nest by individual
# Plot maps of nest for reapeat nester


#---------------------------------------------#
#### 1. Read in Incubation data and format ####
#---------------------------------------------#

## Read in the incubation data set
Inc <- fread("Outputs/Incubation_attempts_for_models.csv")

## set date format
Inc <- Inc %>% mutate(across(c(attempt_end, Green_arrive, attempt_start, Ice_arrive_spring), ymd)) %>% 
       dplyr::select(c(attempt_end, Green_arrive, attempt_start, Ice_arrive_spring, Tag_year, Ringing.location))

##filter the birds that were known to incubate
Inc_br <- Inc %>% filter(is.na(attempt_start) == F)

## set dates to as.date format
Inc_br$attempt_start <- ymd_hms(paste0(Inc_br$attempt_start, "_", "12:00:00"))
Inc_br$attempt_end <- ymd_hms(paste0(Inc_br$attempt_end, "_", "12:00:00"))




#--------------------------------#
#### 2. Read in tracking data ####
#--------------------------------#

## load in GPS tracking data 
GWF.df <- readRDS("Tracking data/Eco_Orn_All_MBformat.RDS")

## use luibridate to parse timestamp
GWF.df$UTC_datetime <- ymd_hms(GWF.df$UTC_datetime)

## create tag_year column
GWF.df$Tag_year <- paste0(GWF.df$device_id, "_", year(GWF.df$UTC_datetime))

## drop the Ecotone tags
GWF.df <- GWF.df %>% 
          mutate(tagtype = substr(device_id, 1, 1)) %>% 
          filter(tagtype == 1)




#-----------------------------------#
#### 3. Calculate the nest sites ####
#-----------------------------------#

#### filter the GPS data to extract periods when birds were known to be incubating
## Join the two data sets
GWF.br <- inner_join(GWF.df, Inc_br, by = "Tag_year")


## filter the incubation windows
GWF.inc <- GWF.br %>% filter(UTC_datetime > attempt_start & UTC_datetime < attempt_end)

## calculate median lat and median long for each tag year
length(unique(GWF.inc$Tag_year))
med_nest_locs <- GWF.inc %>% 
  group_by(Tag_year, Ringing.location, device_id) %>% 
  summarise(breeding_lat = median(Latitude), 
            breeding_long = median(Longitude)) %>% 
  mutate(sub_pop = ifelse(Ringing.location == "WEXF" | Ringing.location == "SESK", "Wexford", "Islay"))

## convert these locaitons into an sf object
NestLocs <- st_as_sf(med_nest_locs, coords = c("breeding_long", "breeding_lat"), 
                    crs = 4326)

## map extents
miny <- min(med_nest_locs$breeding_lat) - 1.25
maxy <-max(med_nest_locs$breeding_lat) + 1.25
minx <-min(med_nest_locs$breeding_long) - 3
maxx <-max(med_nest_locs$breeding_long) + 3




#---------------------------------------#
#### 5. Create base map for plotting ####
#---------------------------------------#

## download global data sets for base maps
countries <- ne_coastline(scale = "medium", returnclass = "sf")
countries2 <- ne_countries(scale = "medium", returnclass = "sf")

## ## Define parameters for reading out plots
## Define device to read plots out as e.g. tiff/jpeg
device <- "png"

## define units for plot size - usually mm
units <- "mm"

## define plot resolution in dpi - 300 usually minimum
dpi <- 500




#--------------------------------#
#### 6. Create nest site maps ####
#--------------------------------#

## Create map of all nesting locations
P1 <- ggplot() + 

  # add the filled in countries
  geom_sf(data = countries2, aes(geometry = geometry), fill = "#BDC3C7", colour = "#BDC3C7", alpha = 0.6) +
  
  # add the nest locations
  geom_sf(data = NestLocs, aes(geometry = geometry, colour = sub_pop),  size = 4, alpha = 0.7) +
  scale_colour_manual(name = "Sub-population",
                      values =c("Islay" = "#0072B2", "Wexford" = "#D55E00"),
                      labels = c("Islay", "Wexford")) +
  
  # set map limits
  coord_sf(xlim = c(minx, maxx),
           ylim = c(miny, maxy), crs = 4326, expand = F) +
  
  annotation_scale(location = "bl", width_hint = 0.2, pad_y = unit(0.1, "in")) +
  annotation_north_arrow(location = "bl", which_north = "true",
                         pad_x = unit(0.1, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_orienteering,
                         height = unit(0.8, "cm"), width = unit(0.7, "cm"),) +
  
  # Plot styling
  theme_light() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.text=element_text(size=12), legend.position = "right",
        axis.title = element_blank(), legend.title= element_text(size=14),
        legend.box.background=element_rect(colour = "#BDC3C7"), legend.box.margin=margin(1,1,1,1),
        legend.spacing.x = unit(0.1, "cm"),  axis.text = element_text(size = 14))


## save the plot as a png
ggsave(plot = P1, 
       filename = "Paper Plots/Supp Fig 5- All nest sites.png",
       device = device,
       units = units, width = 175, height = 300, dpi = dpi,   
)





## Create map of all repeat nesting locations

## keep IDs with repeat years
Dups1 <- NestLocs[duplicated(NestLocs$device_id),]

DupNests <- NestLocs %>% filter(device_id %in% unique(Dups1$device_id))
DupLats <- med_nest_locs %>% filter(device_id %in% unique(Dups1$device_id))

## create the plot
P2 <- ggplot() + 
  
  # add the filled in countries
  geom_sf(data = countries2, aes(geometry = geometry), fill = "#BDC3C7", colour = "#BDC3C7", alpha = 0.6) +
  
  # add the nest locations
  geom_sf(data = DupNests, aes(geometry = geometry, colour = device_id),  size = 4, alpha = 0.7) + 
  scale_color_manual(name = "Tag ID Code", values=met.brewer("Lakota", as.numeric(length(unique(DupNests$device_id))))) +
  
  geom_path(data = DupLats, aes(x = breeding_long, y = breeding_lat, group = device_id), 
            alpha = 0.8, size = 0.5 ) +
  
  # set map limits
  coord_sf(xlim = c(minx, maxx),
           ylim = c(miny, maxy), crs = 4326, expand = F) +
  
  annotation_scale(location = "bl", width_hint = 0.2, pad_y = unit(0.1, "in")) +
  annotation_north_arrow(location = "bl", which_north = "true",
                         pad_x = unit(0.1, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_orienteering,
                         height = unit(0.8, "cm"), width = unit(0.7, "cm"),) +
  
  # Plot styling
  theme_light() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.text=element_text(size=12), legend.position = "right",
        axis.title = element_blank(), legend.title= element_text(size=14),
        legend.box.background=element_rect(colour = "#BDC3C7"), legend.box.margin=margin(1,1,1,1),
        legend.spacing.x = unit(0.1, "cm"), axis.text = element_text(size = 14))

## save the plot as a png
ggsave(plot = P2, 
       filename = "Paper Plots/Supp Fig 6- Repeat nest sites.png",
       device = device,
       units = units, width = 175, height = 300, dpi = dpi,   
)



