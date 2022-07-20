# script to determine occupancy



library(tidyverse)
library(cowplot)
library(ggmap)
#library(rgdal)
library(rgeos)
library(maptools)

#rgeos/rgdal being depricated, so trying sf
#library(sf)

#read in input file to pull in temperature values for each of patches. this does not vary across scenarios.

setwd("./inputs/null")

patchvars_cool <- read.csv("patchvars_cool.csv")
patchvars_desert <- read.csv("patchvars_desert.csv")
patchvars_dry <- read.csv("patchvars_dry.csv")
patchvars <- rbind(patchvars_cool, patchvars_desert, patchvars_dry)

#need to subset for left_join memory issues
patchvars <- patchvars %>%
  select(X, Y, GrowthTemperatureBack)

#then a little renaming
patchvars$XCOORD <- patchvars$X
patchvars$YCOORD <- patchvars$Y
patchvars$Y <- NULL
patchvars$X <- NULL

setwd("../../../outputs/3_site")

#for bringing in files, important to note the directory strucute:
#outputs/3_site/cool_null ; cool_plast; so on
#also make sure no summary files still present

output_files = list.files(pattern = paste(".csv", sep=''), 
                          full.names = TRUE, 
                          recursive = TRUE, 
                          include.dirs = TRUE) %>% 
  map_df(function(x) read_csv(x, col_types = cols(.default = "c")) %>% mutate(filename=paste(dirname((x)),basename(x),sep="/")))


#data_df <- head(output_files) this is useful for testing below when running lots of replicates
data_df <- output_files

#create year and run column from file name after deleting unneccessary rows
data_df$ID <- NULL
data_df$sex <- NULL
data_df$size <- NULL
data_df$mature <- NULL
data_df$newmature <- NULL
data_df$layeggs <- NULL
data_df$capture <- NULL
data_df$infection <- NULL
data_df$Hindex <- NULL
data_df$Species <- NULL
data_df$recapture <- NULL
data_df$SubPatchID <- NULL
data_df$age <- NULL
data_df$CDist <- NULL

data_df <- separate(data = output_files, col = filename, into = c('junk', 'run', 'replicate', 'mc', 'year'), sep = "/")

data_df <- separate(data = data_df, col = year, into = c('ind', 'year'), sep = "d")
data_df <- separate(data = data_df, col = year, into = c('year', 'junk3'), sep = ".cs")
data_df <- separate(data = data_df, col = run, into = c('loc', 'run'), sep = "_")

#remove rows from the initialization, initial pop
data_df <- subset(data_df, data_df$year != -1)
data_df$junk <- NULL
data_df$junk1 <- NULL
data_df$junk2 <- NULL
data_df$junk3 <- NULL
data_df$loc <- as.factor(data_df$loc)
data_df$run <- as.factor(data_df$run)
data_df$year <- as.numeric(data_df$year)
data_df$XCOORD <- as.numeric(data_df$XCOORD)
data_df$YCOORD <- as.numeric(data_df$YCOORD)
data_df$mc <- as.factor(data_df$mc)

saveRDS(data_df, file = "data_df.Rds")

#summarize for maps
for_viz <- data_df %>%
group_by(PatchID, XCOORD, YCOORD, year, loc, run, replicate, mc, .drop = FALSE) %>% 
  summarise(n = n())

for_viz1 <- for_viz %>%
  group_by(PatchID, XCOORD, YCOORD, year, loc, run, .drop = FALSE) %>%
  summarize(m = mean(n), sdev = sd(n))

#bring in temperature column from patchvars
for_viz1 <- merge(x=for_viz1, y=patchvars, by=c("XCOORD", "YCOORD"))

saveRDS(for_viz1, file = "for_viz.Rds")

for_stat <- data_df %>%
  group_by(year, loc, run, replicate, mc, .drop = FALSE) %>% 
  summarise(n = n())

saveRDS(for_stat, file = "for_stat.Rds")

pop <- data_df %>%
  group_by(year, loc, run, replicate, mc) %>%
  tally()

#summary stats more useful for plot
pop1 <- pop %>%
  group_by(year, loc, run, .drop = FALSE) %>%
  summarize(m = mean(n), sdev = sd(n))

saveRDS(pop1, file = "pop1.Rds")

#LOAD summary info
pop <- readRDS('pop1.RDS')
#LOAD by site
for_viz <- readRDS('for_viz.RDS')

#from here down can be done locally (although data_df may be too large)

#copied AOI coordinates from ArcGIS Pro for 
#Dry Creek
#116.4234612?W 43.6655836?N 
#116.0662912?W 43.8163286?N 

#Mann/Keithley
#117.0561113°W 44.3861876°N BL
#116.7206293°W 44.6363213°N TR

#Jacks
#115.8251393°W 42.8174436°N  TR
#116.3219521°W 42.3691728°N  BL

map_dry_creek <- get_map(c(left = -116.4234612, bottom = 43.6655836, top = 43.8163286, right = -116.0662912))

map_cool <- get_map(c(left = -117.0561113, bottom = 44.3861876, top = 44.6363213, right = -116.7206293))

map_desert <- get_map(c(left = -116.3219521, bottom = 42.3691728, top = 42.8174436, right = -115.8251393))

p_dry <- ggmap(map_dry_creek)
p_cool <- ggmap(map_cool)
p_desert <- ggmap(map_desert)

p_dry
p_cool
p_desert

#change UTM to LONG LAT and add back in n
coordinates(for_viz) <- ~XCOORD+YCOORD #similar to SpatialPoints
data_viz_UTM <- SpatialPoints(for_viz, proj4string=CRS("+proj=utm +zone=11 +datum=WGS84"))
data_viz_LONGLAT <- spTransform(data_viz_UTM, CRS("+proj=longlat +datum=WGS84")) #can also be done with spTransform but being depricated
data_viz_ll_df <- as.data.frame(data_viz_LONGLAT)
data_viz_ll_df$m <- for_viz$m
data_viz_ll_df$year <- for_viz$year
data_viz_ll_df$loc <- for_viz$loc
data_viz_ll_df$run <- for_viz$run
data_viz_ll_df$replicate <- for_viz$replicate
data_viz_ll_df$GrowthTemperatureBack <- for_viz$GrowthTemperatureBack

#sep by location
data_viz_ll_df_dry <- data_viz_ll_df %>%
  filter(loc == "dry")
data_viz_ll_df_cool <- data_viz_ll_df %>%
  filter(loc == "cool")
data_viz_ll_df_desert <- data_viz_ll_df %>%
  filter(loc == "desert")

#one problem with points is that it is hard to see due to overlap. remove points
data_viz_ll_df_desert <- data_viz_ll_df_desert %>%
  group_by(run) %>%
  sample_frac(0.2, replace = FALSE)
data_viz_ll_df_cool <- data_viz_ll_df_cool %>%
  group_by(run) %>%
  sample_frac(0.2, replace = FALSE)
data_viz_ll_df_dry <- data_viz_ll_df_dry %>%
  group_by(run) %>%
  sample_frac(0.2, replace = FALSE)

#next, we are going to single out the first (10) and last years (100) of the model
data_viz_ll_df_desert <- data_viz_ll_df_desert %>% 
  filter(year == 10 | year == 100)
data_viz_ll_df_cool <- data_viz_ll_df_cool %>% 
  filter(year == 10 | year == 100)
data_viz_ll_df_dry <- data_viz_ll_df_dry %>% 
  filter(year == 10 | year == 100)


#now to split the temp up based on the input requirements
#some of the inputs have " | " and some have "|" and some have "  |  ". Also, seem to glitch when separating by | so change to "/" before that step.

data_viz_ll_df_desert$GrowthTemperatureBack <- gsub(pattern="|", replacement=" | ", data_viz_ll_df_desert$GrowthTemperatureBack, fixed = TRUE)
data_viz_ll_df_desert$GrowthTemperatureBack <- gsub(pattern="  |  ", replacement=" | ", data_viz_ll_df_desert$GrowthTemperatureBack, fixed = TRUE)
data_viz_ll_df_desert$GrowthTemperatureBack <- gsub(pattern=" | ", replacement="/", data_viz_ll_df_desert$GrowthTemperatureBack, fixed = TRUE)
data_viz_ll_df_desert<- separate(data_viz_ll_df_desert, col = GrowthTemperatureBack, into = c('C10', 'mid', 'C100'), sep = "/")

data_viz_ll_df_dry$GrowthTemperatureBack <- gsub(pattern="|", replacement=" | ", data_viz_ll_df_dry$GrowthTemperatureBack, fixed = TRUE)
data_viz_ll_df_dry$GrowthTemperatureBack <- gsub(pattern="  |  ", replacement=" | ", data_viz_ll_df_dry$GrowthTemperatureBack, fixed = TRUE)
data_viz_ll_df_dry$GrowthTemperatureBack <- gsub(pattern=" | ", replacement="/", data_viz_ll_df_dry$GrowthTemperatureBack, fixed = TRUE)
data_viz_ll_df_dry<- separate(data_viz_ll_df_dry, col = GrowthTemperatureBack, into = c('C10', 'mid', 'C100'), sep = "/")

data_viz_ll_df_cool$GrowthTemperatureBack <- gsub(pattern="|", replacement=" | ", data_viz_ll_df_cool$GrowthTemperatureBack, fixed = TRUE)
data_viz_ll_df_cool$GrowthTemperatureBack <- gsub(pattern="  |  ", replacement=" | ", data_viz_ll_df_cool$GrowthTemperatureBack, fixed = TRUE)
data_viz_ll_df_cool$GrowthTemperatureBack <- gsub(pattern=" | ", replacement="/", data_viz_ll_df_cool$GrowthTemperatureBack, fixed = TRUE)
data_viz_ll_df_cool<- separate(data_viz_ll_df_cool, col = GrowthTemperatureBack, into = c('C10', 'mid', 'C100'), sep = "/")

#change to numeric and delete unwanted column
data_viz_ll_df_desert$mid <- NULL
data_viz_ll_df_desert$C10 <- as.numeric(data_viz_ll_df_desert$C10)
data_viz_ll_df_desert$C100 <- as.numeric(data_viz_ll_df_desert$C100)

data_viz_ll_df_cool$mid <- NULL
data_viz_ll_df_cool$C10 <- as.numeric(data_viz_ll_df_cool$C10)
data_viz_ll_df_cool$C100 <- as.numeric(data_viz_ll_df_cool$C100)

data_viz_ll_df_dry$mid <- NULL
data_viz_ll_df_dry$C10 <- as.numeric(data_viz_ll_df_dry$C10)
data_viz_ll_df_dry$C100 <- as.numeric(data_viz_ll_df_dry$C100)

#next step is to then give the numeric ranges the threshold values
data_viz_ll_df_desert$C10_thresh  <- cut(data_viz_ll_df_desert$C10, breaks = c(-50, 18, 20, 50), labels = c("Cold", "Trigger", "Avoid"))
data_viz_ll_df_desert$C100_thresh <- cut(data_viz_ll_df_desert$C100, breaks = c(-50, 18, 20, 50), labels = c("Cold", "Trigger", "Avoid"))

data_viz_ll_df_dry$C10_thresh  <- cut(data_viz_ll_df_dry$C10, breaks = c(-50, 18, 20, 50), labels = c("Cold", "Trigger", "Avoid"))
data_viz_ll_df_dry$C100_thresh  <- cut(data_viz_ll_df_dry$C100, breaks = c(-50, 18, 20, 50), labels = c("Cold", "Trigger", "Avoid"))

data_viz_ll_df_cool$C10_thresh  <- cut(data_viz_ll_df_cool$C10, breaks = c(-50, 18, 20, 50), labels = c("Cold", "Trigger", "Avoid"))
data_viz_ll_df_cool$C100_thresh  <- cut(data_viz_ll_df_cool$C100, breaks = c(-50, 18, 20, 50), labels = c("Cold", "Trigger", "Avoid"))




#Now to make the maps for the filtered years. see loop on bottom, but problem is getting alignment/only using certain combinations of C and Year.

data_viz_ll_df_desert_10 <- data_viz_ll_df_desert %>%
  filter(year == 10)

JacksCreek10 <- p_desert + 
  geom_path(data = stream_fort, aes(long, lat, group = group), color = "steelblue2") +
  geom_point(data = data_viz_ll_df_desert_10, aes(x = XCOORD, y = YCOORD, size = m, fill = C10_thresh), color = "black", shape = 21) +
  facet_wrap(~run) +
  scale_size_continuous(range = c(2, 6)) +
  scale_fill_manual(values = c("blue", "orange", "red")) +
  scale_alpha_continuous(range = c(0.15, .85)) +
  theme(legend.position = "none") +
  ylab("Latitude") +
  xlab("Longitude") +
  ggtitle(paste("Jack's Creek","10",sep = " ")) +
  theme_bw(base_size = 14)+
  theme(legend.position = "none")

ggsave(
  JacksCreek10, filename=paste("JacksCreek","10",".png",sep=""),dpi = 400, width = 10, height = 8, units = "in")

data_viz_ll_df_desert_100 <- data_viz_ll_df_desert %>%
  filter(year == 100)

JacksCreek100 <-  p_desert + 
  geom_path(data = stream_fort, aes(long, lat, group = group), color = "steelblue2") +
  geom_point(data = data_viz_ll_df_desert_100, aes(x = XCOORD, y = YCOORD, size = m, fill = C100_thresh), color = "black", shape = 21) +
  facet_wrap(~run) +
  scale_size_continuous(range = c(2, 6)) +
  scale_fill_manual(values = c("blue", "orange", "red")) +
  scale_alpha_continuous(range = c(0.15, .85)) +
  theme(legend.position = "none") +
  ylab("Latitude") +
  xlab("Longitude") +
  ggtitle(paste("Jack's Creek","100",sep = " ")) +
  theme_bw(base_size = 14)+
  theme(legend.position = "none")

ggsave(JacksCreek100, filename=paste("JacksCreek","100",".png",sep=""),dpi = 400, width = 10, height = 8, units = "in")


data_viz_ll_df_dry_10 <- data_viz_ll_df_dry %>%
  filter(year == 10)

DryCreek10 <- p_dry + 
  geom_path(data = stream_fort, aes(long, lat, group = group), color = "steelblue2") +
  geom_point(data = data_viz_ll_df_dry_10, aes(x = XCOORD, y = YCOORD, size = m, fill = C10_thresh), color = "black", shape = 21) +
  facet_wrap(~run) +
  scale_size_continuous(range = c(2, 6)) +
  scale_fill_manual(values = c("blue", "orange", "red")) +
  scale_alpha_continuous(range = c(0.15, .85)) +
  theme(legend.position = "none") +
  ylab("Latitude") +
  xlab("Longitude") +
  ggtitle(paste("Dry Creek","10",sep = " ")) +
  theme_bw(base_size = 14)+
  theme(legend.position = "none")
ggsave(
  DryCreek10, filename=paste("DryCreek","10",".png",sep=""),dpi = 400, width = 10, height = 8, units = "in")

data_viz_ll_df_dry_100 <- data_viz_ll_df_dry %>%
  filter(year == 100)

DryCreek100 <-  p_dry + 
  geom_path(data = stream_fort, aes(long, lat, group = group), color = "steelblue2") +
  geom_point(data = data_viz_ll_df_dry_100, aes(x = XCOORD, y = YCOORD, size = m, fill = C100_thresh), color = "black", shape = 21) +
  facet_wrap(~run) +
  scale_size_continuous(range = c(2, 6)) +
  scale_fill_manual(values = c("blue", "orange", "red")) +
  scale_alpha_continuous(range = c(0.15, .85)) +
  theme(legend.position = "none") +
  ylab("Latitude") +
  xlab("Longitude") +
  ggtitle(paste("Dry Creek","100",sep = " ")) +
  theme_bw(base_size = 14)+
  theme(legend.position = "none")

ggsave(DryCreek100, filename=paste("DryCreek","100",".png",sep=""),dpi = 400, width = 10, height = 8, units = "in")



data_viz_ll_df_cool_10 <- data_viz_ll_df_cool %>%
  filter(year == 10)

KMCreek10 <- p_cool + 
  geom_path(data = stream_fort, aes(long, lat, group = group), color = "steelblue2") +
  geom_point(data = data_viz_ll_df_cool_10, aes(x = XCOORD, y = YCOORD, size = m, fill = C10_thresh), color = "black", shape = 21) +
  facet_wrap(~run) +
  scale_size_continuous(range = c(2, 6)) +
  scale_fill_manual(values = c("blue", "orange", "red")) +
  scale_alpha_continuous(range = c(0.15, .85)) +
  theme(legend.position = "none") +
  ylab("Latitude") +
  xlab("Longitude") +
  ggtitle(paste("Keithley/Mann Creek","10",sep = " ")) +
  theme_bw(base_size = 14)+
  theme(legend.position = "none")
ggsave(
  KMCreek10, filename=paste("KMCreek","10",".png",sep=""),dpi = 400, width = 10, height = 8, units = "in")


data_viz_ll_df_cool_100 <- data_viz_ll_df_cool %>%
  filter(year == 100)

KMCreek100 <-  p_cool + 
  geom_path(data = stream_fort, aes(long, lat, group = group), color = "steelblue2") +
  geom_point(data = data_viz_ll_df_cool_100, aes(x = XCOORD, y = YCOORD, size = m, fill = C100_thresh), color = "black", shape = 21) +
  facet_wrap(~run) +
  scale_size_continuous(range = c(2, 6)) +
  scale_fill_manual(values = c("blue", "orange", "red")) +
  scale_alpha_continuous(range = c(0.15, .85)) +
  theme(legend.position = "none") +
  ylab("Latitude") +
  xlab("Longitude") +
  ggtitle(paste("Keithley/Mann","100",sep = " ")) +
  theme_bw(base_size = 14)+
  theme(legend.position = "none")

ggsave(KMCreek100, filename=paste("KMCreek","100",".png",sep=""),dpi = 400, width = 10, height = 8, units = "in")




#### Below is the loop for the maps. The issue here is that the C value needs to change depending on the year value


#for comparing the desert and keithley
#dry_cool_df <- rbind(data_viz_ll_df_cool, data_viz_ll_df_desert)

#dry_cool_df <- dry_cool_df %>%
#  filter(run == "tol-plast")

#dry_cool_df$loc <- plyr::revalue(dry_cool_df$loc, c("cool"="Keithley/Mann", "desert"="Jacks' Creeks"))



#the below loop is useful if wanting every output as a map
#for (i in unique(dry_cool_df$year)) {
#  year = i 
#  
#  data_plot <- dry_cool_df %>%
#    filter(year == i)
#  
#  myplot1 <- p_desert + 
#    geom_path(data = stream_fort, aes(long, lat, group = group), color = "lightblue") +
#    geom_point(data = data_plot, aes(x = XCOORD, y = YCOORD, size = n, alpha = n), color = "black", fill #= "red", shape = 21) +
#    scale_size_continuous(range = c(1, 8)) +
#    scale_alpha_continuous(range = c(0.25, .75)) +
#    theme(legend.position = "none") +
#    ylab("Latitude") +
#    xlab("Longitude") +
#   ggtitle(paste("Jacks' Creeks, Year",i,sep = " ")) +
#    theme_bw(base_size = 14)+
#    theme(legend.position = "none")
#  
#  myplot2 <- p_cool + 
#    geom_path(data = stream_fort, aes(long, lat, group = group), color = "lightblue") +
#    geom_point(data = data_plot, aes(x = XCOORD, y = YCOORD, size = n, alpha = n), color = "black", fill #= "red", shape = 21) +
#    scale_size_continuous(range = c(1, 8)) +
#    scale_alpha_continuous(range = c(0.25, .75)) +
#    theme(legend.position = "none") +
#    ylab("Latitude") +
#    xlab("Longitude") +
#    ggtitle(paste("Keithley/Mann Creeks, Year",i,sep = " ")) +
#    theme_bw(base_size = 14)+
#    theme(legend.position = "none")
#  
#  myplots <- plot_grid(myplot1, myplot2)
#  
#  ggsave(myplots, filename=paste("Cool-Desert",i,".png",sep=""),dpi = 600, width = 10, height = 8, units #= "in")
#}