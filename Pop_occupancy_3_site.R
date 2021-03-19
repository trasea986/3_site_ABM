# script to determine occupancy

setwd("../outputs")

library(tidyverse)
library(cowplot)
library(ggmap)
library(rgdal)
library(rgeos)
library(maptools)


output_files = list.files(pattern = paste(".csv", sep=''), 
                          full.names = TRUE, 
                          recursive = TRUE, 
                          include.dirs = TRUE) %>% 
  map_df(function(x) read_csv(x, col_types = cols(.default = "c")) %>% mutate(filename=paste(dirname((x)),basename(x),sep="/")))




#create year and run column from file name
data_df <- separate(data = output_files, col = filename, into = c('junk', 'run', 'junk1', 'junk2', 'year'), sep = "/")

data_df <- separate(data = data_df, col = year, into = c('ind', 'year'), sep = "d")
data_df <- separate(data = data_df, col = year, into = c('year', 'junk3'), sep = ".cs")
data_df <- separate(data = data_df, col = run, into = c('loc', 'run'), sep = "_")

#remove rows from the initialization, initial pop
data_df <- subset(data_df, data_df$year != -1)
data_df$junk <- NULL
data_df$junk1 <- NULL
data_df$junk2 <- NULL
data_df$junk3 <- NULL
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
data_df$loc <- as.factor(data_df$loc)
data_df$run <- as.factor(data_df$run)
data_df$year <- as.numeric(data_df$year)
data_df$XCOORD <- as.numeric(data_df$XCOORD)
data_df$YCOORD <- as.numeric(data_df$YCOORD)



#summarize for maps

for_viz <- data_df %>%
group_by(PatchID, XCOORD, YCOORD, year, loc, run, .drop = FALSE) %>% 
  summarise(n = n())

for_stat <- data_df %>%
  group_by(year, loc, run, .drop = FALSE) %>% 
  summarise(n = n())

#data_viz_test <- for_viz %>% filter(Year == 10)

#copied coordinates from ArcGIS Pro fro 
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
data_viz_LONGLAT <- spTransform(data_viz_UTM, CRS("+proj=longlat +datum=WGS84"))
data_viz_ll_df <- as.data.frame(data_viz_LONGLAT)
data_viz_ll_df$n <- for_viz$n
data_viz_ll_df$year <- for_viz$year
data_viz_ll_df$loc <- for_viz$loc
data_viz_ll_df$run <- for_viz$run

#sep by location
data_viz_ll_df_dry <- data_viz_ll_df %>%
  filter(loc == "dry")
data_viz_ll_df_cool <- data_viz_ll_df %>%
  filter(loc == "cool")
data_viz_ll_df_desert <- data_viz_ll_df %>%
  filter(loc == "desert")



#add in stream layer data
stream <- readOGR("../data_gis/NorWeST_PredictedStreamTempLines_MiddleSnake.shp")
stream_repro <- spTransform(stream, CRS("+proj=longlat +datum=WGS84"))
stream_fort <- fortify(stream_repro)



#dry creek test
for (i in unique(data_viz_ll_df_dry$year)) {
  year = i 
  
  data_plot <- data_viz_ll_df_dry %>%
    filter(year == i)
  
  myplot <- p_dry + 
    geom_path(data = stream_fort, aes(long, lat, group = group), color = "lightblue") +
  geom_point(data = data_plot, aes(x = XCOORD, y = YCOORD, size = n, alpha = n), color = "black", fill = "red", shape = 21) +
    facet_wrap(~run) +
  scale_size_continuous(range = c(1, 8)) +
  scale_alpha_continuous(range = c(0.25, .75)) +
  theme(legend.position = "none") +
  ylab("Latitude") +
  xlab("Longitude") +
  ggtitle(paste("Dry Creek",i,sep = " ")) +
  theme_bw(base_size = 14)+
  theme(legend.position = "none")
  
  ggsave(myplot, filename=paste("DryCreek",i,".png",sep=""),dpi = 400, width = 10, height = 8, units = "in")
}


#cool creek test
for (i in unique(data_viz_ll_df_cool$year)) {
  year = i 
  
  data_plot <- data_viz_ll_df_cool %>%
    filter(year == i)
  
  myplot <- p_cool + 
    geom_path(data = stream_fort, aes(long, lat, group = group), color = "lightblue") +
    geom_point(data = data_plot, aes(x = XCOORD, y = YCOORD, size = n, alpha = n), color = "black", fill = "red", shape = 21) +
    facet_wrap(~run) +
    scale_size_continuous(range = c(1, 8)) +
    scale_alpha_continuous(range = c(0.25, .75)) +
    theme(legend.position = "none") +
    ylab("Latitude") +
    xlab("Longitude") +
    ggtitle(paste("Keithley/Mann Creek",i,sep = " ")) +
    theme_bw(base_size = 14)+
    theme(legend.position = "none")
  
  ggsave(myplot, filename=paste("Cool",i,".png",sep=""),dpi = 400, width = 10, height = 8, units = "in")
}


#desert creek test
for (i in unique(data_viz_ll_df_desert$year)) {
  year = i 
  
  data_plot <- data_viz_ll_df_desert %>%
    filter(year == i)
  
  myplot <- p_desert + 
    geom_path(data = stream_fort, aes(long, lat, group = group), color = "lightblue") +
    geom_point(data = data_plot, aes(x = XCOORD, y = YCOORD, size = n, alpha = n), color = "black", fill = "red", shape = 21) +
    facet_wrap(~run) +
    scale_size_continuous(range = c(1, 8)) +
    scale_alpha_continuous(range = c(0.25, .75)) +
    theme(legend.position = "none") +
    ylab("Latitude") +
    xlab("Longitude") +
    ggtitle(paste("Jack's Creek",i,sep = " ")) +
    theme_bw(base_size = 14)+
    theme(legend.position = "none")
  
  ggsave(myplot, filename=paste("Desert",i,".png",sep=""),dpi = 400, width = 10, height = 8, units = "in")
}




#for animation comparing the dry and keithley
dry_cool_df <- rbind(data_viz_ll_df_cool, data_viz_ll_df_desert)

dry_cool_df <- dry_cool_df %>%
  filter(run == "sel-plast")

dry_cool_df$loc <- plyr::revalue(dry_cool_df$loc, c("cool"="Keithley/Mann", "desert"="Jacks' Creeks"))

for (i in unique(dry_cool_df$year)) {
  year = i 
  
  data_plot <- dry_cool_df %>%
    filter(year == i)
  
  myplot1 <- p_desert + 
    geom_path(data = stream_fort, aes(long, lat, group = group), color = "lightblue") +
    geom_point(data = data_plot, aes(x = XCOORD, y = YCOORD, size = n, alpha = n), color = "black", fill = "red", shape = 21) +
    scale_size_continuous(range = c(1, 8)) +
    scale_alpha_continuous(range = c(0.25, .75)) +
    theme(legend.position = "none") +
    ylab("Latitude") +
    xlab("Longitude") +
    ggtitle(paste("Jacks' Creeks, Year",i,sep = " ")) +
    theme_bw(base_size = 14)+
    theme(legend.position = "none")
  
  myplot2 <- p_cool + 
    geom_path(data = stream_fort, aes(long, lat, group = group), color = "lightblue") +
    geom_point(data = data_plot, aes(x = XCOORD, y = YCOORD, size = n, alpha = n), color = "black", fill = "red", shape = 21) +
    scale_size_continuous(range = c(1, 8)) +
    scale_alpha_continuous(range = c(0.25, .75)) +
    theme(legend.position = "none") +
    ylab("Latitude") +
    xlab("Longitude") +
    ggtitle(paste("Keithley/Mann Creeks, Year",i,sep = " ")) +
    theme_bw(base_size = 14)+
    theme(legend.position = "none")

  myplots <- plot_grid(myplot1, myplot2)
    
  ggsave(myplots, filename=paste("Cool-Desert",i,".png",sep=""),dpi = 600, width = 10, height = 8, units = "in")
}
