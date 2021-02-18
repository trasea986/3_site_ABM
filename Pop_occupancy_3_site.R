# script to determine occupancy

setwd("D:/OneDrive/GEM3_PostDoc/Agent-Based-Models/3_site_proj/outputs/dry_creek_sel")

library(tidyverse)
library(cowplot)
library(ggmap)
library(rgdal)
library(rgeos)
library(maptools)


data_counts = list.files(pattern = "ind*", 
                         full.names = TRUE, 
                         recursive = TRUE, 
                         include.dirs = TRUE) %>% 
  map_df(function(x) read_csv(x, col_types = cols(.default = "c")) %>% mutate(filename=gsub(".csv","",basename(x))))

#create year column from file name
data_df <- separate(data = data_counts, col = filename, into = c('Rest', 'Year'), sep = "ind")

#remove rows from the initialization, maybe year 1
data_df <- subset(data_df, data_df$Year != -1)
data_df$Year <- as.numeric(data_df$Year)
data_df$XCOORD <- as.numeric(data_df$XCOORD)
data_df$YCOORD <- as.numeric(data_df$YCOORD)

summary <- data_df %>%
  group_by(Year) %>%
  count(PatchID) %>%
  count(Year)


summary_clear <- summary %>% filter(n > 1)

summary_clear$Year <- as.integer(summary_clear$Year)

summary_clear %>%
  group_by(Year) %>%
  summarise(n = n) %>%
  mutate(freq = n / max(n))

#to then summarize

for_viz <- data_df %>%
group_by(PatchID, XCOORD, YCOORD, Year, .drop = FALSE) %>% 
  summarise(n = n())

#data_viz_test <- for_viz %>% filter(Year == 10)

#copied coordinates from ArcGIS Pro
#116.4234612°W 43.6655836°N 
#116.0662912°W 43.8163286°N 

map <- get_map(c(left = -116.4234612, bottom = 43.6655836, top = 43.8163286, right = -116.0662912))

p <- ggmap(map)

#change UTM to LONG LAT and add back in n
coordinates(for_viz) <- ~XCOORD+YCOORD #similar to SpatialPoints
data_viz_UTM <- SpatialPoints(for_viz, proj4string=CRS("+proj=utm +zone=11 +datum=WGS84")) 
data_viz_LONGLAT <- spTransform(data_viz_UTM, CRS("+proj=longlat +datum=WGS84"))
data_viz_ll_df <- as.data.frame(data_viz_LONGLAT)
data_viz_ll_df$n <- for_viz$n
data_viz_ll_df$Year <- for_viz$Year

#remove points for easier visualization



#add in stream layer data
stream <- readOGR("../../data_gis/NorWeST_PredictedStreamTempLines_MiddleSnake.shp")
stream_repro <- spTransform(stream, CRS("+proj=longlat +datum=WGS84"))
stream_fort <- fortify(stream_repro)

data_viz_ll_df <- data_viz_ll_df %>%
  filter(Year != 10)

for (i in unique(data_viz_ll_df$Year)) {

  Year = i 
  
  data_plot <- data_viz_ll_df %>%
    filter(Year == i)
  
  myplot <- p + 
    geom_path(data = stream_fort, aes(long, lat, group = group), color = "lightblue") +
  geom_point(data = data_plot, aes(x = XCOORD, y = YCOORD, size = n, alpha = n), color = "black", fill = "red", shape = 21) +
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