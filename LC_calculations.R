#Script for Least Cost calculations

library(raster)
library(gdistance)
library(geosphere)
library(rgdal)
library(tidyverse)



#bring in raster reclassified file
south_stream <- raster("sstream_f.asc")
cool_dry_streams <- raster("dry_cool_res.asc")
plot(south_stream)
plot(cool_dry_streams)

#bring in points
Patches <- read.csv("P_100m.csv")
cool_patches <- read.csv("cool_basin_temps.csv")
dry_patches <- read.csv("dry_creek_temps.csv")


#create transition object

tr1 <- transition(1/south_stream, transitionFunction=mean, directions=8)
tr1 <- geoCorrection(tr1, type="c")

tr2 <- transition(1/cool_dry_streams, transitionFunction=mean, directions=8)
tr2 <- geoCorrection(tr2, type="c")

#need to now create point file, which is same as above but just the x and y

South_LC_Points <- Patches[-c(1,1)]
#with second sites I'm using tidyverse
cool_patches_LC <- cool_patches %>% select("POINT_X", "POINT_Y")
dry_patches_LC <- dry_patches %>% select("POINT_X", "POINT_Y")

#now to calculate distance matrix. first part is based on just euclidean distance, which I did before but am not including here

#CSF_sp_latlong2 <- CSF_sp_latlong[,c(2,1)]
#eucdist <- distm(SpatialPoints(as.matrix(CSF_sp_latlong2)), fun=distGeo)

#costDist seemed to be weird, next up option is shortest path like the image.
south_distance2 <- costDistance(tr1, as.matrix(South_LC_Points))

dry_distance2 <- costDistance(tr2, as.matrix(dry_patches_LC))
cool_distance2 <- costDistance(tr2, as.matrix(cool_patches_LC))


write.csv(as.matrix(cool_distance2), "cool_dist2.csv")
write.csv(as.matrix(dry_distance2), "dry_dist2.csv")

write.csv(as.matrix(south_distance2), "south_lc_dist2_inv.csv")





#display results section, with a few options I've kept for use later. Note: has external dependencies


plot(south_stream)
#plot(SpatialPoints(CSF_sp), add=TRUE, pch=20, col="red")
#text(CSF_sp[,1]+2, CSF_sp[,2]+2, 1:nrow(CSF_sp))
#too many points, subsetting CSF to only first 50 sites
sp_sub <- South_LC_Points[-c(51:2570),]
plot(SpatialPoints(sp_sub), add=TRUE, pch=20, col="blue")



for (i in 1:50) {
  plot(shortestPath(tr1, as.matrix(sp_sub[1,]), as.matrix(sp_sub[i,]), output="SpatialLines"), add=TRUE)
}

#make a movie
#to install ImageMagick


library(animation)

saveGIF(
  {for (i in 1:50) {
    p <- plot(SA_2016)
    p <- plot(SpatialPoints(CSF_sp), add=TRUE, pch=20, col="red")
    p <- plot(shortestPath(tr1, as.matrix(CSF_sp[1,]), as.matrix(CSF_sp[i,]), output="SpatialLines"), add=TRUE)}
    print(p)
  },
  movie.name = "CSF_2016.gif", interval = .2, nmax =674, ani.width = 600, ani.height = 600,
  outdir = getwd()
)
