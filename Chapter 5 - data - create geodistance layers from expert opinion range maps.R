
#*#*#*#*#*#*#*#*#
# VARIABLES -----
#*#*#*#*#*#*#*#*#

library(raster)
library(gdistance)
library(dismo)
library(rgdal)

setwd("C:/data/modelling")


#*#*#*#*#*#
# DATA ----
#*#*#*#*#*#

#' Read range maps and convert to raster
species.list <- unlist(strsplit(list.files("C:/data/modelling/data/IUCN and SAUP/raster", pattern="*.grd"), split=".grd"))
  
#*#*#*#*#*#*#*#*#*#*#*#*#
# Calculate distance ----
#*#*#*#*#*#*#*#*#*#*#*#*#

#'create cost layer
  data.extent <- raster("C:/data/modelling/data/data.extent.mod.grd")
  costlayer <- transition(data.extent, transitionFunction=function(x){1},8,symm=FALSE)
  costlayer <- geoCorrection(costlayer)
  #proj4string(data.extent) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
  
s=species.list[1]
for(s in species.list){
  print(s)
  spec.raster <- raster(paste("C:/data/modelling/data/IUCN and SAUP/raster/", s, ".grd", sep=""))
    # points from grid cells that are presene in either of the two range maps
  spec.points <- rasterToPoints(spec.raster>0)[rasterToPoints(spec.raster>0)[,3]==1,c(1,2)]
  # calculate cost layer from points
  spec.geodist <- accCost(costlayer, spec.points)
  #plot(spec.geodist)
  writeRaster(spec.geodist, file=paste("C:/Users/LaurensG/Desktop/geodist/distraster/", s, ".grd", sep=""), overwrite=T)
  rm(spec.raster, spec.points, spec.geodist)
}

