# LIBRARIES ####
library(sp)
library(raster)
library(rgdal)
library(dismo)
library(plyr)

#*#*#*#*#*#*#*#*#
# VARIABLES -----
#*#*#*#*#*#*#*#*#

setwd("C:/data/modelling")
species.list <- read.csv("C:/data/modelling/data/speciesinfo.csv", header=T)
data.extent <- raster("C:/data/modelling/data/data.extent.grd")

#limit species with IUCN and SAUP range map
species.expertmap <- unlist(strsplit(list.files("C:/data/modelling/data/IUCN/singleshapefiles", pattern=".shp$"), split=".shp"))
species.saupmap <- unlist(strsplit(list.files("C:/data/modelling/data/SAUP/singleshapefiles", pattern=".shp$"), split=".shp"))
species.list.s <- species.list[species.list$TaxonKey %in% species.expertmap & species.list$TaxonKey %in% species.saupmap,"TaxonKey"]

#get background data
#load("data/data.species.RData")
#data.sample <- data.species[data.species$set=="background",c("x", "y")]


#' The Jaccard coefficient measures similarity between sample sets, and is 
#' defined as the size of the intersection divided by the size of the union of 
#' the sample sets. Function assumes that values in rasters being compared
#' have values [0, 1] (other behaviour is not defined), such as rank rasters
#' produced by Zonation. 
#' 
#limit species with IUCN and SAUP range map
species.expertmap <- unlist(strsplit(list.files("C:/data/modelling/data/IUCN/singleshapefiles", pattern=".shp$"), split=".shp"))
species.saupmap <- unlist(strsplit(list.files("C:/data/modelling/data/SAUP/singleshapefiles", pattern=".shp$"), split=".shp"))
species.list.s <- species.list[species.list$TaxonKey %in% species.expertmap & species.list$TaxonKey %in% species.saupmap,"TaxonKey"]

jaccard=NULL
for(s in 1:length(species.list.s)){
  print(species.list.s[s])
  spec.raster.exp <- raster(paste("data/IUCN/raster/", species.list.s[s], sep=""))>0
  spec.raster.saup <- raster(paste("data/SAUP/raster/", species.list.s[s], sep=""), crs=crs("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))>0
  
  #Calculate Jaccard index for comparison of the two rasters
  spec.rasters <- spec.raster.exp + spec.raster.saup
  spec.jaccard <- length(spec.rasters[spec.rasters==2]) / length(spec.rasters[spec.rasters >=1])
  jaccard <- rbind(jaccard, spec.jaccard)
  rm(spec.raster.exp, spec.raster.saup, spec.rasters)}

rownames(jaccard) <- species.list.s

plot(jaccard)
