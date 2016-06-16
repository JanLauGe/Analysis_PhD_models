library(raster)
library(rasterVis)
library(rgdal)

species.list <- read.csv("G:/HPC/submission/finalmodel/data/species.list.csv")
data.extent <- raster("G:/HPC/submission/finalmodel/data/data.extent.grd")
data.land <- readOGR("C:/data/modelling/data/envi/othershapefiles", layer="land")
setwd("G:/HPC/output/finalmodel")

for(s in 1:length(species.list)){
  try(load(paste("SDMeval_clipped_finalmodel_", species.list[s,"TaxonKey"], ".RData", sep="")))
  
  predBioclim = data.extent
  predMaxent = data.extent
  predBRT = data.extent
  
  for(k in 1:10){
    predBioclim <- predBioclim + results[[k]][[2]][[1]]
    predMaxent <- predMaxent + results[[k]][[2]][[2]]
    predBRT <- predBRT + results[[k]][[2]][[3]]
  }
  
  rm(results)
}

tiff("C:/Users/LaurensG/Desktop/examplemap.tiff", height=10, width=16, unit="cm", res=300)
levelplot(predMaxent, margin=F, contour=F, colorkey=list(space="right"), pretty=T, par.settings=BuRdTheme, main=paste(species.list[s,"TaxonName"], " - ", "Maxent models")) +
  layer(sp.polygons(data.land, fill="black"))
dev.off()



library(gridExtra)
library(rgdal)
library(maptools)
library(plyr)
library(rasterVis)
library(raster)

species.list <- read.csv("G:/HPC/submission/finalmodel/data/species.list.csv")
data.extent <- raster("G:/HPC/submission/finalmodel/data/data.extent.grd")
data.land <- readOGR("C:/data/modelling/data/envi/othershapefiles", layer="land")
setwd("G:/HPC/output/finalmodel")

TaxonKeys=c(600143,601342,605367,600656,600138)
load("C:/data/modelling/data/data.species.RData")

s=5
for(s in 1:5){
  TaxonKey = TaxonKeys[s]
  records.species <- data.species[data.species$TaxonKey==TaxonKey,c("x","y")]
  
  Plot <- levelplot(data.extent, margin=F, contour=F, colorkey=F, pretty=T, par.settings=BuRdTheme, main=unique(data.species[data.species$TaxonKey==TaxonKey,"TaxonName"])) + 
    layer(sp.polygons(data.land, fill="black"), columns=1) +
    layer(sp.points(SpatialPoints(get(paste("records.species",s,sep=""))), pch=20, cex=.5, col="red"))
  
  tiff(paste("C:/Users/LaurensG/Desktop/plots/record maps/plotmap",s,".tiff",sep=""), height=10, width=15, unit="cm", res=300)
    Plot
  dev.off()
}
grid.arrange(Plot1, Plot2, Plot3, Plot4, Plot5, nrow=4, ncol=2)

