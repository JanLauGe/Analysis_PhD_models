#Create plot of environmental predictor variables

library(raster)
library(rasterVis)
library(rgdal)
library(gridExtra)
library(RColorBrewer)

setwd("C:/data/modelling/")
data.envi <- stack("C:/data/modelling/data/data.envi.grd")
data.land <- readOGR("C:/data/modelling/data/envi/othershapefiles", layer="land")

cols <- brewer.pal(10, "Greens")

library(colorspace)
PlotTheme=rasterTheme(region=sequential_hcl(10, power=2.2))
levelplot(data.envi[[1]], par.settings=PlotTheme)

PlotTheme=rasterTheme(colorRamp(colors=c("yellow", "green"), space="Lab", interpolate="linear"))


var1 <- 
  levelplot(data.envi[[1]]/1000, margin=F, contour=F, colorkey=list(space="right"), 
  par.settings=BuRdTheme, main="Bathymetry", xlab=NULL, 
  scales=list(x=list(draw=FALSE))) + layer(sp.polygons(data.land, fill="black"))

var2 <- 
  levelplot(data.envi[[6]]/1000, main="Primary Productivity",
  par.settings=rasterTheme(region=rainbow(12)[4:12]),
  at=c(0,0.2,0.4,0.8,1.6,2.4,4,8), pretty=T, margin=F, contour=F, xlab=NULL, ylab=NULL, 
  scales=list(draw=FALSE)) + layer(sp.polygons(data.land, fill="black"))

var3 <- 
  levelplot(data.envi[[2]], margin=F, contour=F, colorkey=list(space="right"),
  #par.settings=rasterTheme(region=brewer.pal('Greens', n=9)[3:9]),
  par.settings=rasterTheme(region=topo.colors(10)),
  main="Sea Surface Temperature", xlab=NULL,
  scales=list(x=list(draw=FALSE))) + layer(sp.polygons(data.land, fill="black"))

var4 <- 
  levelplot(data.envi[[3]], margin=F, contour=F, colorkey=list(space="right"), 
  par.settings=BuRdTheme, main="Sea Bottom Temperature", xlab=NULL, ylab=NULL, 
  scales=list(draw=FALSE)) + layer(sp.polygons(data.land, fill="black"))

var5 <- 
  levelplot(data.envi[[4]], margin=F, contour=F, colorkey=list(space="right"), 
  par.settings=BuRdTheme, main="Salinity", xlab=NULL, 
  scales=list(x=list(draw=FALSE))) + layer(sp.polygons(data.land, fill="black"))

var6 <- 
  levelplot(data.envi[[5]], margin=F, contour=F, colorkey=list(space="right"), 
  par.settings=BuRdTheme, main="Ice Concentration", ylab=NULL, 
  scales=list(y=list(draw=FALSE))) + layer(sp.polygons(data.land, fill="black"))

var7 <- 
  levelplot(data.envi[[7]], margin=F, contour=F, colorkey=list(space="right"), 
  par.settings=BuRdTheme, main="Distance from Land") +
  layer(sp.polygons(data.land, fill="black"))

grid.arrange(var1, var2, var3, var4, var5, var6, var7, nrow=4, ncol=2)
####
#
#
#
tiff("C:/Users/LaurensG/Desktop/examplemap.tiff", height=10, width=16, unit="cm", res=300)
levelplot(predMaxent, margin=F, contour=F, colorkey=list(space="right"), pretty=T, par.settings=BuRdTheme, main=paste(species.list[s,"TaxonName"], " - ", "Maxent models")) +
  layer(sp.polygons(data.land, fill="black"))
dev.off()



