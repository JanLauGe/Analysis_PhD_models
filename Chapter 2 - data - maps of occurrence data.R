library(sp)
library(raster)

occcount <- function(data.spec, data.extent){
  data.species.points <- SpatialPoints(cbind(data.species$x, data.species$y))
  occcount <- rasterize(data.species.points, data.extent, fun="count")
  data.final <- rasterToPoints(occcount, spatial=T)
  return(data.final)}

data.extent <- raster("C:/data/modelling/data/envi/data.extent.grd")

#Final set of occurrences
load("C:/data/modelling/data/species/data.species.RData")
final <- occcount(data.spec, data.extent)
plot(final)
