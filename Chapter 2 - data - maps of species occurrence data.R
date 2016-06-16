library(sp)
library(raster)
setwd("C:/data/modelling")
data.envi <- stack("data/envi/#AquaMaps.grd")
data.extent <- raster("data/envi/extent.grd")
binraster <- aggregate(data.extent, fact=5, fun="sum")
values(binraster) <- 0

# Records per cell ----

#OBIS records per cell
  load("data/species/data.obis.RData")
  spsdo <- SpatialPoints(cbind(data.obis$longitude, data.obis$latitude))
  binnedspsdo <- rasterize(spsdo, data.envi[[1]], fun="count", background=0)
  plot(binnedspsdo)
  rpcobis <- aggregate(binnedspsdo, fact=5, fun="sum")
  plot(rpcobis)
  writeRaster(rpcobis, "C:/data/modelling/output/arcGIS/species data/rpcobis.asc", overwrite=T)

#GBIF records per cell
  load("data/species/data.gbif.RData")
  spsdg <- SpatialPoints(cbind(data.gbif$longitude, data.gbif$latitude))
  binnedspsdg <- rasterize(spsdg, data.envi[[1]], fun="count", background=0)
  plot(binnedspsdg)
  rpcgbif <- aggregate(binnedspsdg, fact=5, fun="sum")
  plot(rpcgbif)
  writeRaster(rpcgbif, "C:/data/modelling/output/arcGIS/species data/rpcgbif.asc", overwrite=T)
  
#species data records per cell
  load("data/species/data.species.RData")
  spsd <- SpatialPoints(cbind(data.species$x, data.species$y))
  binnedspsd <- rasterize(spsd, data.envi[[1]], fun="count", background=0)
  plot(binnedspsd)
  rpc <- aggregate(binnedspsd, fact=5, fun="sum")
  plot(rpc)
  writeRaster(rpc, "C:/data/modelling/output/arcGIS/species data/rpc.asc", overwrite=T)

# Species per cell ----

#obis species per cell
  load("data/species/data.obis.RData")
  species="Engraulis ringens"
  spcobis=binraster
  for(species in unique(data.obis$TaxonName)){
    specrecso <- data.obis[data.obis$TaxonName==species,]
    spspecrecso <- SpatialPoints(cbind(specrecso$longitude, specrecso$latitude))
    spspecrecsobinned <- rasterize(spspecrecso, binraster, fun="count", background=0)
    spspecrecsobinned <- spspecrecsobinned>0
    spcobis <- spcobis + spspecrecsobinned
  }
  plot(spcobis)
  writeRaster(spcobis, file="C:/data/modelling/output/arcGIS/species data/spcobis.asc", overwrite=T)

#gbif species per cell
  load("data/species/data.gbif.RData")
  data.gbif <- data.gbif[!data.gbif$latitude==0 | !data.gbif$latitude==0,]
  data.gbif[data.gbif$latitude < 1 & data.gbif$latitude > -1 & data.gbif$longitude < 2 & data.gbif$longitude > -2, c("latitude","longitude")]

  species="Engraulis ringens"
  spcgbif=binraster
  for(species in unique(data.gbif$TaxonName)){
    specrecsg <- data.gbif[data.gbif$TaxonName==species,]
    spspecrecsg <- SpatialPoints(cbind(specrecsg$longitude, specrecsg$latitude))
    spspecrecsgbinned <- rasterize(spspecrecsg, binraster, fun="count", background=0)
    spspecrecsgbinned <- spspecrecsgbinned>0
    spcgbif <- spcgbif + spspecrecsgbinned
  }
  plot(spcgbif)
  writeRaster(spcgbif, file="C:/data/modelling/output/arcGIS/species data/spcgbif.asc", overwrite=T)

#species data species per cell
  load("data/species/data.species.RData")
  data.species <- data.species[!data.species$TaxonName=="none",]
  species="Engraulis ringens"
  spc=binraster
  for(species in unique(data.species$TaxonName)){
    specrecs <- data.species[data.species$TaxonName==species,]
    spspecrecs <- SpatialPoints(cbind(specrecs$x, specrecs$y))
    spspecrecsbinned <- rasterize(spspecrecs, binraster, fun="count", background=0)
    spspecrecsbinned <- spspecrecsbinned>0
    spc <- spc + spspecrecsbinned
  }
  plot(spc)
  writeRaster(spc, file="C:/data/modelling/output/arcGIS/species data/spc.asc", overwrite=T)

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#
# date recorded analysis ----
obisdates <- data.obis$datecollected[!data.obis$datecollected==""]
obisdates1 <- obisdates

#make "month/day/year" dates"year/month/day"
for(obisdaten in 1:length(obisdates)){
  obisdate <- obisdates[obisdaten]
  
  if(length(grep("*/*/....", obisdate))==1){
    print("reformatting")
    obisdates[obisdaten] <- paste(strsplit(obisdate, "/")[[1]][3], "-", strsplit(obisdate, "/")[[1]][1], "-", strsplit(obisdate, "/")[[1]][2], sep="")
  }else if(length(grep("..../../..", obisdate))==0){
    print(paste("non american date - ", obisdaten))    
  }
}
obisdates2 <- obisdates

#as good as it gets. some of these will be wrong!
for(obisdaten in 1:length(obisdates)){
  obisdate <- obisdates[obisdaten]
  
  if(length(grep("..-..-....", obisdate))==1 & strsplit(obisdate, "-")[[1]][2]>12){
    print("reformatting")
    obisdates[obisdaten] <- paste(strsplit(obisdate, "-")[[1]][3], "-", strsplit(obisdate, "-")[[1]][1], "-", strsplit(obisdate, "/")[[1]][2], sep="")
  }else if(length(grep("..../../..", obisdate))==0){
    print(paste("non american date - ", obisdaten))    
  }
}
obisdates3 <- obisdates

#
obisdates.new <- as.Date(obisdates3)
months <- as.integer(format(obisdates.new, "%m"))
summary(months)
table(months)
plot(months)

#################

summary(data.gbif$month)
table(data.gbif$month)
plot(data.gbif$month)







#make records "date / date / date"
obisdates <- gsub("-", "/", obisdates)
obisdates4 <- obisdates


#make records "dd/mm/yyyy"
errors=NULL
for(obisdaten in 1:length(obisdates)){
  obisdate <- obisdates[obisdaten]
  
  if(length(grep("./../....", obisdate))==1 & nchar(obisdate) < 10){
    print("reformatting")
    obisdates[obisdaten] <- paste("0", obisdate, sep="")
  }else if(length(grep(".././....", obisdate))==1){
    print("reformatting")
    obisdates[obisdaten] <- paste(strsplit(obisdate, "/")[[1]][1], "/0", strsplit(obisdate, "/")[[1]][2], "/", strsplit(obisdate, "/")[[1]][3], sep="")
  }else if(length(grep("././....", obisdate))==1){
    print("reformatting")
    obisdates[obisdaten] <- paste("0", strsplit(obisdate, "/")[[1]][1], "/0", strsplit(obisdate, "/")[[1]][2], "/", strsplit(obisdate, "/")[[1]][3], sep="")
  }else if(length(grep("../../....", obisdate))==1){
    print("already formatted")
  }else{
    print(paste("#something wrong here - ", obisdaten))
    errors <- c(errors, obisdaten)
  }
}
obisdates5 <- obisdates


as.Date(obisdates)



