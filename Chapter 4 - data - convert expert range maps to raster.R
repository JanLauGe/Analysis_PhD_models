species.list.s <- species.list.s3
path <- "C:/data/modelling/data/IUCN and SAUP/"
s=1


for(s in 1:length(species.list.s1)){
  IUCN <- readOGR("data/IUCN/singleshapefiles", species.list.s1[s])
  raster.IUCN <- raster(paste("data/IUCN/raster/", species.list.s1[s], sep=""))
    
  png(paste("data/IUCN/maps/", species.list.s1[s], ".png", sep=""))
  plot(raster.IUCN)
  plot(IUCN, border="red", add=TRUE)
  dev.off()
}

for(s in 1:length(species.list.s2)){
  SAUP <- readOGR("data/SAUP/singleshapefiles", species.list.s2[s])
  raster.SAUP <- raster(paste("data/SAUP/raster/", species.list.s2[s], sep=""))
  
  png(paste("data/SAUP/maps/", species.list.s2[s], ".png", sep=""))
  plot(raster.SAUP)
  plot(SAUP, border="blue", add=TRUE)
  dev.off()
}


#600217

species.list.s <- species.list.s3
#species.list.s <- species.list[species.list$TaxonKey %in% as.numeric(species.saupmap),"TaxonKey"]
path <- "C:/data/modelling/data/SAUP/"
s=1

for(s in 1:length(species.list.s)){
  if(!file.exists(paste(path, "raster/", species.list.s[s], ".grd", sep=""))){
    poly <- readOGR(paste(path, "singleshapefiles", sep=""), species.list.s[s])
    
    poly.l <- rasterize(as(poly, 'SpatialLines'), data.extent, background=0)
    poly.p <- rasterize(poly, data.extent, background=0)
    spec.raster <- (poly.l + poly.p) > 0
    #spec.raster <- rasterize(poly, data.extent, background=0)
    spec.raster <- mask(spec.raster > 0, data.extent)
    writeRaster(spec.raster, paste("data/SAUP/raster/", species.list.s[s], sep=""), overwrite=T)
    
    png(paste("data/SAUP/maps/", species.list.s[s], ".png", sep=""))
    plot(spec.raster)
    plot(poly, border="blue", add=TRUE)
    dev.off()
    
  }else if(file.exists(paste(path, "raster/", species.list.s[s], ".grd", sep=""))){
    print(paste(species.list.s[s], "already converted"))
  }else{
    print("error")
  }
  rm(poly, spec.raster, poly.p, poly.l)
}


#Combine SAUP and IUCN raster

s=1
all.metrics=NULL
data.area=raster("C:/data/modelling/data/envi/AquaMaps_oarea.grd")
for(s in 1:length(species.list.s)){
  if(file.exists(paste("C:/data/modelling/data/SAUP/singleshapefiles/", species.list.s[s], ".shp", sep=""))){
    SAUP <- readOGR("data/SAUP/singleshapefiles", species.list.s[s])
  }
  if(file.exists(paste("C:/data/modelling/data/IUCN/singleshapefiles/", species.list.s[s], ".shp", sep=""))){
    IUCN <- readOGR("data/IUCN/singleshapefiles", species.list.s[s])
  }
  
  if(file.exists(paste("C:/data/modelling/data/SAUP/raster/", species.list.s[s], ".grd", sep=""))){
    raster.SAUP <- raster(paste("data/SAUP/raster/", species.list.s[s], sep=""))
  }
  if(file.exists(paste("C:/data/modelling/data/IUCN/raster/", species.list.s[s], ".grd", sep=""))){
    raster.IUCN <- raster(paste("data/IUCN/raster/", species.list.s[s], sep=""))
  }
  
  if(exists("raster.IUCN") & exists("raster.SAUP")){
    spec.raster <- (raster.SAUP > 0) + (raster.IUCN > 0)
    
    # Compare the two rasters
    rsSAUP <- sum(data.area[raster.SAUP > 0])
    rsIUCN <- sum(data.area[raster.IUCN > 0])
    rsBOTH <- sum(data.area[spec.raster > 0])
    
    # Calculate the intersection of the two rasters, this is given by adding 
    # the binary rasters together -> 2 indicates intersection
    combination <- (raster.SAUP > 0) + (raster.IUCN > 0)
    intersection <- combination == 2
    
    # Union is all the area covered by the both rasters
    union <- combination >= 1
    
    jaccard <- freq(intersection, value=1) / freq(union, value=1)
    metric <- c("TaxonKey"=species.list.s[s], "jaccard"=jaccard, "rsIUCN"=rsIUCN, "rsSAUP"=rsSAUP, "rsBOTH"=rsBOTH)
    all.metrics <- rbind(all.metrics, metric)
    rm(jaccard, combination, intersection, union)
    
  }else if(exists("raster.IUCN")){
    spec.raster <- raster.IUCN > 0
    rsIUCN <- sum(data.area[spec.raster > 0])
    metric <- c("TaxonKey"=species.list.s[s], "jaccard"=NA, "rsIUCN"=rsIUCN, "rsSAUP"=NA, "rsBOTH"=NA)
    all.metrics <- rbind(all.metrics, metric)
    
  }else if(exists("raster.SAUP")){
    spec.raster <- raster.SAUP > 0
    rsSAUP <- sum(data.area[spec.raster > 0])
    metric <- c("TaxonKey"=species.list.s[s], "jaccard"=NA, "rsIUCN"=NA, "rsSAUP"=rsSAUP, "rsBOTH"=NA)
    all.metrics <- rbind(all.metrics, metric)
    
  }else{
    print("error")
  }
  
  writeRaster(spec.raster, paste("C:/data/modelling/data/IUCN and SAUP/raster/", species.list.s[s], sep=""), overwrite=T)

  png(paste("C:/data/modelling/data/IUCN and SAUP/maps/", species.list.s[s], ".png", sep=""))
  plot(spec.raster)
  if(exists("raster.IUCN")){
    plot(IUCN, border="red", add=TRUE)}
  if(exists("raster.SAUP")){
    plot(SAUP, border="blue", add=TRUE)}
  dev.off()
  
  rm(SAUP, IUCN, raster.SAUP, raster.IUCN, spec.raster, metric, rsIUCN, rsSAUP, rsBOTH)
}
write.csv(all.metrics, file="C:/Users/LaurensG/Desktop/expert_rangemaps_metrics.csv")
  
#   
# 
# s=1
# for(s in 1:length(species.list.s)){
#   try(SAUP <- readOGR("data/SAUP/singleshapefiles", species.list.s[s]))
#   try(IUCN <- readOGR("data/IUCN/singleshapefiles", species.list.s[s]))
#   
#   try(spec.raster1 <- rasterize(SAUP, data.extent, background=0, getCover=T))
#   try(spec.raster1 <- mask(spec.raster1, data.extent))
#   try(spec.raster2 <- rasterize(IUCN, data.extent, background=0, getCover=T))
#   try(spec.raster2 <- mask(spec.raster2, data.extent))
#   
#   if(exists("spec.raster1") & exists("spec.raster2")){
#     spec.raster <- spec.raster1 + spec.raster2
#     
#   }else if(exists("spec.raster1")){
#     spec.raster <- spec.raster1
#     
#     png(paste("data/IUCN/maps/", species.list.s[s], ".png", sep=""))
#     plot(spec.raster1)
#     plot(SAUP, add=TRUE)
#     dev.off()
#     
#   }else if(exists("spec.raster2")){
#     spec.raster <- spec.raster2
#   }
#   spec.raster <- mask(spec.raster, data.extent)
#   
#   
#   png(paste("data/IUCN and SAUP/maps/", species.list.s[s], ".png", sep=""))
#   plot(spec.raster)
#   try(plot(IUCN, border="red", add=TRUE))
#   try(plot(SAUP, border="blue", add=TRUE))
#   dev.off()
#   
#   #plot(spec.raster1)
#   #plot(spec.raster2)
#   #plot(spec.raster > 0)
#   
#   try(writeRaster(spec.raster, paste("data/IUCN and SAUP/raster/", species.list.s[s], sep=""), overwrite=T))
#   try(writeRaster(spec.raster1, paste("data/SAUP/raster/", species.list.s[s], sep=""), overwrite=T))
#   try(writeRaster(spec.raster2, paste("data/IUCN/raster/", species.list.s[s], sep=""), overwrite=T))
# 
#   rm(SAUP, IUCN, spec.raster, spec.raster1, spec.raster2)
# }

