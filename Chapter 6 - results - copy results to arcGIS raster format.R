
for(s in c(600143, 601342, 605367, 600656, 600138)){ 
  
  if(file.exists(paste("C:/data/modelling/data/SAUP/raster/", s, ".grd", sep=""))){
    SAUP <- raster(paste("C:/data/modelling/data/SAUP/raster/", s, ".grd", sep=""))
    writeRaster(SAUP, file=paste("C:/data/modelling/output/arcGIS/example maps occ exp base geod/", s, " SAUP expert range.asc", sep=""))
  }
  if(file.exists(paste("C:/data/modelling/data/IUCN/raster/", s, ".grd", sep=""))){
    IUCN <- raster(paste("C:/data/modelling/data/IUCN/raster/", s, ".grd", sep=""))
    writeRaster(IUCN, file=paste("C:/data/modelling/output/arcGIS/example maps occ exp base geod/", s, " IUCN expert range.asc", sep=""))
  }
  
  write.csv(data.species.aggregated[data.species.aggregated$TaxonKey==s,], file=paste("C:/data/modelling/output/arcGIS/example maps occ exp base geod/", s, " unclipped.csv", sep=""))
  
  load(paste("C:/data/modelling/output/baseline/SDMeval_fulld_", s, ".RData", sep=""))
  writeRaster(results[[2]][[2]], file=paste("C:/data/modelling/output/arcGIS/example maps occ exp base geod/", s, " maxent baseline baseline.asc", sep=""))
    
  load(paste("C:/data/modelling/output/baseline/SDMeval.clipped_", s, ".RData", sep=""))
  writeRaster(results[[2]][[2]], file=paste("C:/data/modelling/output/arcGIS/example maps occ exp base geod/", s, " maxent baseline clipped.asc", sep=""))
  
  load(paste("C:/data/modelling/output/geodist2/SDMeval_geodist_", s, ".RData", sep=""))
  writeRaster(results[[2]][[2]], file=paste("C:/data/modelling/output/arcGIS/example maps occ exp base geod/", s, " maxent baseline geodist.asc", sep=""))
}

