
# install.packages("XML", dependencies = TRUE)
# download.file("http://www.omegahat.org/Prerelease/XMLSchema_0.8-0.tar.gz", "XMLSchema")
# install.packages("XMLSchema", type="source", repos = NULL)
# download.file("http://www.omegahat.org/Prerelease/SSOAP_0.91-0.tar.gz", "SSOAP")
# install.packages("SSOAP", type="source", repos = NULL)
# library(SSOAP)


#*#*#*#*#*#*#*#*#
# FUNCTIONS ----

  #Search GBIF for occurrence records
  get.gbif.occs <- function(p.species, fields){
    #Get GBIF backbone ID of the species by searching for name
    gbif.id <- name_backbone(name=p.species, kingdom="Animalia")$speciesKey
    #Download occurrence records for the species
    gbif.data <- occ_search(taxonKey=gbif.id, hasCoordinate=T, spatialIssues=F, return="data", fields=fields, limit=1000000)
    return(gbif.data)}

  #Sample background and create background points with probability correcting for cell area
  sample.background <- function(data.extent, p.number.of.background.points){
    data.background.points <- data.frame(randomPoints(data.extent, n=p.number.of.background.points, lonlatCorrection=T))
    data.background.points <- data.frame(TaxonKey=0, TaxonName="none", x=data.background.points$x, y=data.background.points$y, set="background", occ=0, stringsAsFactors=F)
    return(data.background.points)}


#*#*#*#*#*#*#
# HEAD ----

  # load libraries
    libraries = c("XML", "XMLSchema", "SSOAP", "rgbif", "plyr", "raster", "dismo")
    #try(install.packages(libraries, configure.args=c(repos="http://cran.us.r-project.org")))
    try(lapply(libraries, require, character.only=T))
    rm(libraries)

  # set working directory
    setwd("C:/data/modelling/data/species/")

  

#*#*#*#*#*#*#
# BODY ----

# read data species list
  species.list <- read.csv("C:/data/modelling/data/species/speciesinfo.new.csv", header=T)

# #~Get WoRMS data ----
# # connect to WoRMS
#   webconnect = processWSDL("http://www.marinespecies.org/aphia.php?p=soap&wsdl=1")
#   WoRMS = genSOAPClientInterface(, webconnect)
# 
# # get original species id's
#   results=NULL
#   for(s in 1:length(species.list[,1])){
#     # find Aphia ID of record
#     #TaxonName = species.list[s,2]
#     #AphiaID <- WoRMS@functions$getAphiaID(TaxonName,marine_only=0,('http://www.marinespecies.org/aphia.php?p=soap'))
#     AphiaID <- species.list$AphiaID[s]
#     AphiaInfo <- WoRMS@functions$getAphiaRecordByID(AphiaID,('http://www.marinespecies.org/aphia.php?p=soap'))
#     result <- as.data.frame(rbind(unlist(lapply(slotNames(AphiaInfo), function(x) slot(AphiaInfo, x)))))
#     result <- as.data.frame(c(species.list[s,1:2], result))
#     results <- rbind.fill(results, result)
#     print(c(species.list[s,1], AphiaID, AphiaInfo@match_type))
#   }
#   colnames(results) <- c("TaxonKey", "TaxonName", slotNames(AphiaInfo)[c(1:12,14:25)])
#   rm(result, AphiaID, AphiaInfo, WoRMS, webconnect, s, "http://aphia/v1.0")
#   data.worms <- results
#   write.csv(results, file="speciesinfo.csv")
#   save(data.worms, file="data.worms.RData")
# 
#   data.worms <- read.csv("speciesinfo.csv")
#   species.list <- data.frame("TaxonKey"=data.worms$TaxonKey, "TaxonName"=data.worms$valid_name, "oldTaxonName"=data.worms$TaxonName, data.worms[,4:28])
#   write.csv(species.list, file="C:/data/modelling/data/species/specieslistworms.csv", row.names=F)
# 
# species.list <- read.csv(file="C:/data/modelling/data/species/specieslistworms.csv", header=T)
# 
# 
# #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#
# #~Get GBIF occurence records based on species name ----
#   s=1
#   data.gbif=NULL
#   options(warn=1)
#   fields=c("species","decimalLongitude","decimalLatitude",
#            "scientificName","taxonRank","occurrenceStatus","occurrenceRemarks",
#            "year","month","depth","depthAccuracy","protocol","samplingProtocol","issues","coordinateAccuracy","geodeticDatum",
#            "georeferenceSources","georeferenceProtocol","georeferenceVerificationStatus","georeferenceRemarks")
# 
# #not found: 989, 1024, 1069, 1077
# #skipping 370 as it takes ages(?)
#   for(s in 307){#c(307, 338, 370, 627, 1032)){#1078:length(data.worms$valid_name)){
#     raw.gbif <- get.gbif.occs(as.character(data.worms$valid_name[s]), fields)
#     raw.gbif <- cbind(species.list[s,1:2], raw.gbif, row.names=NULL)
#     save(raw.gbif, file=paste("GBIF files/", species.list[s,1], ".RData", sep=""))
#     data.gbif <- rbind.fill(data.gbif, raw.gbif)
#     print(c(species.list[s,1], "done"))
#   }
#   colnames(data.gbif) <- c(colnames(data.gbif)[1:6], "longitude", "latitude", colnames(data.gbif)[9:22])
#   save(data.gbif, file="data.gbif.RData", overwrite=T)
# 

#~Get OBIS occurence records based on species name ----
  s=1
  data.obis=NULL
  notfound=NULL
  options(warn=1)
  
  for(s in 235:length(species.list[,1])){
    filepath <- paste("C:/data/modelling/data/species/OBIS files/", species.list[s,1],".csv", sep="")
    if (file.exists(filepath)){
      raw.obis <- read.csv(filepath, header=T, sep=",", quote='"', dec=".", encoding="latin1", stringsAsFactor=F)
      raw.obis <- cbind(species.list[s,1:2],raw.obis, row.names=NULL)
      data.obis <- rbind.fill(data.obis, raw.obis)
      print(paste(species.list[s,1]," done"))
      rm(raw.obis)}
    else{
      print(paste(species.list[s,1],"not found"))
      notfound <- c(notfound, species.list[s,1])
    }
  }
  #not found: 600427 622620 690077 690082 690143 690176 690230 690288 690411 690445 690614 690618 690624
  save(data.obis, file="data.obis.RData", overwrite=T)

#*#*#*#*#*#*#*#*#*#*#*#*#
#~Process species data ----
  load("C:/data/modelling/data/species/data.obis.RData")
  summary(data.obis)
  
  load("C:/data/modelling/data/species/data.gbif.RData")
  summary(data.gbif)
  
  #how many names in data per name in species records?
  gbifnames <- unique(data.gbif[,c(1,2,5)])
  sort(table(gbifnames[,2]))
  obisnames <- unique(data.obis[,c(1,2,7)])
  sort(table(obisnames[,2]))
  
  #Join data tables
  data.gbif <- cbind(data.gbif, "source"="gbif")
  data.obis <- cbind(data.obis, "source"="obis")
  data.species <- rbind(data.obis[,c("TaxonKey", "TaxonName", "latitude", "longitude", "source")], data.gbif[,c("TaxonKey", "TaxonName", "latitude", "longitude", "source")])
  save(data.species, file="C:/data/modelling/data/species/1_data.species.initial.RData")

#~Clean data
  load("C:/data/modelling/data/species/1_data.species.initial.RData")
  #Remove rows with NA [9,165,440 no change](previously 10,029,039)
  data.species <- data.species[complete.cases(data.species),]
  #Remove rows without species.name [no change]
  data.species <- data.species[data.species$TaxonName != "",]
  #Remove rows with invalid coordinates [no change]
  data.species <- data.species[!data.species$latitude > 90 & !data.species$latitude< -90,]
  data.species <- data.species[!data.species$longitude > 180 & !data.species$longitude< -180,]
  #Get unique species records 9,165,440 -> 2,809,749]
  data.species <- unique(data.species[,1:4])
  #Remove any records that have 0 degree longitude OR 0 degree latitude! [3,178,513 -> 3,177,577]
  #data.species <- data.species[!data.species$latitude==0 & !data.species$longitude==0,]
  #drop unused levels (species names) from factor data
  #write.csv(data.species, file="C:/data/modelling/data/species/data.species.cleaned.csv", row.names=F)
  #data.species <- read.csv("C:/data/modelling/data/species/data.species.cleaned.csv", header=T)


#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#
#~Aggregate species data to grid cells ----
  data.extent <- raster("C:/data/modelling/data/envi/extent05deg.grd")
  data.species.aggregated=NULL
  s=1
  for(s in 1:length(species.list$TaxonKey)){
    l.specnum = as.character(species.list$TaxonKey[s])
    l.specname = as.character(species.list$TaxonName[s])
    print(c(l.specnum, l.specname))
    
    l.species <- data.species[data.species$TaxonKey==l.specnum,]
    records.raster.species <- rasterize(cbind(l.species$longitude, l.species$latitude), data.extent, field=1, fun="sum")
    #plot(records.raster.species, xlim=c(-180,180), ylim=c(-90,90), col="red")
    #points(SpatialPoints(l.species[,4:3]), xlim=c(-180,180), ylim=c(-90,90))
    points=rasterToPoints(records.raster.species, spatial=F)
    
    if(dim(points)[1]==1){
      result <- data.frame(TaxonKey=l.specnum, TaxonName=l.specname, x=points[1], y=points[2], occ=1, set="none", stringsAsFactors=F)
    }else{
      result <- data.frame(TaxonKey=l.specnum, TaxonName=l.specname, points[,1:2], occ=1, set="none", stringsAsFactors=F)
    }
    result$set[sample(x=length(result[,1]), size=ceiling(length(result[,1])/2), replace=F)] <- "training"
    result$set[result$set=="none"] <- "test"
    data.species.aggregated <- rbind(data.species.aggregated, result)
  }
  #tidy up
  try(rm(s, l.specnum, l.specname, l.species, records.raster.species, points, result))
  save(data.species.aggregated, file="C:/data/modelling/data/species/2_data.species.aggregated.RData")


# #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#
# # FAO AREAS ----
# 
# #~Get FAO areas from FishBase ----
#   #install.packages("devtools")
#   library(devtools)
#   #devtools::install_github("ropensci/rfishbase@rfishbase2.0")
#   devtools::install_github("ropensci/rfishbase")
#   library("rfishbase")
#   
#   s=1
#   fao=NULL
#   notfound=NULL
#   #Find commercial species on fishbase
#   for(s in 849:length(species.list$TaxonName)){
#     try(s.fao <- faoareas(species_list=as.vector(species.list$TaxonName[s])))
#     if(exists("s.fao")){
#       print(species.list$TaxonKey[s])
#       s.fao <- cbind(TaxonKey=as.numeric(species.list$TaxonKey[s]), TaxonName=as.character(species.list$TaxonName[s]), s.fao)
#       fao <- c(fao, list(s.fao))
#     }else{
#       print(c(species.list$TaxonKey[s], "not found"))
#       notfound <- c(notfound, species.list$TaxonKey[s])
#       s.fao <- cbind(TaxonKey=as.numeric(species.list$TaxonKey[s]), TaxonName=as.character(species.list$TaxonName[s]), Error="not found")
#       fao <- c(fao, list(s.fao))
#     }
#     rm(s.fao)
#   }
#   data.FAOareas <- fao
#   #not found: 600238, 600664, 601455, 602420, 604791, 606422, 607046, 690003, 690005, 690009, 690010, 690011, 690013, 690016, 690018, 690022, 690023, 690024, 690029, 690030, 690031, 690032, 690033, 690036, 690039, 690042, 690043, 690047, 690049, 690050, 690051, 690052, 690053, 690054, 690055, 690058, 690059, 690060, 690064, 690065, 690067, 690069, 690070, 690071, 690075, 690077, 690078, 690081, 690082, 690085, 690087, 690088, 690089, 690090, 690091, 690092, 690100, 690105, 690109, 690115, 690119, 690120, 690122, 690123, 690124, 690125, 690126, 690127, 690133, 690143, 690146, 690151, 690156, 690157, 690158, 690161, 690163, 690166, 690167, 690168, 690169, 690176, 690180, 690182, 690185, 690189, 690190, 690194, 690196, 690198, 690199, 690201, 690203, 690205, 690206, 690215, 690217, 690218, 690219, 690227, 690228, 690230, 690232, 690241, 690242, 690249, 690253, 690257, 690259, 690265, 690268, 690269, 690270, 690271, 690272, 690273, 690274, 690275, 690279, 690282, 690283, 690284, 690285, 690286, 690287, 690288, 690296, 690302, 690303, 690305, 690309, 690313, 690314, 690318, 690320, 690321, 690322, 690325, 690327, 690328, 690332, 690334, 690335, 690336, 690337, 690338, 690341, 690343, 690354, 690356, 690364, 690370, 690374, 690377, 690378, 690379, 690381, 690383, 690387, 690392, 690398, 690399, 690400, 690403, 690404, 690406, 690411, 690420, 690423, 690425, 690429, 690430, 690432, 690433, 690434, 690440, 690444, 690445, 690455, 690470, 690595, 690598, 690599, 690606, 690608, 690609, 690610, 690611, 690612, 690613, 690614, 690615, 690616, 690617, 690618, 690619, 690620, 690621, 690622, 690624, 690625, 690628, 690629, 690631, 690649, 690651, 690658, 690665, 690674, 690675, 690676, 690677, 690678, 690679, 690680, 690681, 690682, 690683, 690684, 690685, 690686, 690687, 690688, 690689, 690690, 690691, 690692, 690693, 690694, 690750, 690751, 999004, 999005, 999021, 999022
#   save(data.faoareas, file="C:/data/modelling/data/species/data.faoareas.RData")

#~Clip specdata by FAO area ----
  load("C:/data/modelling/data/species/data.faoareas.RData")
  data.faoraster <- raster("C:/data/modelling/data/envi/faoraster.grd")
  load("C:/data/modelling/data/species/2_data.species.aggregated.RData")
  #species.list[species.list$TaxonName=="Peprilus simillimus",1:5]

  s=1
  data.species.clipped=NULL
  ##dataframe of fao distributions
  #species.list.faoareas <- matrix(data=NA, nrow=0, ncol=22)
  #colnames(species.list.faoareas) <- c("TaxonKey", "TaxonName", sort(unique(values(data.faoraster))))
  #species.list.faoareas <- as.data.frame(species.list.faoareas)
  #read all species data and clip it by FAO areas if data is available
  for(s in 1:length(species.list$TaxonKey)){
    if(data.faoareas[data.faoareas$TaxonKey==species.list$TaxonKey[s], c("NASAUP")]=="x" & data.faoareas[data.faoareas$TaxonKey==species.list$TaxonKey[s], c("NAFISHBASE")]=="x"){
      print(c("no fao data for", species.list$TaxonKey[s], as.character(species.list$TaxonName[s])))
      s.data.faoareas <- as.integer(c(18, 27, 21, 61, 67, 37, 77, 31, 34, 51, 57, 71, 87, 41, 47, 81, 58, 48, 88))
      s.data.species <- data.species.aggregated[data.species.aggregated$TaxonKey==species.list$TaxonKey[s],]
      s.selectedfao <- extract(data.faoraster, s.data.species[,c("x","y")])
      result <- s.data.species[s.selectedfao %in% s.data.faoareas,]
      data.species.clipped <- rbind(data.species.clipped, result)
      
      png(paste("C:/data/modelling/data/species/clip/", species.list$TaxonKey[s] ,".png", sep=""), width=800, height=500, units="px")
      plot(data.faoraster, main=paste(species.list$TaxonKey[s], species.list$TaxonName[s]))
      points(s.data.species[,c("x","y")], pch=1)
      points(s.data.species[,c("x","y")], pch=20, col="red")
      dev.off()
    }else{
      s.data.faoareas <- unique(as.numeric(data.faoareas[data.faoareas$TaxonKey==species.list$TaxonKey[s],3:28]))
      s.data.species <- data.species.aggregated[data.species.aggregated$TaxonKey==species.list$TaxonKey[s],]
      s.selectedfao <- extract(data.faoraster, s.data.species[,c("x","y")])
      result <- s.data.species[s.selectedfao %in% s.data.faoareas,]
      data.species.clipped <- rbind(data.species.clipped, result)
      
      png(paste("C:/data/modelling/data/species/clip/", species.list$TaxonKey[s] ,".png", sep=""), width=800, height=500, units="px")
      plot(data.faoraster, main=paste(species.list$TaxonKey[s], species.list$TaxonName[s]))
      points(s.data.species[,c("x","y")], pch=1)
      points(result[,c("x","y")], pch=20, col="red")
      dev.off()
      
      print(c(species.list$TaxonKey[s], as.character(species.list$TaxonName[s])))
      
    }
    #species.list.faoareas <- rbind.fill(species.list.faoareas, cbind(TaxonKey=species.list$TaxonKey[s], TaxonName=as.character(species.list$TaxonName[s]), as.data.frame(t(as.data.frame(s.data.faoareas, row.names=s.data.faoareas)))))
    rm(s, s.data.faoareas, s.selectedfao)
  }
  save(data.species.clipped, file="C:/data/modelling/data/species/3_data.species.clipped.RData")
  #write.csv(species.list.faoareas, file="C:/data/modelling/data/species/species.list.faoareas.csv")


#~Add background data
  load("C:/data/modelling/data/species/3_data.species.clipped.RData")
  a <- sample.background(data.extent, 10000)
  data.species.final <- rbind(data.species.clipped, a)
  save(data.species.final, file="C:/data/modelling/data/species/4_data.species.final.RData")

  load("C:/data/modelling/data/species/4_data.species.final.RData")
  data.species <- data.species.final
  save(data.species, file="C:/data/modelling/data/species/data.species.RData")



#*#*#*#*#*#*#*#*#*#*#*#*#*#
#~Modify species list ----

  load("C:/data/modelling/data/species/1_data.species.initial.RData")
  load("C:/data/modelling/data/species/4_data.species.final.RData")
  species.list <- read.csv("C:/data/modelling/data/species/speciesinfo.csv", header=T)

  #numer of records per species
  rownames(species.list) <- species.list$TaxonKey
  species.list <- merge(species.list, cbind(table(data.species.final$TaxonKey)), by="row.names", all.x=TRUE)
  #number of cells with species recorded
  rownames(species.list) <- species.list$TaxonKey
  species.list <- merge(species.list, cbind(table(data.species$TaxonKey)), by="row.names", all.x=TRUE)
  #drop "row.names" column
  species.list <- species.list[3:length(species.list)]
  #correct column names for record and cell count
  colnames(species.list) <- c(colnames(species.list)[3:length(colnames(species.list))-2], "ncells", "nrecords")
  
  #assign core to each species
  species.list <- cbind(species.list, "Core"=sample(x=1:16, size=1086, replace=T))
  #rearrange columns
  species.list <- species.list[,c("TaxonKey", "TaxonName", "Core", "ncells", "nrecords", "oldTaxonName", "AphiaID", "url", "scientificname", "authority", "rank", "status", "unacceptreason", "valid_AphiaID", "valid_name", "valid_authority", "comment", "kingdom", "phylum", "order", "family", "genus", "citation", "lsid", "isMarine", "isBrackish", "isFreshwater", "isTerrestrial", "isExtinct", "match_type", "modified")]
  
  write.csv(species.list, file="C:/data/modelling/data/species/speciesinfo.new.csv", row.names=F)
  species.list <- read.csv("C:/data/modelling/data/species/speciesinfo.new.csv", header=T)



