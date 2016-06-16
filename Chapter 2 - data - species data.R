
# install.packages("XML", dependencies = TRUE)
# download.file("http://www.omegahat.org/Prerelease/XMLSchema_0.8-0.tar.gz", "XMLSchema")
# install.packages("XMLSchema", type="source", repos = NULL)
# download.file("http://www.omegahat.org/Prerelease/SSOAP_0.91-0.tar.gz", "SSOAP")
# install.packages("SSOAP", type="source", repos = NULL)
# library(SSOAP)

# FUNCTIONS ----
  get.gbif.occs <- function(p.species, fields){
    #Get GBIF backbone ID of the species by searching for name
    gbif.id <- name_backbone(name=p.species, kingdom="Animalia")$speciesKey
    #Download occurrence records for the species
    gbif.data <- occ_search(taxonKey=gbif.id, hasCoordinate=T, spatialIssues=F, return="data", fields=fields, limit=200000)
    return(gbif.data)}

# HEAD ----
# load libraries
  libraries = c("XML", "XMLSchema", "SSOAP", "rgbif", "plyr")
  #try(install.packages(libraries, configure.args=c(repos="http://cran.us.r-project.org")))
  try(lapply(libraries, require, character.only=T))
  rm(libraries)
# set working directory
  setwd("C:/data/modelling/data/species/")
# read data species list
  species.list <- read.csv("C:/data/modelling/data/species/speciesinfo.new.csv", header=T)
  

# BODY ----

# connect to WoRMS
  webconnect = processWSDL("http://www.marinespecies.org/aphia.php?p=soap&wsdl=1")
  WoRMS = genSOAPClientInterface(, webconnect)

# get original species id's
  results=NULL
  for(s in 1:length(species.list[,1])){
    # find Aphia ID of record
    #TaxonName = species.list[s,2]
    #AphiaID <- WoRMS@functions$getAphiaID(TaxonName,marine_only=0,('http://www.marinespecies.org/aphia.php?p=soap'))
    AphiaID <- species.list$AphiaID[s]
    AphiaInfo <- WoRMS@functions$getAphiaRecordByID(AphiaID,('http://www.marinespecies.org/aphia.php?p=soap'))
    result <- as.data.frame(rbind(unlist(lapply(slotNames(AphiaInfo), function(x) slot(AphiaInfo, x)))))
    result <- as.data.frame(c(species.list[s,1:2], result))
    results <- rbind.fill(results, result)
    print(c(species.list[s,1], AphiaID, AphiaInfo@match_type))
  }
  colnames(results) <- c("TaxonKey", "TaxonName", slotNames(AphiaInfo)[c(1:12,14:25)])
  rm(result, AphiaID, AphiaInfo, WoRMS, webconnect, s, "http://aphia/v1.0")
  data.worms <- results
  write.csv(results, file="speciesinfo.csv")
  save(data.worms, file="data.worms.RData")

  load("data.worms.RData")
  data.worms <- read.csv("speciesinfo.csv")


# Get GBIF occurence records based on species name ----
  s=1
  data.gbif=NULL
  options(warn=1)
  fields=c("species","decimalLongitude","decimalLatitude",
           "scientificName","taxonRank","occurrenceStatus","occurrenceRemarks",
           "year","month","depth","depthAccuracy","protocol","samplingProtocol","issues","coordinateAccuracy","geodeticDatum",
           "georeferenceSources","georeferenceProtocol","georeferenceVerificationStatus","georeferenceRemarks")

#skipping 370 as it takes ages(?)
  for(s in 371:length(data.worms$valid_name)){
    raw.gbif <- get.gbif.occs(as.character(data.worms$valid_name[s]), fields)
    raw.gbif <- cbind(species.list[s,1:2], raw.gbif, row.names=NULL)
    save(raw.gbif, file=paste("GBIF files/", species.list[s,1], ".RData", sep=""))
    data.gbif <- rbind.fill(data.gbif, raw.gbif)
    print(c(species.list[s,1], "done"))
  }
  save(data.gbif, file="GBIF.RData", overwrite=T)


# Get OBIS occurence records based on species name ----
  s=1
  data.obis=NULL
  options(warn=1)
  
  for(s in 1:length(species.list[,1])){
    filepath <- paste("C:/data/modelling/data/species/OBIS files/", species.list[s,1],".csv", sep="")
    if (file.exists(filepath)){
      raw.obis <- read.csv(filepath, header=T, sep=",", quote='"', dec=".", encoding="latin1", stringsAsFactor=F)
      raw.obis <- cbind(species.list[s,1:2],raw.obis, row.names=NULL)
      data.obis <- rbind.fill(data.obis, raw.obis)
      print(paste(species.list[s,1]," done"))
      rm(raw.obis)}
    else{
      print(paste(species.list[s,1],"not found"))}
  }
  #not found: 600427
  save(data.obis, file="OBIS.RData", overwrite=T)


##############
##############

load("data.obis.RData")
load("data.gbif.RData")

# Join data tables
  data.species <- rbind(data.obis[,c("TaxonKey", "TaxonName", "latitude", "longitude")], data.gbif[,c("TaxonKey", "TaxonName", "latitude", "longitude")])

# ~ clean data ----
    #Remove rows with NA [10,624,626 no change]
    data.species <- data.species[complete.cases(data.species),]
    #Remove rows without species.name [no change]
    data.species <- data.species[data.species$TaxonName != "",]
    #Remove rows with invalid coordinates [no change]
    data.species <- data.species[!data.species$latitude > 90 & !data.species$latitude< -90,]
    data.species <- data.species[!data.species$longitude > 180 & !data.species$longitude< -180,]
    
    #Get unique species records [10,624,626 -> 3,395,324]
    data.species <- unique(data.species)
    
    #Remove any records that have 0 degree longitude OR 0 degree latitude! [3,395,324 -> 3,394.397]
    data.species <- data.species[!data.species$latitude==0 & !data.species$longitude==0,]

  save(data.species, file="C:/data/modelling/data/specdata/data.species.RData")
  load("C:/data/modelling/data/species/data.species.RData")


# in species list but not in species data: 63 - Auxis rochei rochei, Auxis thazard thazard, Sarda chiliensis chiliensis, Salmo trutta trutta, Oncorhynchus masou masou, Salvelinus alpinus alpinus, Osmerus mordax mordax, Merluccius gayi gayi, Emmelichthys nitidus nitidus, Centroscymnus cryptacanthus, Mullus barbatus barbatus, Tylosurus crocodilus crocodilus, Scomberesox saurus saurus, Etrumeus teres, Clupea pallasii pallasii, Diplodus sargus sargus, Gasterosteus aculeatus aculeatus, Diplodus cervinus cervinus, Diplodus argenteus argenteus, Pseudopentaceros richardsoni, Trachyscorpia cristulata cristulata, Hoplostethus mediterraneus mediterraneus, Dipturus linteus, Chrysophrys auratus, Acanthopagrus schlegelii schlegelii, Lepidonotothen mizops, Lepidonotothen nudifrons, Lepidonotothen larseni, Patagonotothen brevicauda brevicauda, Penaeus merguiensis, Penaeus stylirostris, Penaeus kerathurus, Mytilus edulis platensis, Aulacomya atra, Cerastoderma edule, Penaeus brevirostris, Penaeus chinensis, Microcosmus vulgaris, Holthuispenaeopsis atlantica, Anadara sativa, Penaeus indicus, Haliporoides sibogae sibogae, Haliporoides triarthrus triarthrus, Penaeus japonicus, Penaeus aztecus, Penaeus duorarum, Penaeus setiferus, Venerupis corrugata, Penaeus brasiliensis, Penaeus paulensis, Aristaeopsis edwardsiana, Lithodes santolla, Penaeus notialis, Palinurus gilchristi, Lithodes maja, Nototodarus sloani, Penaeus latisulcatus, Penaeus vannamei, Penaeus californiensis, Mizuhopecten yessoensis, Lithodes aequispinus, Zygochlamys delicatula, Charonia lampas rubicunda
unique(data.species$TaxonName)[!unique(data.species$TaxonName) %in% species.list$TaxonName]
# Zero records: 1 - Pecten yessoensis (0)
species.list[!species.list$TaxonName %in% unique(data.species$TaxonName),]
# Less than 50 records: 179 - Carangoides bajad (49), Palinurus mauritanicus (49), Somniosus rostratus (49), Abudefduf luridus (48), Marsupenaeus japonicus (47), Paralonchurus peruanus (47), Penaeus chinensis (47), Sardinella fimbriata (46), Sarda chiliensis (45), Bathyraja irrasa (44), Gymnothorax unicolor (44), Oncorhynchus masou masou (44), Amblyraja taaf (42), Larimichthys crocea (42), Raja brachyura (42), Hyperoglyphe bythites (41), Loxechinus albus (41), Metapenaeus joyneri (41), Chionobathyscus dewitti (40), Bythaelurus canescens (39), Glaucosoma hebraicum (39), Penaeus brevirostris (39), Acanthurus sohal (38), Epinephelides armatus (38), Meretrix lusoria (38), Mytilus planulatus (38), Rhabdosargus haffara (38), Penaeus vannamei (38), Chaceon bicolor (37), Coregonus oxyrinchus (37), Mytilus coruscus (37), Pleuroncodes planipes (37), Ibacus ciliatus (36), Scolopsis taeniata (36), Turbo cornutus (36), Penaeus merguiensis (36), Aristeus varidens (35), Sicyonia ingentis (35), Sarda chiliensis chiliensis (35), Merluccius gayi gayi (35), Rhinobatos planiceps (34), Concholepas concholepas (33), Joturus pichardi (33), Laemonema longipes (33), Rastrelliger brachysoma (33), Sardinella brasiliensis (33), Diplodus argenteus argenteus (32), Panulirus homarus (31), Prolatilus jugularis (31), Sardinella longiceps (31), Zidona dufresnei (31), Chirocentrus nudus (30), Hyporhamphus sajori (30), Liza saliens (30), Oncorhynchus masou (30), Sagmariasus verreauxi (30), Acanthopagrus schlegelii (29), Fenneropenaeus indicus (28), Glossanodon semifasciatus (28), Miichthys miiuy (28), Doryteuthis gahi (27), Genypterus maculatus (27), Isacia conceptionis (27), Larimichthys polyactis (27), Paralomis granulosa (27), Pogonophryne permitini (27), Scomberomorus niphonius (27), Solenocera agassizii (27), Diplodus cervinus (26), Lampris immaculatus (26), Nematopalaemon schmitti (26), Ostrea lurida (26), Tenualosa toli (26), Mizuhopecten yessoensis (26), Parapristipoma octolineatum (25), Psettodes bennettii (25), Aphanopus intermedius (24), Haliporoides diomedeae (24), Arctoscopus japonicus (23), Jasus frontalis (23), Paralomis formosa (23), Paralomis spinosissima (23), Acipenser medirostris (22), Acipenser ruthenus (22), Metanephrops andamanicus (22), Panulirus longipes (22), Paralithodes brevipes (22), Pseudopleuronectes herzensteini (22), Penaeus stylirostris (22), Haliotis gigantea (21), Lithodes santolla (21), Takifugu vermicularis (20), Dipturus linteus (20), Patagonotothen brevicauda brevicauda" (20), Artemia salina (19), Centroscymnus owstonii (19), Ethmidium maculatum (18), Farfantepenaeus brevirostris (18), Litopenaeus vannamei (18), Penaeus penicillatus (18), Trachinotus mookalee (18), Acipenser gueldenstaedtii (17), Plectropomus pessuliferus (17), Farfantepenaeus californiensis (16), Huso huso (16), Normanichthys crockeri (16), Scyllarides latus (16), Centroscymnus cryptacanthus (16), Acipenser stellatus (15), Bathyraja meridionalis (15), Clupanodon thrissa (15), Cynoscion analis (15), Polititapes rhomboides (15), Scomberomorus lineolatus (15), Haliporoides triarthrus (14), Penaeus indicus (14), Grammoplites suppositus (13), Lithodes antarcticus (13), Spisula sachalinensis (13), Hyporhamphus ihi (12), Semele solida (12), Brevoortia aurea (11), Palinurus delagoae (11), Acetes japonicus (10), Cheilopogon agoo (10), Choerodon rubescens (10), Clupeonella cultriventris (10), Genypterus chilensis (10), Heterocarpus reedi (10), Panulirus regius (10), Rajella lintea (10), Tawera elliptica (10), Nototodarus sloani (10), Alepocephalus umbriceps (9), Cervimunida johni (9), Heterocarpus vicarius (9), Mithrax armatus (9), Nototodarus sloanii (9), Patagonotothen brevicauda (9), Pollicipes pollicipes (9), Scarus persicus (9), Solen vagina (9), Cheilodactylus variegatus (8), Chione stutchburyi (8), Cilus gilberti (8), Fenneropenaeus merguiensis (8), Gerres nigri (8), Pleurogrammus azonus (8), Spisula ovalis (8), Haliporoides triarthrus triarthrus (8), Apostichopus japonicus (7), Choromytilus chorus (7), Hypoptychus dybowskii (7), Jasus paulensis (6), Mesodesma donacium (6), Panulirus cygnus (6), Microcosmus vulgaris (6), Metanephrops mozambicus (5), Eleginus nawaga (4), Pyura chilensis (4), Brevoortia pectinata (3), Clupea bentincki (3), Alosa immaculata (2), Anadara inaequivalvis (2), Litopenaeus stylirostris (2), Nibea mitsukurii (2), Paralomis aculeata (2), Parapenaeopsis atlantica (2), Tapes pullastra (2), Anadara sativa (2), Charonia lampas (1), Farfantepenaeus paulensis (1), Leukoma thaca (1), Liza klunzingeri (1), Merluccius gayi (1), Mytilus chilensis (1), Spaniblennius riodourensis (1), Holthuispenaeopsis atlantica (1)
sort(table(data.species$TaxonName)[table(data.species$TaxonName)<50], decreasing=T)




#########
a <- sample.background(data.envi[[1]], 10000)[,1:2]
for(s in 1:length(specieslist[,1])){
  #get training data
  p <- data.species[data.species$TaxonKey==specieslist[s,1] & data.species$set=="training",3:4]
  #get test data
  e <- data.species[data.species$TaxonKey==specieslist[s,1] & data.species$set=="test",3:4]
  