#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#
#R code modelling commercial fish species distributions #
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#

#PhD project
#Jan Laurens Geffert
#Department of Geography
#University of Cambridge
#jlg57@cam.ac.uk


#*#*#*#*#*#*#*#*#
## LIBRARIES ####
#*#*#*#*#*#*#*#*#

x = c("sp", "raster", "dismo", "rJava", "gbm", "plyr", "ENMeval", "vegan", "MASS", "psych", "rasterVis")
#try(install.packages(x, configure.args=c(repos="http://cran.us.r-project.org")))
try(lapply(x, require, character.only=T))
rm(x)


#*#*#*#*#*#*#*#*#
## FUNCTIONS ####
#*#*#*#*#*#*#*#*#

# Create output directory
f.create.output.dir <- function(output.dir){
  if(file.exists(output.dir)==F){
    dir.create(output.dir, recursive=T)
    x <- ("#### output directory created")
  }else if(file.exists(output.dir)==T){
    x <- ("#### output directory already exists")
  }else{
    x <- ("#### error creating output directory!")}
  return(x)}

# Function to run maxent model with settings
f.MAXENTrun <- function(envidata, p, a, p.output.dir){
  m <- maxent(x=envidata, p=p, a=a,
              removeDuplicates=T, overwrite=T, path=paste(p.output.dir, "/maxent", sep=""),
              #See Maxent.jar programme -> help for a list and summary of these commands
              args=c("removeduplicates=true", "linear=TRUE", "quadratic=TRUE", "product=TRUE", "threshold=TRUE", "hinge=TRUE",
                     "responsecurves=true", "pictures=false", "plots=true", "writeplotdata=true",
                     "askoverwrite=false", "skipifexists=false", "randomseed=false", "randomtestpoints=0", "jackknife=false",
                     "writeclampgrid=false", "writemess=false", "warnings=true", "allowpartialdata=false",
                     "betamultiplier=1.0"))    
  return(m)}

# Run a Boosted Regression Tree
# p and a should be data.frames with columns x and y as dec coordinates, envidata environmental raster stacks 
f.BRTrun <- function(p, a, envidata, learning.rate, bag.fraction, data.extent){
  r <- rbind(cbind(p, "occ"=1), cbind(a, "occ"=0))
  brtdata <- cbind(r, extract(envidata, r[,1:2]))
  nichemodel <- gbm.step(data=brtdata, gbm.x=4:length(colnames(brtdata)), gbm.y=3, family="bernoulli", tree.complexity=3, learning.rate=learning.rate, bag.fraction=bag.fraction)
  prediction.values <- predict.gbm(nichemodel, as.data.frame(values(envidata)), n.trees=nichemodel$gbm.call$best.trees, type="response")
  prediction <- envidata[[1]]
  values(prediction) <- prediction.values
  prediction <- mask(prediction, data.extent)
  return(list(nichemodel, prediction))}


#*#*#*#*#*#*#*#*#
## VARIABLES ####
#*#*#*#*#*#*#*#*#

# General
# setwd("/home/jlg57/sb1")
setwd("C:/Users/LaurensG/Desktop/geodist")

# Set max memory limit
# memory.limit(4095)
# For reproducible results
set.seed(1)
# backup graphical parameters
# p.par <- par()
# Set Java virtual machine memory
# options(java.parameters = "-Xmx4g" )

# Select output directory
p.output.dir = paste("R output/",Sys.Date(), "/", sep="")
try(f.create.output.dir(p.output.dir))
#p.testfraction = 0.5

# Set Java virtual machine memory
# options(java.parameters = "-Xmx4g" )


#*#*#*#*#*#*#*#*#*#
## LOAD DATA  #####
#*#*#*#*#*#*#*#*#*#

#species.list <- read.csv("data/speciesinfo.csv", header=T)
#Subset for individual core on HPC nodes
#species.list <- species.list[species.list$Core==1,]
#species.list <- species.list[1:3,]

#species.list <- read.csv("C:/data/modelling/data/speciesinfo.csv", header=T)
#species.iucnmap <- unlist(strsplit(list.files("C:/data/modelling/data/IUCN/singleshapefiles", pattern=".shp$"), split=".shp"))
#species.saupmap <- unlist(strsplit(list.files("C:/data/modelling/data/SAUP/singleshapefiles", pattern=".shp$"), split=".shp"))
#' limit species with IUCN OR SAUP range map
#species.list <- species.list[species.list$TaxonKey %in% species.iucnmap | species.list$TaxonKey %in% species.saupmap,"TaxonKey"]

species.list <- unlist(strsplit(list.files("C:/Users/LaurensG/Desktop/geodist/distraster/", pattern=".grd$"), split=".grd"))

load("data/data.species.RData")
data.extent <- raster("data/data.extent.grd")


#*#*#*#*#*#*#*#*#
## ANALYSIS #####
#*#*#*#*#*#*#*#*#

#limit to species with more than 50 records
#species.list <- species.list[species.list$ncells>=50,]

#get background data
a <- data.species[data.species$set=="background",c("x","y")]


s=1
#done till 597
#128, 147, 373, 598, 601, 605, 628 not enough records 
#179 did n
for(s in 1:length(species.list)){
  
  data.envi <- stack("C:/Users/LaurensG/Desktop/geodist/data/data.envi.grd", paste("C:/Users/LaurensG/Desktop/geodist/distraster/", species.list[s], ".grd", sep=""))
  names(data.envi[[8]]) <- "geodist"
  #get training data
  p <- data.species[data.species$TaxonKey==species.list[s] & data.species$set=="training",c("x","y")]
  #get test data
  e <- data.species[data.species$TaxonKey==species.list[s] & data.species$set=="test",c("x","y")]
  #delete records with NA values
  #p <- p[complete.cases(extract(data.envi[[1]], SpatialPoints(p))),]
  #e <- e[complete.cases(extract(data.envi[[1]], SpatialPoints(e))),]
  print(as.character(species.list[s]))
  
  # ~ Run the models ####
  # Bioclim
  BIOCLIM.niche <- bioclim(p=p, x=data.envi)
  BIOCLIM.pred <- predict(BIOCLIM.niche, data.envi)
  
  # Maximum entropy model
  MAXENT.niche <- f.MAXENTrun(p=p, a=a, envidata=data.envi, p.output.dir=paste(p.output.dir,species.list[s],sep="/"))
  MAXENT.pred <- predict(MAXENT.niche, data.envi)
  
  # Boosted regression trees model
  try(BRT <- f.BRTrun(p, a, data.envi, learning.rate=0.01, bag.fraction=0.5, data.extent))
  try(BRT.niche <- BRT[[1]])
  try(BRT.pred <- BRT[[2]])
  
  # ~ Evaluate models without clipping ----
  if(exists("BRT.pred")){
    
    BIOCLIM.eval <- evaluate(p=e, a=a, model=BIOCLIM.niche, x=data.envi)
    MAXENT.eval <- evaluate(p=e, a=a, model=MAXENT.niche, x=data.envi)
    BRT.eval <- evaluate(p=e, a=a, model=BRT.niche, x=data.envi, n.trees=BRT.niche$gbm.call$best.trees)
    
    # ~ Save model outputs
    results <- list(list(BIOCLIM.niche, MAXENT.niche, BRT.niche), 
                    list(BIOCLIM.pred, MAXENT.pred, BRT.pred), 
                    list(BIOCLIM.eval, MAXENT.eval, BRT.eval))
    save(results, file=paste(p.output.dir, "/SDMeval_geodist_", species.list[s], ".RData", sep=""))
  }else{
    
    BRT.niche="not done"
    BRT.pred="not done"
    BRT.eval="not done"
    
    BIOCLIM.eval <- evaluate(p=e, a=a, model=BIOCLIM.niche, x=data.envi)
    MAXENT.eval <- evaluate(p=e, a=a, model=MAXENT.niche, x=data.envi)
    
    # ~ Save model outputs
    results <- list(list(BIOCLIM.niche, MAXENT.niche, BRT.niche), 
                    list(BIOCLIM.pred, MAXENT.pred, BRT.pred), 
                    list(BIOCLIM.eval, MAXENT.eval, BRT.eval))
    save(results, file=paste(p.output.dir, "/SDMeval_geodist_", species.list[s], ".RData", sep=""))
  }
  
  # ~ Plot species records and predicted distribution ----
  try(f.create.output.dir(paste(p.output.dir, "maps/", sep="")))
  
  png(paste(p.output.dir, "maps/", species.list[s] ,"_occs.png", sep=""), width=10350, height=5400, units="px")
  plot(MAXENT.pred, col=rgb(0,0,255, maxColorValue=255), colNA="black", axes=F, legend=F)#, addfun=c(points(points(p, pch=3, col="red", cex=5)), points(e, pch=3, col="white", cex=5)))
  points(p, pch=19, col="red", cex=4)
  points(e, pch=19, col="white", cex=4)
  dev.off()
  
  png(paste(p.output.dir, "maps/", species.list[s] ,"_BIOCLIM_geodist.png", sep=""), width=3450*3, height=1800*3, units="px")
  par(pty="m")
  plot(BIOCLIM.pred, col=c(rgb(0,0,255, maxColorValue=255), rgb(54,97,255, maxColorValue=255), rgb(56,172,255, maxColorValue=255), rgb(0,255,255, maxColorValue=255), rgb(145,255,180, maxColorValue=255), rgb(210,255,105, maxColorValue=255), rgb(255,255,0, maxColorValue=255), rgb(255,183,0, maxColorValue=255), rgb(255,111,0, maxColorValue=255), rgb(255,0,0, maxColorValue=255)), colNA="black", breaks=c(1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0), axes=F, legend=F)#, addfun=c(points(points(p, pch=3, col="red", cex=5)), points(e, pch=3, col="white", cex=5)))
  dev.off()
  png(paste(p.output.dir, "maps/", species.list[s] ,"_MAXENT_geodist.png", sep=""), width=3450*3, height=1800*3, units="px")
  par(pty="m")
  plot(MAXENT.pred, col=c(rgb(0,0,255, maxColorValue=255), rgb(54,97,255, maxColorValue=255), rgb(56,172,255, maxColorValue=255), rgb(0,255,255, maxColorValue=255), rgb(145,255,180, maxColorValue=255), rgb(210,255,105, maxColorValue=255), rgb(255,255,0, maxColorValue=255), rgb(255,183,0, maxColorValue=255), rgb(255,111,0, maxColorValue=255), rgb(255,0,0, maxColorValue=255)), colNA="black", breaks=c(1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0), axes=F, legend=F)#, addfun=c(points(points(p, pch=3, col="red", cex=5)), points(e, pch=3, col="white", cex=5)))
  dev.off()
  png(paste(p.output.dir, "maps/", species.list[s] ,"_BRT_geodist.png", sep=""), width=3450*3, height=1800*3, units="px")
  par(pty="m")
  try(plot(BRT.pred, col=c(rgb(0,0,255, maxColorValue=255), rgb(54,97,255, maxColorValue=255), rgb(56,172,255, maxColorValue=255), rgb(0,255,255, maxColorValue=255), rgb(145,255,180, maxColorValue=255), rgb(210,255,105, maxColorValue=255), rgb(255,255,0, maxColorValue=255), rgb(255,183,0, maxColorValue=255), rgb(255,111,0, maxColorValue=255), rgb(255,0,0, maxColorValue=255)), colNA="black", breaks=c(1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0), axes=F, legend=F))#, addfun=c(points(points(p, pch=3, col="red", cex=5)), points(e, pch=3, col="white", cex=5)))
  dev.off()
  
  print(paste("###", s, "###"))
  timestamp()
}

print("This core is done!")


