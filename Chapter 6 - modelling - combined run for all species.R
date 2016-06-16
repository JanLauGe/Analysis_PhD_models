#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#
#R code modelling commercial fish species distributions #
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#

#PhD project
#Jan Laurens Geffert
#Department of Geography
#University of Cambridge
#jlg57@cam.ac.uk

#' Finalmodel
#' This code combines the modifications of the previous chapters to modell all commercially important marine fish species
#' If there is at least n records. Where available expert range map is used, otherwise FAO area is used instead

# Set core
core=1
# Set minimum number of records
n=10
# Set number of k-partitions
k=10


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
f.MAXENTrun <- function(p, a, data.envi, p.output.dir){
  m <- maxent(x=data.envi, p=p, a=a,
              removeDuplicates=T, overwrite=T, path=paste(p.output.dir, "/maxent", sep=""),
              #See Maxent.jar programme -> help for a list and summary of these commands
              args=c("removeduplicates=TRUE", "linear=TRUE", "quadratic=TRUE", "product=TRUE", "threshold=TRUE", "hinge=TRUE",
                     "responsecurves=TRUE", "pictures=FALSE", "plots=FALSE", "writeplotdata=TRUE",
                     "askoverwrite=FALSE", "writeclampgrid=FALSE", "writemess=FALSE", "warnings=TRUE",
                     "allowpartialdata=false", "betamultiplier=1.0"))    
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


#http://www.earthskysea.org/!r/maxentAic.r
f.MAXENT.aic <- function(data.envi, p, a, e, p.output.dir=NULL, beta.test=c(0.5, 1:20)){
  aicFrame=NULL  
  for(this.beta in beta.test){
    #print(paste(this.beta))
    rnumb <- round(runif(1) * 10^12)
    f.create.output.dir(paste(p.output.dir, '/', rnumb, sep=''))
    #train model
    trialModel <- maxent(x=data.envi, p=p, a=a,
                         removeDuplicates=T, overwrite=T, path=paste(p.output.dir, "/", rnumb, sep=""),
                         #See Maxent.jar programme -> help for a list and summary of these commands
                         args=c("removeduplicates=true", "linear=TRUE", "quadratic=TRUE", "product=TRUE", "threshold=TRUE", "hinge=TRUE",
                                "responsecurves=true", "pictures=false", "plots=true", "writeplotdata=true",
                                "askoverwrite=false", "skipifexists=false", "randomseed=false", "randomtestpoints=0", "jackknife=false",
                                "writeclampgrid=false", "writemess=false", "warnings=true", "allowpartialdata=false",
                                paste( 'betamultiplier=', this.beta, sep='' )))  
    
    #predict to training (and maybe test presences)
    predPresLikelihood <- predict(object=trialModel, x=data.envi, na.rm=TRUE, args='outputformat=raw')
    predPresLikelihood <- extract(predPresLikelihood, e)
    
    #predict to raster
    thisRaster <- predict(object=trialModel, x=data.envi, filename=paste(p.output.dir, '/', rnumb, '/mxTuning_beta_', round(100 * this.beta), sep=''), na.rm=TRUE, format='raster', overwrite=TRUE, args='outputformat=raw')
    mapSum <- cellStats(thisRaster, 'sum')
    #calculate log likelihood
    logLik <- sum(log(predPresLikelihood/mapSum), na.rm=TRUE)
    #delete temps
    rm(thisRaster); gc()
    Sys.sleep(1)
    unlink(paste(p.output.dir, '/', rnumb, sep=''), recursive=TRUE)
    
    ## calculate number of parameters
    K=0
    for(thisLambda in Lambdas <- trialModel@lambdas){ # for each line in lambda object
      commaPos <- gregexpr( text=thisLambda, pattern=',') # get location of commas
      if(length(commaPos[[1]]) > 1) { # if there is >1 comma in this line (this is not a parameter line)
        paramValue <- as.numeric( substr(x=thisLambda, start=commaPos[[1]][1]+1, stop=commaPos[[1]][2]-1 ) ) # convert string between first two commas to numeric
        if(paramValue !=0) K <- K + 1}} # increment number of parameters
    
    ## calculate AICc
    AICc <- -2 * logLik + 2 * K + (2 * K * (K + 1)) / (nrow(e) - K - 1)
    thisAicFrame <- data.frame(beta=this.beta, n=nrow(e), logLik=logLik, K=K, AICc=AICc)
    #print(thisAicFrame)
    aicFrame <- rbind(aicFrame, thisAicFrame)}
  #remove models with more parameters than data points that have more than 0 parameters
  #aicFrame <- subset(aicFrame, aicFrame$n >= aicFrame$K & aicFrame$K > 0) 
  aicFrame <- aicFrame[order(aicFrame$AICc, aicFrame$beta),]
  final.beta <- aicFrame[1,"beta"]
  
  finalModel <- maxent(x=data.envi, p=p, a=a,
                       removeDuplicates=T, overwrite=T, path=paste(p.output.dir, "/maxent", sep=""),
                       #See Maxent.jar programme -> help for a list and summary of these commands
                       args=c("removeduplicates=true", "linear=TRUE", "quadratic=TRUE", "product=TRUE", "threshold=TRUE", "hinge=TRUE",
                              "responsecurves=true", "pictures=false", "plots=true", "writeplotdata=true",
                              "askoverwrite=false", "skipifexists=false", "randomseed=false", "randomtestpoints=0", "jackknife=false",
                              "writeclampgrid=false", "writemess=false", "warnings=true", "allowpartialdata=false",
                              paste( 'betamultiplier=', final.beta, sep='' )))
  return(list(aicFrame, finalModel))}


#Compare bag fractions for BRT using AICc
f.BRT.aic <- function(p, a, e, data.envi, learning.rate, bag.test=c(seq(0.1,0.9,0.1)), data.extent=data.extent){
  aicFrame=NULL  
  for(this.bag in bag.test){
    rnumb <- round(runif(1) * 10^12)
    f.create.output.dir(paste(p.output.dir, '/', rnumb, sep=''))
    
    #train model
    trialModel <- f.BRTrun(p, a, data.envi, learning.rate=learning.rate, bag.fraction=this.bag, data.extent=data.extent)
    #predict to training (and maybe test presences)
    predPresLikelihood <- trialModel[[2]]
    predPresLikelihood <- extract(predPresLikelihood, e)
    mapSum <- cellStats(trialModel[[2]], 'sum')
    
    #calculate log likelihood
    logLik <- sum(log(predPresLikelihood/mapSum), na.rm=TRUE)
    #calculate number of parameters
    K=trialModel[[1]]$gbm.call$best.trees
    #calculate AICc
    AICc <- -2 * logLik + 2 * K + (2 * K * (K + 1)) / (nrow(e) - K - 1)
    thisAicFrame <- data.frame(bag=this.bag, n=nrow(p), logLik=logLik, K=K, AICc=AICc)
    print(thisAicFrame)
    aicFrame <- rbind(aicFrame, thisAicFrame)}
  aicFrame <- aicFrame[order(aicFrame$AICc, aicFrame$bag),]
  
  final.bag <- aicFrame[1,"bag"]
  finalModel <- f.BRTrun(p, a, data.envi, learning.rate=learning.rate, bag.fraction=final.bag, data.extent=data.extent)
  return(list(aicFrame, finalModel))}


#*#*#*#*#*#*#*#*#
## VARIABLES ####
#*#*#*#*#*#*#*#*#

# General
# Set wd
setwd("C:/Users/LaurensG/Desktop/finalmodel")
# setwd("/home/jlg57/sb1")
# Set output dir
p.output.dir = paste(getwd(), "/outputs", sep="")
try(f.create.output.dir(p.output.dir))
# For reproducible results
set.seed(1)
# backup graphical parameters
p.par <- par()


#*#*#*#*#*#*#*#*#*#
## LOAD DATA  -----
#*#*#*#*#*#*#*#*#*#

# species info
species.list <- read.csv("data/species.list.csv")#unlist(strsplit(list.files("C:/Users/LaurensG/Desktop/geodist/distraster/", pattern=".grd$"), split=".grd"))
# raster extent
data.extent <- raster("data/data.extent.grd")
# environmental variables
data.envi <- stack("data/data.envi.grd")
# FAO raster
data.faoraster <- raster("data/data.faoraster.grd")
# FAO info
data.faoareas <- read.csv("data/data.faoareas.csv")

# species data as data.species
load("data/data.species.RData")
# background data as a
load("data/data.background.sb1.RData")

#*#*#*#*#*#*#*#*#*#*#
## PRE-ANALYSIS -----
#*#*#*#*#*#*#*#*#*#*#

#limit to species with n records
count <- data.frame("TaxonKey"=aggregate(data.species[,2], by=list(data.species$TaxonKey), FUN="length")[,1],
                    "ncells"=aggregate(data.species[,2], by=list(data.species$TaxonKey), FUN="length")[,2])
species.list <- merge(species.list, count, by="TaxonKey", all.x=T)
species.list <- species.list[species.list$ncells>=n,]
species.list <- species.list[species.list$Core==core,]

# RUN FOR SPECIES ----
s=1
for(s in 1:length(species.list$TaxonKey)){
  summ=data.frame("TaxonKey"=as.numeric(NULL), "k"=as.numeric(NULL), "beta"=as.numeric(NULL), "npar"=as.numeric(NULL), "AICc"=as.numeric(NULL), "rangemap"=as.logical(NULL), 
                  "BIO.auc"=as.numeric(NULL), "MAX.auc"=as.numeric(NULL), "BRT.auc"=as.numeric(NULL), "BIO.cor"=as.numeric(NULL), "MAX.cor"=as.numeric(NULL), "BRT.cor"=as.numeric(NULL))
  print(paste("### ", species.list[s,"TaxonKey"], " - ", as.character(species.list[s,"TaxonName"])))
  if(file.exists(paste(p.output.dir, "/SDMeval_finalmodel_", species.list[s,"TaxonKey"], ".RData", sep=""))){
    print(paste(species.list[s,"TaxonKey"], "~ already done"))
  }else{
    print(paste(species.list[s,"TaxonKey"], "~ starting"))
    # reset aic
    MAXENT.aic=NULL
      
    # Add expert range if exists
    if(file.exists(paste("data/distraster/", species.list[s,"TaxonKey"], ".grd", sep=""))){
      data.envi.s <- stack(data.envi, paste("data/distraster/", species.list[s,"TaxonKey"], ".grd", sep=""))
      names(data.envi.s[[8]]) <- "rangemap"
      print(paste(species.list[s,"TaxonKey"], "~ with range map"))
      rangemap.s=TRUE
    }else{
      data.envi.s <- data.envi
      print(paste(species.list[s,"TaxonKey"], "~ no range map"))
      rangemap.s=FALSE
    }
  
    #get training data
    p <- data.species[data.species$TaxonKey==species.list[s,"TaxonKey"],c("x","y")]
    
    # ~ k-fold partitioning of species data ----
    pk <- kfold(p, k)
    krun=1
    for(krun in 1:k){
      print(paste(species.list[s,"TaxonKey"], "#k", krun))
      train <- p[!pk==krun,]
      test <- p[pk==krun,]
    
      # ~ Run the models ----
      # Bioclim
      print(paste(species.list[s,"TaxonKey"], "~ running BIOCLIM..."))
      BIOCLIM.niche <- bioclim(p=train, x=data.envi.s)
      BIOCLIM.pred <- predict(BIOCLIM.niche, data.envi.s)
      print(paste(species.list[s,"TaxonKey"], "~ done"))
      
      # Maximum entropy model
      print(paste(species.list[s,"TaxonKey"], "~ running MAXENT..."))
      MAXENT.niche <- f.MAXENT.aic(p=train, e=test, a=a, data.envi=data.envi.s, p.output.dir=paste(p.output.dir,species.list[s,"TaxonKey"],sep="/"), beta.test=c(0.5, 1, 1.5, 2, 2.5, 3:20))
      MAXENT.pred <- predict(MAXENT.niche[[2]], data.envi.s)
      # Get optimal AICc value
      MAXENT.aic <- c(MAXENT.aic, list(MAXENT.niche[[1]]))
      beta.s <- MAXENT.aic[[1]][1,1]
      AICc.s <- MAXENT.aic[[1]][1,5]
      npar.s <- MAXENT.aic[[1]][1,4]
      print(paste(species.list[s,"TaxonKey"], "~ done"))
      
      # Boosted regression trees model
      print(paste(species.list[s,"TaxonKey"], "~ running BRT..."))
      # sink to make output silent
      f=file()
      sink(file=f)
      BRT <- f.BRT.aic(p=train, a=a, e=train, data.envi=data.envi.s, learning.rate=0.01, bag.test=0.5, data.extent)
      sink()
      close(f)
      print(paste(species.list[s,"TaxonKey"], "~ done"))
      BRT.niche <- BRT[[2]][[1]]
      BRT.pred <- BRT[[2]][[2]]
      
      # ~ Clip models without range maps ----
      if(rangemap.s==FALSE){
        faoareas.s <- data.faoareas[data.faoareas$TaxonKey==species.list[s,1],3:28]
        if(exists("faoareas.s")){
          print(paste(species.list[s,"TaxonKey"], "~ clipping by FAO area"))    
          faoraster.s <- data.faoraster %in% faoareas.s
          faoraster.s[!faoraster.s %in% c(1)] <- NA
          
          BIOCLIM.pred <- mask(BIOCLIM.pred, faoraster.s, updatevalue=0)
          MAXENT.pred <- mask(MAXENT.pred, faoraster.s, updatevalue=0)
          BRT.pred <- mask(BRT.pred, faoraster.s, updatevalue=0)
        }else if(!exists(faoareas.s)){
          print(paste(species.list[s,"TaxonKey"], "!FAO area info not found"))
          break
        }else{
          print(paste(species.list[s,"TaxonKey"], "!FAO info error"))}
        
      }else if(rangemap.s==TRUE){
        print(paste(species.list[s,"TaxonKey"], "~ no clipping"))
      }else{
        print(paste(species.list[s,"TaxonKey"], "!Error with range map"))
        break}
      
      # ~ Evaluate models without clipping ----
      BIOCLIM.eval <- evaluate(p=test, a=a, model=BIOCLIM.niche, x=data.envi.s)
      MAXENT.eval <- evaluate(p=test, a=a, model=MAXENT.niche[[2]], x=data.envi.s)
      BRT.eval <- evaluate(p=test, a=a, model=BRT.niche, x=data.envi.s, n.trees=BRT.niche$gbm.call$best.trees)
      
      # ~ Save model outputs ----
      assign(paste("results", krun, sep=""), list(list(BIOCLIM.niche, MAXENT.niche[[2]], BRT.niche), 
                                              list(BIOCLIM.pred, MAXENT.pred, BRT.pred), 
                                              list(BIOCLIM.eval, MAXENT.eval, BRT.eval)))
      summ.s <- data.frame("TaxonKey"=species.list[s,"TaxonKey"], "k"=krun, "beta"=beta.s, "npar"=npar.s, "AICc"=AICc.s, "rangemap"=rangemap.s, 
                           "BIO.auc"=BIOCLIM.eval@auc, "MAX.auc"=MAXENT.eval@auc, "BRT.auc"=BRT.eval@auc, "BIO.cor"=BIOCLIM.eval@cor[[1]], "MAX.cor"=MAXENT.eval@cor[[1]], "BRT.cor"=BRT.eval@cor[[1]])
      summ <- rbind(summ, summ.s)
    }
    results <- list(results1, results2, results3, results4, results5)
    save(results, file=paste(p.output.dir, "/SDMeval_finalmodel_", species.list[s,"TaxonKey"], ".RData", sep=""), overwrite=T)
    write.csv(summ, file=paste(p.output.dir, "/AICc_finalmodels_", species.list[s,"TaxonKey"], ".csv", sep=""))
    print(paste(species.list[s,"TaxonKey"], "FINISHED"))
  }
rm(faoraster.s, summ.s)
}

print("This core is done!")


