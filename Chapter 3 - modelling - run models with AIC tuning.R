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


  # Sample background and create background points with probability correcting for cell area
  f.sample.background <- function(data.extent, p.number.of.background.points){
    data.background.points <- data.frame(randomPoints(data.extent, n=p.number.of.background.points, lonlatCorrection=T))
    data.background.points <- data.frame(x=data.background.points$x, y=data.background.points$y, set="background", occ=0)
    return(data.background.points)}


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
  f.BRTrun <- function(p, a, envidata, learning.rate=0.01, bag.fraction=0.5, data.extent){
    r <- rbind(cbind(p, "occ"=1), cbind(a, "occ"=0))
    brtdata <- cbind(r, extract(envidata, r[,1:2]))
    nichemodel <- gbm.step(data=brtdata, gbm.x=4:length(colnames(brtdata)), gbm.y=3, family="bernoulli", tree.complexity=3, learning.rate=learning.rate, bag.fraction=bag.fraction)
    prediction.values <- predict.gbm(nichemodel, as.data.frame(values(envidata)), n.trees=nichemodel$gbm.call$best.trees, type="response")
    prediction <- envidata[[1]]
    values(prediction) <- prediction.values
    prediction <- mask(prediction, data.extent)
    return(list(nichemodel, prediction))}


  #http://www.earthskysea.org/!r/maxentAic.r
  f.maxent.aic <- function(data.envi, p, a, e, p.output.dir=NULL, beta.test=c(0.5, 1:20)){
    aicFrame=NULL  
    for(this.beta in beta.test){
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
        print(thisAicFrame)
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
  f.BRT.aic <- function(data.envi, p, a, e, learning.rate, bag.test=c(seq(0.1,0.9,0.1)), data.extent=data.extent){
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
  # setwd("/home/jlg57")
  setwd("C:/data/modelling")

  # Set max memory limit
  # memory.limit(4095)
  # For reproducible results
  #set.seed(1)
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

  species.list <- read.csv("data/speciesinfo.csv", header=T)
  
  load("data/data.species.RData")
  data.envi <- stack("data/data.envi.grd")
  data.extent <- raster("data/data.extent.grd")
  
#*#*#*#*#*#*#*#*#
## ANALYSIS #####
#*#*#*#*#*#*#*#*#

  core=1
  species.list <- species.list[species.list$Core==core,]

  #limit to species with more than 50 records
  species.list <- species.list[species.list$ncells>=50,]

  #get background data
  a <- data.species[data.species$set=="background",c("x","y")]
  #a <- data.species[data.species$set=="background"]
  
  s=1
  for(s in 1:length(species.list$TaxonName)){
    #get training data
    p <- data.species[data.species$TaxonKey==species.list[s,1] & data.species$set=="training",c("x","y")]
    #get test data
    e <- data.species[data.species$TaxonKey==species.list[s,1] & data.species$set=="test",c("x","y")]
    #delete records with NA values
    #p <- p[complete.cases(extract(data.envi[[1]], SpatialPoints(p))),]
    #e <- e[complete.cases(extract(data.envi[[1]], SpatialPoints(e))),]
    print(as.character(species.list[s,2]))
    
    # ~ Run the models ####
    # Bioclim
    BIOCLIM.niche <- bioclim(p=p, x=data.envi)
    BIOCLIM.pred <- predict(BIOCLIM.niche, data.envi)
    
    # Maximum entropy model
    MAXENT.aic <- f.maxent.aic(data.envi, p, a, e, p.output.dir, beta.test=c(0.5, 1, 1.5, 2, 2.5, 3:10))
    MAXENT.niche <- MAXENT.aic[[2]]
    MAXENT.pred <- predict(MAXENT.niche, data.envi)
    
    # Boosted regression trees model
    try(BRT.aic <- f.BRT.aic(data.envi, p, a, e, learning.rate=0.01, bag.test=c(seq(0.1,0.9,0.1)), data.extent=data.extent))
    try(BRT.niche <- BRT.aic[[2]][[1]])  
    try(BRT.pred <- BRT.aic[[2]][[2]])

    # ~ Evaluate models ####
    if(exists("BRT.pred")){

      BIOCLIM.eval <- evaluate(p=e, a=a, model=BIOCLIM.niche, x=data.envi)
      MAXENT.eval <- evaluate(p=e, a=a, model=MAXENT.niche, x=data.envi)
      BRT.eval <- evaluate(p=e, a=a, model=BRT.niche, x=data.envi, n.trees=BRT.niche$gbm.call$best.trees)
      
      # ~ Save model outputs
      results <- list(list(BIOCLIM.niche, MAXENT.niche, BRT.niche), 
                      list(BIOCLIM.pred, MAXENT.pred, BRT.pred), 
                      list(BIOCLIM.eval, MAXENT.eval, BRT.eval),
                      list(MAXENT.aic[[1]],BRT.aic[[1]]))
      save(results, file=paste(p.output.dir, "/SDMeval_aic_", species.list[s,1], ".RData", sep=""))
    }else{
      
      BRT.niche="not done"
      BRT.pred="not done"
      BRT.eval="not done"
      
      BIOCLIM.eval <- evaluate(p=e, a=a, model=BIOCLIM.niche, x=data.envi)
      MAXENT.eval <- evaluate(p=e, a=a, model=MAXENT.niche, x=data.envi)
      
      # ~ Save model outputs
      results <- list(list(BIOCLIM.niche, MAXENT.niche, BRT.niche), 
                      list(BIOCLIM.pred, MAXENT.pred, BRT.pred), 
                      list(BIOCLIM.eval, MAXENT.eval, BRT.eval),
                      list(MAXENT.aic[[1]],BRT.aic[[1]]))
      save(results, file=paste(p.output.dir, "/SDMeval_aic_", species.list[s,1], ".RData", sep=""))}
    
    # ~ Plot species records and predicted distribution
    try(f.create.output.dir(paste(p.output.dir, "maps/", sep="")))
    
    png(paste(p.output.dir, "maps/", species.list[s,1] ,"_occs.png", sep=""), width=10350, height=5400, units="px")
    plot(MAXENT.pred, col=rgb(0,0,255, maxColorValue=255), colNA="black", axes=F, legend=F)#, addfun=c(points(points(p, pch=3, col="red", cex=5)), points(e, pch=3, col="white", cex=5)))
    points(p, pch=19, col="red", cex=4)
    points(e, pch=19, col="white", cex=4)
    dev.off()
    
    png(paste(p.output.dir, "maps/", species.list[s,1] ,"_BIOCLIM_pred.png", sep=""), width=3450*3, height=1800*3, units="px")
    par(pty="m")
    plot(BIOCLIM.pred, col=c(rgb(0,0,255, maxColorValue=255), rgb(54,97,255, maxColorValue=255), rgb(56,172,255, maxColorValue=255), rgb(0,255,255, maxColorValue=255), rgb(145,255,180, maxColorValue=255), rgb(210,255,105, maxColorValue=255), rgb(255,255,0, maxColorValue=255), rgb(255,183,0, maxColorValue=255), rgb(255,111,0, maxColorValue=255), rgb(255,0,0, maxColorValue=255)), colNA="black", breaks=c(1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0), axes=F, legend=F)#, addfun=c(points(points(p, pch=3, col="red", cex=5)), points(e, pch=3, col="white", cex=5)))
    dev.off()
    png(paste(p.output.dir, "maps/", species.list[s,1] ,"_MAXENT_pred.png", sep=""), width=3450*3, height=1800*3, units="px")
    par(pty="m")
    plot(MAXENT.pred, col=c(rgb(0,0,255, maxColorValue=255), rgb(54,97,255, maxColorValue=255), rgb(56,172,255, maxColorValue=255), rgb(0,255,255, maxColorValue=255), rgb(145,255,180, maxColorValue=255), rgb(210,255,105, maxColorValue=255), rgb(255,255,0, maxColorValue=255), rgb(255,183,0, maxColorValue=255), rgb(255,111,0, maxColorValue=255), rgb(255,0,0, maxColorValue=255)), colNA="black", breaks=c(1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0), axes=F, legend=F)#, addfun=c(points(points(p, pch=3, col="red", cex=5)), points(e, pch=3, col="white", cex=5)))
    dev.off()
    png(paste(p.output.dir, "maps/", species.list[s,1] ,"_BRT_pred.png", sep=""), width=3450*3, height=1800*3, units="px")
    par(pty="m")
    try(plot(BRT.pred, col=c(rgb(0,0,255, maxColorValue=255), rgb(54,97,255, maxColorValue=255), rgb(56,172,255, maxColorValue=255), rgb(0,255,255, maxColorValue=255), rgb(145,255,180, maxColorValue=255), rgb(210,255,105, maxColorValue=255), rgb(255,255,0, maxColorValue=255), rgb(255,183,0, maxColorValue=255), rgb(255,111,0, maxColorValue=255), rgb(255,0,0, maxColorValue=255)), colNA="black", breaks=c(1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0), axes=F, legend=F))#, addfun=c(points(points(p, pch=3, col="red", cex=5)), points(e, pch=3, col="white", cex=5)))
    dev.off()
    
    rm(BRT.niche, BRT.pred)
  }
  
  print("This core is done!")
