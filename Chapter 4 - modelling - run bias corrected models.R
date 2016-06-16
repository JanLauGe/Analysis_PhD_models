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
  setwd("C:/Users/LaurensG/Desktop/sdm")

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

  species.list <- read.csv("data/speciesinfo2.csv", header=T)
  #Subset for individual core on HPC nodes
  #species.list <- species.list[species.list$Core==1,]
  #species.list <- species.list[1:3,]

  load("data/data.faoareas.RData")
  load("data/data.species.RData")
  data.envi <- stack("data/data.envi.grd")
  data.faoraster <- raster("data/data.faoraster.grd")
  data.extent <- raster("data/data.extent.grd")


#*#*#*#*#*#*#*#*#
## ANALYSIS #####
#*#*#*#*#*#*#*#*#

  #limit to species with more than 50 records
  #species.list <- species.list[species.list$ncells>=50,]
  species.list <- species.list[species.list$TaxonKey %in% c(690434, 690681, 603574, 600487, 605537, 690024, 690105),]

  #get background data
  #a <- f.sample.background(data.envi[[1]], 10000)[,1:2]
  #a <- data.species[data.species$set=="background",c("x","y")]
  
# #~Create biased background data from all species data
#   plot(data.species[!data.species$set=="background",c("x","y")])
#   bg <- data.species[!data.species$set=="background",c("x","y")]
#   bgselect <- sample(1:length(bg[,1]), size=10000)
#   a <- bg[bgselect,]
#   
#   plot(a)
#   save(a, file="data/data.background.sb1.RData")
#   rm(bg, bgselect)

  #load("C:/Users/LaurensG/Desktop/sampling bias/sb1/data/data.background.sb1.RData")
  load("data/data.background.sb1.RData")

  s=1
  for(s in 1:length(species.list$TaxonKey)){
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
    MAXENT.niche <- f.MAXENTrun(p=p, a=a, envidata=data.envi, p.output.dir=paste(p.output.dir,species.list[s,1],sep="/"))
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
      save(results, file=paste(p.output.dir, "/SDMeval_fulld_sb1_", species.list[s,1], ".RData", sep=""))
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
      save(results, file=paste(p.output.dir, "/SDMeval_sb1_", species.list[s,1], ".RData", sep=""))
}

# ~ Plot species records and predicted distribution ----
    try(f.create.output.dir(paste(p.output.dir, "maps/", sep="")))
    
    png(paste(p.output.dir, "maps/", species.list[s,1] ,"_occs.png", sep=""), width=10350, height=5400, units="px")
    plot(MAXENT.pred, col=rgb(0,0,255, maxColorValue=255), colNA="black", axes=F, legend=F)#, addfun=c(points(points(p, pch=3, col="red", cex=5)), points(e, pch=3, col="white", cex=5)))
    points(p, pch=19, col="red", cex=4)
    points(e, pch=19, col="white", cex=4)
    dev.off()
    
    png(paste(p.output.dir, "maps/", species.list[s,1] ,"_BIOCLIM_pred_sb1.png", sep=""), width=3450*3, height=1800*3, units="px")
    par(pty="m")
    plot(BIOCLIM.pred, col=c(rgb(0,0,255, maxColorValue=255), rgb(54,97,255, maxColorValue=255), rgb(56,172,255, maxColorValue=255), rgb(0,255,255, maxColorValue=255), rgb(145,255,180, maxColorValue=255), rgb(210,255,105, maxColorValue=255), rgb(255,255,0, maxColorValue=255), rgb(255,183,0, maxColorValue=255), rgb(255,111,0, maxColorValue=255), rgb(255,0,0, maxColorValue=255)), colNA="black", breaks=c(1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0), axes=F, legend=F)#, addfun=c(points(points(p, pch=3, col="red", cex=5)), points(e, pch=3, col="white", cex=5)))
    dev.off()
    png(paste(p.output.dir, "maps/", species.list[s,1] ,"_MAXENT_pred_sb1.png", sep=""), width=3450*3, height=1800*3, units="px")
    par(pty="m")
    plot(MAXENT.pred, col=c(rgb(0,0,255, maxColorValue=255), rgb(54,97,255, maxColorValue=255), rgb(56,172,255, maxColorValue=255), rgb(0,255,255, maxColorValue=255), rgb(145,255,180, maxColorValue=255), rgb(210,255,105, maxColorValue=255), rgb(255,255,0, maxColorValue=255), rgb(255,183,0, maxColorValue=255), rgb(255,111,0, maxColorValue=255), rgb(255,0,0, maxColorValue=255)), colNA="black", breaks=c(1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0), axes=F, legend=F)#, addfun=c(points(points(p, pch=3, col="red", cex=5)), points(e, pch=3, col="white", cex=5)))
    dev.off()
    png(paste(p.output.dir, "maps/", species.list[s,1] ,"_BRT_pred_sb1.png", sep=""), width=3450*3, height=1800*3, units="px")
    par(pty="m")
    try(plot(BRT.pred, col=c(rgb(0,0,255, maxColorValue=255), rgb(54,97,255, maxColorValue=255), rgb(56,172,255, maxColorValue=255), rgb(0,255,255, maxColorValue=255), rgb(145,255,180, maxColorValue=255), rgb(210,255,105, maxColorValue=255), rgb(255,255,0, maxColorValue=255), rgb(255,183,0, maxColorValue=255), rgb(255,111,0, maxColorValue=255), rgb(255,0,0, maxColorValue=255)), colNA="black", breaks=c(1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0), axes=F, legend=F))#, addfun=c(points(points(p, pch=3, col="red", cex=5)), points(e, pch=3, col="white", cex=5)))
    dev.off()


#~Evaluate models after clipping ----
  
  # Clip predictions
    s.faoareas <- unique(as.numeric(data.faoareas[data.faoareas$TaxonKey==species.list$TaxonKey[s],3:28]))
    s.faoraster <- data.faoraster %in% s.faoareas
    s.faoraster[!s.faoraster %in% c(1)] <- NA

    BIOCLIM.pred <- mask(BIOCLIM.pred, s.faoraster, updatevalue=0)
    MAXENT.pred <- mask(MAXENT.pred, s.faoraster, updatevalue=0)
    try(BRT.pred <- mask(BRT.pred, s.faoraster, updatevalue=0))


  if(exists("BRT.pred")){
    BIOCLIM.eval <- evaluate(p=extract(BIOCLIM.pred, e), a=extract(BIOCLIM.pred, a))
    MAXENT.eval <- evaluate(p=extract(MAXENT.pred, e), a=extract(MAXENT.pred, a))
    BRT.eval <- evaluate(p=extract(BRT.pred, e), a=extract(BRT.pred, a), n.trees=BRT.niche$gbm.call$best.trees)
    
    # ~ Save model outputs
    results <- list(list(BIOCLIM.niche, MAXENT.niche, BRT.niche), 
                    list(BIOCLIM.pred, MAXENT.pred, BRT.pred), 
                    list(BIOCLIM.eval, MAXENT.eval, BRT.eval))
    save(results, file=paste(p.output.dir, "/SDMeval.clipped_sb1_", species.list[s,1], ".RData", sep=""))
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
    save(results, file=paste(p.output.dir, "/SDMeval_sb1_", species.list[s,1], ".RData", sep=""))}
  
# ~ Plot species records and clipped predicted distribution ----
    try(f.create.output.dir(paste(p.output.dir, "maps/", sep="")))
    
    png(paste(p.output.dir, "maps/", species.list[s,1] ,".clipped_BIOCLIM_pred_sb1.png", sep=""), width=3450*3, height=1800*3, units="px")
    par(pty="m")
    plot(BIOCLIM.pred, col=c(rgb(0,0,255, maxColorValue=255), rgb(54,97,255, maxColorValue=255), rgb(56,172,255, maxColorValue=255), rgb(0,255,255, maxColorValue=255), rgb(145,255,180, maxColorValue=255), rgb(210,255,105, maxColorValue=255), rgb(255,255,0, maxColorValue=255), rgb(255,183,0, maxColorValue=255), rgb(255,111,0, maxColorValue=255), rgb(255,0,0, maxColorValue=255)), colNA="black", breaks=c(1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0), axes=F, legend=F)#, addfun=c(points(points(p, pch=3, col="red", cex=5)), points(e, pch=3, col="white", cex=5)))
    dev.off()
    png(paste(p.output.dir, "maps/", species.list[s,1] ,".clipped_MAXENT_pred_sb1.png", sep=""), width=3450*3, height=1800*3, units="px")
    par(pty="m")
    plot(MAXENT.pred, col=c(rgb(0,0,255, maxColorValue=255), rgb(54,97,255, maxColorValue=255), rgb(56,172,255, maxColorValue=255), rgb(0,255,255, maxColorValue=255), rgb(145,255,180, maxColorValue=255), rgb(210,255,105, maxColorValue=255), rgb(255,255,0, maxColorValue=255), rgb(255,183,0, maxColorValue=255), rgb(255,111,0, maxColorValue=255), rgb(255,0,0, maxColorValue=255)), colNA="black", breaks=c(1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0), axes=F, legend=F)#, addfun=c(points(points(p, pch=3, col="red", cex=5)), points(e, pch=3, col="white", cex=5)))
    dev.off()
    png(paste(p.output.dir, "maps/", species.list[s,1] ,".clipped_BRT_pred_sb1.png", sep=""), width=3450*3, height=1800*3, units="px")
    par(pty="m")
    try(plot(BRT.pred, col=c(rgb(0,0,255, maxColorValue=255), rgb(54,97,255, maxColorValue=255), rgb(56,172,255, maxColorValue=255), rgb(0,255,255, maxColorValue=255), rgb(145,255,180, maxColorValue=255), rgb(210,255,105, maxColorValue=255), rgb(255,255,0, maxColorValue=255), rgb(255,183,0, maxColorValue=255), rgb(255,111,0, maxColorValue=255), rgb(255,0,0, maxColorValue=255)), colNA="black", breaks=c(1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0), axes=F, legend=F))#, addfun=c(points(points(p, pch=3, col="red", cex=5)), points(e, pch=3, col="white", cex=5)))
    dev.off()
}
  
  print("This core is done!")


###
a <- as.data.frame(all.eval1)
#new table
bla <- rbind(data.frame("TaxonID"=as.numeric(rownames(a)), "auc"=as.numeric(a$BIO.bas.auc), "model"=as.factor("BIOCLIM"), "run"=as.factor("baseline")),
             data.frame("TaxonID"=as.numeric(rownames(a)), "auc"=as.numeric(a$BIO.cli.auc), "model"=as.factor("BIOCLIM"), "run"=as.factor("clipped")),
             data.frame("TaxonID"=as.numeric(rownames(a)), "auc"=as.numeric(a$BIO.geo.auc), "model"=as.factor("BIOCLIM"), "run"=as.factor("expert")),
             data.frame("TaxonID"=as.numeric(rownames(a)), "auc"=as.numeric(a$MAX.bas.auc), "model"=as.factor("MAXENT"), "run"=as.factor("baseline")),
             data.frame("TaxonID"=as.numeric(rownames(a)), "auc"=as.numeric(a$MAX.cli.auc), "model"=as.factor("MAXENT"), "run"=as.factor("clipped")),
             data.frame("TaxonID"=as.numeric(rownames(a)), "auc"=as.numeric(a$MAX.geo.auc), "model"=as.factor("MAXENT"), "run"=as.factor("expert")),
             data.frame("TaxonID"=as.numeric(rownames(a)), "auc"=as.numeric(a$BRT.bas.auc), "model"=as.factor("BRT"), "run"=as.factor("baseline")),
             data.frame("TaxonID"=as.numeric(rownames(a)), "auc"=as.numeric(a$BRT.cli.auc), "model"=as.factor("BRT"), "run"=as.factor("clipped")),
             data.frame("TaxonID"=as.numeric(rownames(a)), "auc"=as.numeric(a$BRT.geo.auc), "model"=as.factor("BRT"), "run"=as.factor("expert"))
)

library(ggplot2)

p <- ggplot(data=bla)
p + geom_boxplot(aes(x=model, y=auc, z=run, fill=model), notch=T, width = 0.5, outlier.shape=21, position=position_dodge(0.5)) + 
  theme_bw() + theme(axis.text.x=element_text(color = "black", size=11, angle=0, vjust=.8, hjust=0.8)) +
  labs(title = "") + ylab("AUC") + xlab("") + scale_x_discrete(labels=c("BIOCLIM" = c("bli","bla","blubb"), "MAXENT" = "Dose 1", "BRT" = "Dose 2"))

