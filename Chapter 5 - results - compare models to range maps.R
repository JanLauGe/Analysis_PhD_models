# LIBRARIES ####
library(sp)
library(raster)
library(rgdal)
library(dismo)

#################
##VARIABLES #####
#*#*#*#*#*#*#*#*#

  setwd("C:/data/modelling")
  species.list <- read.csv("C:/data/modelling/data/speciesinfo.csv", header=T)
  data.extent <- raster("C:/data/modelling/data/data.extent.grd")
  
  #limit species selection to species modelled and with range map available
  #list.files("C:/data/modelling/output/bias1", pattern="_fulld_sb1")
  #species.modelled <- unlist(strsplit(unlist(strsplit(list.files("C:/data/modelling/output/bias1", pattern="_fulld_sb1"), split="SDMeval_fulld_sb1_")), split=".RData"))
  #species.expertmap <- unlist(strsplit(list.files("C:/data/modelling/data/iucn/SAUP/singleshapefiles", pattern=".shp$"), split=".shp"))
  #species.saupmap <- unlist(strsplit(list.files("C:/data/modelling/data/SAUP/singleshapefiles", pattern=".shp$"), split=".shp"))

  species.modelled <- unlist(strsplit(unlist(strsplit(list.files("C:/data/modelling/output/bias1", pattern="_fulld_sb1"), split="SDMeval_fulld_sb1_")), split=".RData"))
  species.iucnmap <- unlist(strsplit(list.files("C:/data/modelling/data/IUCN/raster/", pattern=".grd$"), split=".grd"))
  species.saupmap <- unlist(strsplit(list.files("C:/data/modelling/data/SAUP/raster", pattern=".grd$"), split=".grd"))


  #Species with model predictions and IUCN range map
  species.list.s <- species.list[species.list$TaxonKey %in% species.modelled & species.list$TaxonKey %in% species.iucnmap,"TaxonKey"]
  #Species with model predictions and SAUP range map
  species.list.s <- species.list[species.list$TaxonKey %in% species.modelled & species.list$TaxonKey %in% species.saupmap,"TaxonKey"]
  #Species with model predictions and SAUP AND IUCN range map
  species.list.s <- species.list[species.list$TaxonKey %in% species.modelled & (species.list$TaxonKey %in% species.saupmap | species.list$TaxonKey %in% species.iucnmap),"TaxonKey"]

  load("data/data.species.RData")
  data.sample <- data.species[data.species$set=="background",c("x", "y")]

#*#*#*#*#*#*#*#*#*#*#
# IUCN RANGE MAP ####
#*#*#*#*#*#*#*#*#*#*#

  all.auc.con=NULL
  all.cor.con=NULL
  all.exp.cor=NULL
  all.auc.exp=NULL
  all.cor.exp=NULL
  #path="data/IUCN/SAUP/"
  path="data/SAUP/"
  
  # SAUP for 600231 broken
  # Read IUCN range maps and convert to polygons
  for(s in 78:length(species.list.s)){
    if(file.exists(paste(path, "raster/", species.list.s[s], ".grd", sep=""))){
      print(paste(species.list.s[s], " already exists!", sep=""))
      spec.raster <- raster(paste(path, "raster/", species.list.s[s], sep=""))
    }else{
      print(species.list.s[s])
      spec.poly <- readOGR(paste(path, "singleshapefiles", sep=""), species.list.s[s])
      spec.raster <- rasterize(spec.poly, data.extent, background=0, getCover=T)
      spec.raster <- mask(spec.raster, data.extent)
      writeRaster(spec.raster, paste(path, "raster/", species.list.s[s], sep=""), overwrite=F)
    }
    spec.raster <- spec.raster>0
    
  # Read model output
    load(paste("C:/data/modelling/output/baseline/SDMeval_fulld_", species.list.s[s], ".RData", sep=""))
    spec.model.normal <- results
    load(paste("C:/data/modelling/output/bias1/SDMeval_fulld_sb1_", species.list.s[s], ".RData", sep=""))
    spec.model.bias <- results
  
  # normalize HSV from raster ----
    BIO = spec.model.normal[[2]][[1]]
    MAX = spec.model.normal[[2]][[2]]
    MAX.b = spec.model.bias[[2]][[2]]
    BRT = spec.model.normal[[2]][[3]]
    BRT.b = spec.model.bias[[2]][[3]]
    EO = spec.raster
    for(m in c("BIO", "MAX", "MAX.b", "BRT", "BRT.b")){
      hsv <- get(m)
      assign(m, (hsv - cellStats(hsv, min)) / (cellStats(hsv, max) - cellStats(hsv, min)))
    }
  
  
  # Get conventional AUC and COR values for comparison
    spec.auc.con <- cbind("BIO.auc"=spec.model.normal[[3]][[1]]@auc, "MAX.auc"=spec.model.normal[[3]][[2]]@auc, "MAX.b.auc"=spec.model.bias[[3]][[2]]@auc, "BRT.auc"=spec.model.normal[[3]][[3]]@auc, "BRT.b.auc"=spec.model.bias[[3]][[3]]@auc)
    spec.cor.con <- cbind("BIO.auc"=spec.model.normal[[3]][[1]]@cor, "MAX.auc"=spec.model.normal[[3]][[2]]@cor, "MAX.b.auc"=spec.model.bias[[3]][[2]]@cor, "BRT.auc"=spec.model.normal[[3]][[3]]@cor, "BRT.b.auc"=spec.model.bias[[3]][[3]]@cor)
    row.names(spec.auc.con) <- species.list.s[s]
    row.names(spec.cor.con) <- species.list.s[s]
  
  #create presence and absence points from expert range map for evaluation ----
    species.sample <- cbind(data.sample, exp=extract(spec.raster, data.sample))
    species.sample <- species.sample[complete.cases(species.sample),]
    species.sample <- SpatialPointsDataFrame(species.sample[,1:2], species.sample)
    #plot(species.sample, col=species.sample@data[,3])
  
    p <- species.sample[species.sample@data[,3]==1,]
    a <- species.sample[species.sample@data[,3]==0,]
  
    spec.auc.exp=NULL
    spec.cor.exp=NULL
    for(m in c("BIO", "MAX", "MAX.b", "BRT", "BRT.b")){
      assign(paste(m, ".p", sep=""), extract(get(m), p))
      assign(paste(m, ".a", sep=""), extract(get(m), a))
      eval.auc <- evaluate(get(paste(m, ".p", sep="")), get(paste(m, ".a", sep="")))@auc
      eval.cor <- evaluate(get(paste(m, ".p", sep="")), get(paste(m, ".a", sep="")))@cor
      spec.auc.exp <- cbind(spec.auc.exp, eval.auc)
      spec.cor.exp <- cbind(spec.cor.exp, eval.cor)
    }
  colnames(spec.auc.exp) <- c("BIO", "MAX", "MAX.b", "BRT", "BRT.b")
  colnames(spec.cor.exp) <- c("BIO", "MAX", "MAX.b", "BRT", "BRT.b")
    
  all.auc.exp <- rbind(all.auc.exp, spec.auc.exp)
  all.cor.exp <- rbind(all.cor.exp, spec.cor.exp)
  all.auc.con <- rbind(all.auc.con, spec.auc.con)
  all.cor.con <- rbind(all.cor.con, spec.cor.con)
  }

##########


opar <- par()
par(mfrow=c(1,2))

# compare AUC values for standard and TGBA
  png(file="C:/Users/LaurensG/Desktop/boxplot1_auccor_bias_saup_wiucn.png", width=40, height=30, units="cm", res=300)
  par(mfrow=c(2,2))
  boxplot(all.auc.con, names=c("BIO", "MAX standard", "MAX TGBA", "BRT standard", "BRT TGBA"))
  boxplot(all.cor.con, names=c("BIO", "MAX standard", "MAX TGBA", "BRT standard", "BRT TGBA"))

  boxplot(all.auc.exp, names=c("BIO", "MAX standard", "MAX TGBA", "BRT standard", "BRT TGBA"))
  boxplot(all.cor.exp, names=c("BIO", "MAX standard", "MAX TGBA", "BRT standard", "BRT TGBA"))
  dev.off()

  t.test(all.auc.exp[,2], all.auc.exp[,3], paired=T)
  t.test(all.cor.exp[,2], all.cor.exp[,3], paired=T)
  t.test(all.auc.exp[,4], all.auc.exp[,5], paired=T)
  t.test(all.cor.exp[,4], all.cor.exp[,5], paired=T)

  speciesinfo <- read.csv("C:/data/modelling/data/speciesinfo2.csv")
  rownames(speciesinfo) <- speciesinfo$TaxonKey
  rownames(all.auc.exp) <- species.list.s

  speciesinfo2 <- cbind(species.list[species.list$TaxonKey %in% rownames(all.auc.exp),c(1,4,5)], speciesinfo[rownames(speciesinfo) %in% rownames(all.auc.exp),c(2,4)], all.auc.exp, all.cor.exp)
  colnames(speciesinfo2) <- c(colnames(speciesinfo2[1:10]), "cor.BIO", "cor.MAX", "cor.MAX.b", "cor.BRT", "cor.BRT.b")
  speciesinfo2 <- speciesinfo2[speciesinfo2$lifestyle %in% c("bathydemersal", "bathypelagic", "bethopelagic", "demersal", "pelagic-neritic", "pelagic-oceanic", "reef-associated"),]
  speciesinfo2 <- droplevels(speciesinfo2)

t.test(speciesinfo2[speciesinfo2$lifestyle=="bathydemersal",8] - speciesinfo2[speciesinfo2$lifestyle=="bathydemersal",7])
t.test(speciesinfo2[speciesinfo2$lifestyle=="bathypelagic",8] - speciesinfo2[speciesinfo2$lifestyle=="bathypelagic",7])
t.test(speciesinfo2[speciesinfo2$lifestyle=="demersal",8] - speciesinfo2[speciesinfo2$lifestyle=="demersal",7])
t.test(speciesinfo2[speciesinfo2$lifestyle=="pelagic-neritic",8] - speciesinfo2[speciesinfo2$lifestyle=="pelagic-neritic",7])
t.test(speciesinfo2[speciesinfo2$lifestyle=="pelagic-oceanic",8] - speciesinfo2[speciesinfo2$lifestyle=="pelagic-oceanic",7])
t.test(speciesinfo2[speciesinfo2$lifestyle=="reef-associated",8] - speciesinfo2[speciesinfo2$lifestyle=="reef-associated",7])


diff <- as.data.frame(cbind(as.numeric(speciesinfo2[,8] - speciesinfo2[,7]), speciesinfo2$lifestyle))


as.float(diff[diff$lifestyle=="bathydemersal",1])
#Difference in AUC
boxplot(speciesinfo2[,8] - speciesinfo2[,7] ~ speciesinfo2$lifestyle)
#Difference in COR
boxplot(speciesinfo2[,13] - speciesinfo2[,12] ~ speciesinfo2$lifestyle)


  par(mfrow=c(2,1))
  boxplot(speciesinfo2[,7] ~ speciesinfo2$lifestyle)
  boxplot(speciesinfo2[,8] ~ speciesinfo2$lifestyle)
  dev.off()

  par(mfrow=c(2,1))
  boxplot(speciesinfo2[,10] ~ speciesinfo2$lifestyle)
  boxplot(speciesinfo2[,11] ~ speciesinfo2$lifestyle)
  dev.off()

###########

#############
# compare HSVs average, inside, outside of expert range
  colnames(all.mean) <- c("ma.BIO", "ma.MAX", "ma.MAX.b", "ma.BRT", "ma.BRT.b",
                          "mi.BIO", "mi.MAX", "mi.MAX.b", "mi.BRT", "mi.BRT.b",
                          "mo.BIO", "mo.MAX", "mo.MAX.b", "mo.BRT", "mo.BRT.b",
                          "va.BIO", "va.MAX", "va.MAX.b", "va.BRT", "va.BRT.b",
                          "vi.BIO", "vi.MAX", "vi.MAX.b", "vi.BRT", "vi.BRT.b",
                          "vo.BIO", "vo.MAX", "vo.MAX.b", "vo.BRT", "vo.BRT.b")
  all.mean <- as.data.frame(all.mean)

  speciesinfo <- read.csv("C:/data/modelling/data/speciesinfo2.csv")
  rownames(speciesinfo) <- speciesinfo$TaxonKey
  specinfo <- cbind(speciesinfo[rownames(speciesinfo) %in% rownames(all.mean),], all.mean)


bla <- aggregate(cbind("mean"=ma.MAX, "mean.b"=ma.MAX.b, "in"=mi.MAX, "in.b"=mi.MAX.b, "inr"=mi.MAX/ma.MAX, "inr.b"=mi.MAX.b/ma.MAX.b, "out"=mo.MAX/ma.MAX, "out.b"=mo.MAX.b/ma.MAX.b) ~ lifestyle, mean, data=specinfo)
bla <- aggregate(cbind("mean"=ma.MAX, "mean.b"=ma.MAX.b, "in"=mi.MAX/ma.MAX, "in.b"=mi.MAX.b/ma.MAX.b, "out"=mo.MAX/ma.MAX, "out.b"=mo.MAX.b/ma.MAX.b) ~ lifestyle, mean, data=specinfo)


#match based on taxonkey


  boxplot(all.mean[,c("ma.MAX", "ma.MAX.b", "ma.BRT", "ma.BRT.b")])

  boxplot( 
    all.mean[,7], all.mean[,8],
    all.mean[,9], all.mean[,10],
    names=c("MAXENT standard", "MAXENT TGBA", "BRT standard", "BRT TGBA"),
    title="Mean HSV inside")

  boxplot(
    all.mean[,12], all.mean[,13], 
    all.mean[,14], all.mean[,15],
    names=c("MAXENT standard", "MAXENT TGBA", "BRT standard", "BRT TGBA"),
    title="Mean HSV outside")

boxplot(
  all.mean[,7] / all.mean[,2], all.mean[,8] / all.mean[,3], 
  all.mean[,9] / all.mean[,4], all.mean[,10] / all.mean[,5],
  names=c("MAXENT standard", "MAXENT TGBA", "BRT standard", "BRT TGBA"),
  title="Mean HSV inside/total")

boxplot(
  all.mean[,12] / all.mean[,2], all.mean[,13] / all.mean[,3], 
  all.mean[,14] / all.mean[,4], all.mean[,15] / all.mean[,5],
  names=c("MAXENT standard", "MAXENT TGBA", "BRT standard", "BRT TGBA"),
  title="Mean HSV outside/total")

  boxplot(all.mean[,c("mi.MAX", "mi.MAX.b", "mo.MAX", "mo.MAX.b", "vi.MAX", "vi.MAX.b", "vo.MAX", "vo.MAX.b")])


#test difference
t.test(all.mean[,c("mi.MAX")], all.mean[,c("mi.MAX.b")], paired=F)
t.test(all.mean[,c("mo.MAX")], all.mean[,c("mo.MAX.b")], paired=F)
t.test(all.mean[,c("vi.MAX")], all.mean[,c("vi.MAX.b")], paired=F)
t.test(all.mean[,c("vo.MAX")], all.mean[,c("vo.MAX.b")], paired=F)

boxplot(all.mean[,c("mi.MAX", "mi.MAX.b", "mo.MAX", "mo.MAX.b")], xlab="intersect category", ylab="predicted habitat suitability", names=c("standard inside", "TGBG inside", "standard outside", "TGBG outside"))

#mean for maxent and maxent.b
mean.s <- mean(all.mean[,c("mi.MAX")] + all.mean[,c("mo.MAX")])
mean.b <- boxplot(all.mean[,c("mi.MAX", "mi.MAX.b", "mo.MAX", "mo.MAX.b")])


boxplot(log(mean.out/mean.in))
boxplot(log(var.out/var.in))


#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#
# plot models in environmental space ----

  data.envi <- stack("data/data.envi.grd")

  # Read model output
  species.list <- c(600143, 601342, 605367, 600656, 600138)
  species.rasters <- data.envi
  for(s in species.list){
  load(paste("C:/data/modelling/output/baseline/SDMeval_fulld_", s, ".RData", sep=""))
  spec.model.normal <- results
  load(paste("C:/data/modelling/output/bias1/SDMeval_fulld_sb1_", s, ".RData", sep=""))
  spec.model.bias <- results
  
  print(s)
  rasterpred <- stack(spec.model.normal[[2]][[2]], spec.model.bias[[2]][[2]])
  names(rasterpred) <- c(paste(s, ".standard", sep=""), paste(s, ".TGBA", sep=""))
  species.rasters <- stack(species.rasters, rasterpred)
  }

    #600143=Thunnus albacares
    #601342=Pleuronectes platessa
    #605367=Epinephelus areolatus
    #600656=Centroscyllium fabricii
    #600138=Somniosus microcephalus

  data.values <- values(species.rasters)
  data.values <- data.values[order(data.values[,8]), ]
  
  # Plot values in environmental space
  par(mfrow=c(2,3))
  rbPal <- colorRampPalette(c("grey", "black"))
  s1="X600143.standard"
  s2="X600143.TGBA"
  data.values <- data.values[order(data.values[,8]), ]
  plot(data.values[,c("landdist", "sstanmean")], pch=20, cex=data.values[,s1]*2, col=rbPal(10)[as.numeric(cut(data.values[,s1], breaks=10))])
  plot(data.values[,c("sstanmean", "depthmean")], pch=20, cex=data.values[,s1]*2, col=rbPal(10)[as.numeric(cut(data.values[,s1], breaks=10))])
  plot(data.values[,c("sstanmean", "primprod")], pch=20, cex=data.values[,s1]*2, col=rbPal(10)[as.numeric(cut(data.values[,s1], breaks=10))])
  plot(data.values[,c("landdist", "sstanmean")], pch=20, cex=data.values[,s2]*2, col=rbPal(10)[as.numeric(cut(data.values[,s2], breaks=10))])
  plot(data.values[,c("sstanmean", "depthmean")], pch=20, cex=data.values[,s2]*2, col=rbPal(10)[as.numeric(cut(data.values[,s2], breaks=10))])
  plot(data.values[,c("sstanmean", "primprod")], pch=20, cex=data.values[,s2]*2, col=rbPal(10)[as.numeric(cut(data.values[,s2], breaks=10))])
  
  s1="X601342.standard"
  s2="X601342.TGBA"
  data.values <- data.values[order(data.values[,s1]), ]
  plot(data.values[,c("landdist", "sstanmean")], pch=20, cex=data.values[,s1]*2, col=rbPal(10)[as.numeric(cut(data.values[,s1], breaks=10))])
  plot(data.values[,c("sstanmean", "depthmean")], pch=20, cex=data.values[,s1]*2, col=rbPal(10)[as.numeric(cut(data.values[,s1], breaks=10))])
  plot(data.values[,c("sstanmean", "primprod")], pch=20, cex=data.values[,s1]*2, col=rbPal(10)[as.numeric(cut(data.values[,s1], breaks=10))])
  plot(data.values[,c("landdist", "sstanmean")], pch=20, cex=data.values[,s2]*2, col=rbPal(10)[as.numeric(cut(data.values[,s2], breaks=10))])
  plot(data.values[,c("sstanmean", "depthmean")], pch=20, cex=data.values[,s2]*2, col=rbPal(10)[as.numeric(cut(data.values[,s2], breaks=10))])
  plot(data.values[,c("sstanmean", "primprod")], pch=20, cex=data.values[,s2]*2, col=rbPal(10)[as.numeric(cut(data.values[,s2], breaks=10))])

  s1="X605367.standard"
  s2="X605367.TGBA"
  data.values <- data.values[order(data.values[,s1]), ]
  plot(data.values[,c("landdist", "sstanmean")], pch=20, cex=data.values[,s1]*2, col=rbPal(10)[as.numeric(cut(data.values[,s1], breaks=10))])
  plot(data.values[,c("sstanmean", "depthmean")], pch=20, cex=data.values[,s1]*2, col=rbPal(10)[as.numeric(cut(data.values[,s1], breaks=10))])
  plot(data.values[,c("sstanmean", "primprod")], pch=20, cex=data.values[,s1]*2, col=rbPal(10)[as.numeric(cut(data.values[,s1], breaks=10))])
  plot(data.values[,c("landdist", "sstanmean")], pch=20, cex=data.values[,s2]*2, col=rbPal(10)[as.numeric(cut(data.values[,s2], breaks=10))])
  plot(data.values[,c("sstanmean", "depthmean")], pch=20, cex=data.values[,s2]*2, col=rbPal(10)[as.numeric(cut(data.values[,s2], breaks=10))])
  plot(data.values[,c("sstanmean", "primprod")], pch=20, cex=data.values[,s2]*2, col=rbPal(10)[as.numeric(cut(data.values[,s2], breaks=10))])

  s1="X600656.standard"
  s2="X600656.TGBA"
  data.values <- data.values[order(data.values[,s1]), ]
  plot(data.values[,c("landdist", "sstanmean")], pch=20, cex=data.values[,s1]*2, col=rbPal(10)[as.numeric(cut(data.values[,s1], breaks=10))])
  plot(data.values[,c("sstanmean", "depthmean")], pch=20, cex=data.values[,s1]*2, col=rbPal(10)[as.numeric(cut(data.values[,s1], breaks=10))])
  plot(data.values[,c("sstanmean", "primprod")], pch=20, cex=data.values[,s1]*2, col=rbPal(10)[as.numeric(cut(data.values[,s1], breaks=10))])
  plot(data.values[,c("landdist", "sstanmean")], pch=20, cex=data.values[,s2]*2, col=rbPal(10)[as.numeric(cut(data.values[,s2], breaks=10))])
  plot(data.values[,c("sstanmean", "depthmean")], pch=20, cex=data.values[,s2]*2, col=rbPal(10)[as.numeric(cut(data.values[,s2], breaks=10))])
  plot(data.values[,c("sstanmean", "primprod")], pch=20, cex=data.values[,s2]*2, col=rbPal(10)[as.numeric(cut(data.values[,s2], breaks=10))])

  s1="X600138.standard"
  s2="X600138.TGBA"
  data.values <- data.values[order(data.values[,s1]), ]
  plot(data.values[,c("landdist", "sstanmean")], pch=20, cex=data.values[,s1]*2, col=rbPal(10)[as.numeric(cut(data.values[,s1], breaks=10))])
  plot(data.values[,c("sstanmean", "depthmean")], pch=20, cex=data.values[,s1]*2, col=rbPal(10)[as.numeric(cut(data.values[,s1], breaks=10))])
  plot(data.values[,c("sstanmean", "iceconann")], pch=20, cex=data.values[,s1]*2, col=rbPal(10)[as.numeric(cut(data.values[,s1], breaks=10))])
  plot(data.values[,c("landdist", "sstanmean")], pch=20, cex=data.values[,s2]*2, col=rbPal(10)[as.numeric(cut(data.values[,s2], breaks=10))])
  plot(data.values[,c("sstanmean", "depthmean")], pch=20, cex=data.values[,s2]*2, col=rbPal(10)[as.numeric(cut(data.values[,s2], breaks=10))])
  plot(data.values[,c("sstanmean", "iceconann")], pch=20, cex=data.values[,s2]*2, col=rbPal(10)[as.numeric(cut(data.values[,s2], breaks=10))])

#
  
  ###########################
  ###########################
  
  
  
  

  #MODIFIED EVAL FUNCTION ####
  eval2 <- function (p, a, model, x, tr, ...){
    if (!missing(x)) {
      p <- predict(model, data.frame(extract(x, p)))
      a <- predict(model, data.frame(extract(x, a)), ...)
    }
    else if (is.vector(p) & is.vector(a)) {
    }
    else {
      p <- predict(model, data.frame(p), ...)
      a <- predict(model, data.frame(a), ...)
    }
    
    p <- na.omit(p)
    a <- na.omit(a)
    np <- length(p)
    na <- length(a)
    if (na == 0 | np == 0){
      stop("cannot evaluate a model without absence and presence data that are not NA")}
    if (missing(tr)){
      if (length(p) > 1000){
        tr <- as.vector(quantile(p, 0:1000/1000))
      }
      else {
        tr <- p}
      if (length(a) > 1000){
        tr <- c(tr, as.vector(quantile(a, 0:1000/1000)))
      }
      else {
        tr <- c(tr, a)}
      tr <- sort(unique(round(tr, 8)))
      tr <- c(tr - 1e-04, tr[length(tr)] + c(0, 1e-04))
    }
    else {
      tr <- sort(as.vector(tr))}
    
    N <- na + np
    xc <- new("ModelEvaluation")
    xc@presence = p
    xc@absence = a
    R <- sum(rank(c(p, a))[1:np]) - (np * (np + 1)/2)
    xc@auc <- R/(na * np)
    cr <- try(cor.test(c(p, a), c(rep(1, length(p)), rep(0, length(a)))), 
              silent = TRUE)
    if (class(cr) != "try-error") {
      xc@cor <- cr$estimate
      xc@pcor <- cr$p.value
    }
    res <- matrix(ncol = 4, nrow = length(tr))
    colnames(res) <- c("tp", "fp", "fn", "tn")
    xc@t <- tr
    for (i in 1:length(tr)) {
      res[i, 1] <- length(p[p >= tr[i]])
      res[i, 2] <- length(a[a >= tr[i]])
      res[i, 3] <- length(p[p < tr[i]])
      res[i, 4] <- length(a[a < tr[i]])
    }
    xc@confusion = res
    a = res[, 1]
    b = res[, 2]
    c = res[, 3]
    d = res[, 4]
    xc@np <- as.integer(np)
    xc@na <- as.integer(na)
    xc@prevalence = (a + c)/N
    xc@ODP = (b + d)/N
    xc@CCR = (a + d)/N
    xc@TPR = a/(a + c)
    xc@TNR = d/(b + d)
    xc@FPR = b/(b + d)
    xc@FNR = c/(a + c)
    xc@PPP = a/(a + b)
    xc@NPP = d/(c + d)
    xc@MCR = (b + c)/N
    xc@OR = (a * d)/(c * b)
    prA = (a + d)/N
    prY = (a + b)/N * (a + c)/N
    prN = (c + d)/N * (b + d)/N
    prE = prY + prN
    xc@kappa = (prA - prE)/(1 - prE)
    return(xc)
  }

