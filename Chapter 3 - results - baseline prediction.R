library(sp)
library(raster)
library(dismo)
library(rgdal)

#function to read results from Rfile
f.load.results <- function(TaxonKey, fileprefix="SDMeval_fulld_", filedir="C:/data/modelling/output/baseline/"){
  try(load(paste(filedir, fileprefix, TaxonKey, ".RData", sep="")))
  if(exists("results")==F){
    print(c(species.list$TaxonKey[s], as.character(species.list$TaxonName[s]), "not found"))
  }else if(exists("results")==T){
    preds <- stack(results[[2]])
    names(preds) <- c("BIOCLIM", "MAXENT", "BRT")
    evals <- results[[3]]   
    return(list(preds, evals))
  }else{
    print(c(species.list$TaxonKey[s], as.character(species.list$TaxonName[s]), "ERROR!!!"))
  }}


#read data species list
species.list <- read.csv("C:/data/modelling/data/speciesinfo.csv", header=T)
#setwd("C:/data/modelling/output/baseline")
setwd("C:/data/modelling/output/bias1")

#filedir="C:/data/modelling/output/baseline/"
filedir="C:/data/modelling/output/bias1/"

#identify species modelled
#files <- unlist(strsplit(unlist(strsplit(list.files()[grep("SDMeval_fulld_", list.files())], "SDMeval_fulld_")),".RData"))
files <- unlist(strsplit(unlist(strsplit(list.files()[grep("SDMeval_fulld_", list.files())], "SDMeval_fulld_sb1_")),".RData"))
species.list <- species.list[species.list$TaxonKey %in% files,]

#loop over species
s=1
preds=NULL
evals=NULL
results=NULL
models <- c("BIOCLIM", "MAXENT", "BRT")

for(s in 1:length(species.list$TaxonKey)){
  try(loaded <- f.load.results(TaxonKey=species.list$TaxonKey[s], fileprefix="SDMeval_fulld_sb1_", filedir="C:/data/modelling/output/bias1/"))
  if(exists("loaded")){
    thresholds <- lapply(loaded[[2]], function(x) threshold(x))
    thresholds <- cbind(thresholds[[1]], thresholds[[2]], thresholds[[3]])
    names(thresholds) <- c("BIOCLIM.kappa", "BIOCLIM.spec_sens", "BIOCLIM.no_omission", "BIOCLIM.prevalence", "BIOCLIM.equal_sens_spec", "BIOCLIM.sensitivity",
                           "MAXENT.kappa", "MAXENT.spec_sens", "MAXENT.no_omission", "MAXENT.prevalence", "MAXENT.equal_sens_spec", "MAXENT.sensitivity",
                           "BRT.kappa", "BRT.spec_sens", "BRT.no_omission", "BRT.prevalence", "BRT.equal_sens_spec", "BRT.sensitivity")
    try(values <- cbind(TaxonKey=species.list$TaxonKey[s],
                       auc.BIOCLIM=loaded[[2]][[1]]@auc, auc.MAXENT=loaded[[2]][[2]]@auc, auc.BRT=loaded[[2]][[3]]@auc,
                       cor.BIOCLIM=loaded[[2]][[1]]@cor, cor.MAXENT=loaded[[2]][[2]]@cor, cor.BRT=loaded[[2]][[3]]@cor,
                       thresholds))
    results <- rbind(results, values)
    print(paste(species.list$TaxonKey[s], "done"))
    rm(loaded)
  }else{
    print(paste("###", species.list$TaxonKey[s], "not found"))
  }}

results.bl <- results
#results.sb <- results

#Change names to bind to another results object
colnames(results.sb) <- unlist(lapply(colnames(results.sb), function(x) paste("sb.", x, sep="")))
results.sb <- results.sb[results.sb$sb.TaxonKey %in% results.bl$TaxonKey,]
comb <- cbind(results.bl, results.sb)
save(comb, file="output/results.bl.sb.RData")

#Plot some graphs ----
boxplot(comb[,c("auc.BIOCLIM", "auc.MAXENT", "sb.auc.MAXENT", "auc.BRT", "sb.auc.BRT")])
boxplot(comb[,c("cor.BIOCLIM", "cor.MAXENT", "sb.cor.MAXENT", "cor.BRT", "sb.cor.BRT")])

#test significance
t.test(comb$sb.auc.BRT, comb$sb.auc.MAXENT, paired=T, var.equal=F, na.action="na.omit", alternative="less")

var.test(results$auc.BIOCLIM, results$auc.MAXENT)
var.test(results$auc.BRT, results$auc.MAXENT)
var.test(results$auc.BIOCLIM, results$auc.BRT)
t.test(results$auc.BIOCLIM, results$auc.MAXENT, paired=T, var.equal=F, na.action="na.omit", alternative="less")
t.test(results$auc.BRT, results$auc.MAXENT, paired=T, var.equal=F, na.action="na.omit", alternative="less")
t.test(results$auc.BIOCLIM, results$auc.BRT, paired=T, var.equal=F, na.action="na.omit", alternative="less")

var.test(results$cor.BIOCLIM, results$cor.MAXENT)
var.test(results$cor.BRT, results$cor.MAXENT)
var.test(results$cor.BIOCLIM, results$cor.BRT)
t.test(results$cor.BIOCLIM, results$cor.MAXENT, paired=T, var.equal=F, na.action="na.omit", alternative="less")
t.test(results$cor.BRT, results$cor.MAXENT, paired=T, var.equal=F, na.action="na.omit", alternative="less")
t.test(results$cor.BIOCLIM, results$cor.BRT, paired=T, var.equal=F, na.action="na.omit", alternative="less")




# Compare to IUCN range maps ----
  #Get land shapefile and extentraster to calculate coverage of land area per cell
  #land <- readOGR("data_other/shapefiles", "land")
  #rasterize(land, data.extentraster, )
  
  # Read IUCN range maps
  setwd("C:/data/modelling")
  data.extent <- raster("data/data.extent.grd")
  shapefiles <- list.files("data/iucn/SAUP/singleshapefiles", pattern=c(".shp$"), full.names=F)
  for(s in 1:length(shapefiles)){
    TaxonKey <- as.numeric(strsplit(shapefiles[s],"[.]")[[1]][1])
    spec.poly <- readOGR("data/iucn/SAUP/singleshapefiles", TaxonKey)
    spec.raster <- rasterize(spec.poly, data.extent, background=0, getCover=T)
    spec.raster <- mask(spec.raster, data.extent)
    assign(paste("s", TaxonKey, sep=""), spec.raster)
    writeRaster(get(paste("s", TaxonKey, sep="")), file=paste("data/iucn/raster/", TaxonKey, sep=""), format="raster", overwrite=T)
    
    if(TaxonKey %in% species.list$TaxonKey){
      #try(loaded <- f.load.results(TaxonKey=species.list$TaxonKey[s], fileprefix="SDMeval_fulld_", filedir="C:/data/modelling/output/baseline/"))
      png(paste(TaxonKey, ".poly.png", sep=""), width=10350, height=5400, units="px")
      plot(data.extent)
      plot(spec.poly, add=T, col="red", border="black")
      dev.off()
    }
  }


#Read IUCN rasters
library(dismo)
setwd("C:/data/modelling")
load("data/data.species.RData")
data.envi <- stack("data/data.envi.grd")
data.extent <- raster("data/data.extent.grd")
#data.faoraster <- raster("data/data.faoraster.grd")

files <- unlist(strsplit(list.files("data/iucn/raster", pattern=c(".grd$"), full.names=F), split=".grd"))
files <- files[files %in% comb$TaxonKey]
a <- data.species[data.species$set=="background", c("x", "y")]
for(thisfile in files){
  try(loaded1 <- f.load.results(TaxonKey=as.numeric(thisfile), fileprefix="SDMeval_fulld_", filedir="C:/data/modelling/output/baseline/"))
  try(loaded2 <- f.load.results(TaxonKey=as.numeric(thisfile), fileprefix="SDMeval_fulld_sb1_", filedir="C:/data/modelling/output/bias1/"))
  data.range <- raster(paste("data/iucn/raster/", thisfile, ".grd", sep=""))
  
  par(mfrow=c(3,1))
  plot(loaded1[[1]][[2]])
  plot(loaded2[[1]][[2]])
  plot(data.range)
  
  par(mfrow=c(1,1))
  bla1 <- values(mask(loaded1[[1]][[1]], data.range, maskvalue=0, inver=F))
  bla2 <- values(mask(loaded2[[1]][[1]], data.range, maskvalue=0, inver=T))
  bla3 <- values(mask(loaded1[[1]][[2]], data.range, maskvalue=0, inver=F))
  bla4 <- values(mask(loaded2[[1]][[2]], data.range, maskvalue=0, inver=F))
  bla5 <- values(mask(loaded1[[1]][[2]], data.range, maskvalue=0, inver=T))
  bla6 <- values(mask(loaded2[[1]][[2]], data.range, maskvalue=0, inver=T))
  boxplot(bla1, bla2, bla3, bla4, names=c("normal", "bias", "normal", "bias"))
  boxplot(bla3, bla4, bla5, bla6)
  
  plot(values(loaded2[[1]][[2]]), values(data.range))
  
  #p <- data.species[data.species$TaxonKey==thisfile & data.species$set=="test", c("x", "y")]
  #evaluate(p, a, data.range)
