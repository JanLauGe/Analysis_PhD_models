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
setwd("G:/HPC/output/aic")

#filedir="C:/data/modelling/output/baseline/"
filedir="G:/HPC/output/aic/"

#identify species modelled
#files <- unlist(strsplit(unlist(strsplit(list.files()[grep("SDMeval_fulld_", list.files())], "SDMeval_fulld_")),".RData"))
#files <- unlist(strsplit(unlist(strsplit(list.files()[grep("SDMeval_fulld_", list.files())], "SDMeval_fulld_sb1_")),".RData"))
files <- unlist(strsplit(unlist(strsplit(list.files()[grep("SDMeval_fulld_aic_", list.files())], "SDMeval_fulld_aic_")),".RData"))
species.list <- species.list[species.list$TaxonKey %in% files,]

#loop over species
s=1
preds=NULL
evals=NULL
results=NULL
models = c("BIOCLIM", "MAXENT", "BRT")
fileprefix = "SDMeval_fulld_aic_"
filedir="G:/HPC/output/aic/"
for(s in 1:length(species.list$TaxonKey)){
  try(loaded <- f.load.results(TaxonKey=species.list$TaxonKey[s], fileprefix=fileprefix, filedir=filedir))
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

results.aic <- results
#results.sb <- results

#Change names to bind to another results object
colnames(results.aic) <- unlist(lapply(colnames(results.aic), function(x) paste("aic.", x, sep="")))
write.csv(results.aic, file="C:/data/modelling/output/results_aic_tuning.csv")
results.aic <- read.csv("C:/data/modelling/output/results_aic_tuning.csv")

file.list <- list.files("G:/HPC/output/aic/", pattern="AIC_summary_")
results <- data.frame(NULL)
for (file in file.list){
  result <- read.csv(paste("G:/HPC/output/aic/", file, sep=""))
  results <- rbind(results, result)
}
results.aic2 <- results
write.csv(results.aic2, file="C:/data/modelling/output/results_aic_tuning2.csv")


######
results.aic1 <- read.csv("C:/data/modelling/output/results_aic_tuning.csv")
results.aic2 <- read.csv("C:/data/modelling/output/results_aic_tuning2.csv")

#Merge with other results
specinfo <- read.csv("C:/Users/LaurensG/Desktop/specinfo.csv")
merge(results,)


