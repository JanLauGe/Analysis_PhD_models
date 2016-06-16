# LIBRARIES ####
library(sp)
library(raster)
library(rgdal)
library(dismo)

#################
##VARIABLES #####
#*#*#*#*#*#*#*#*#

  setwd("C:/data/modelling")
  data.extent <- raster("C:/data/modelling/data/data.extent.grd")
  species.list1 <- unlist(strsplit(unlist(strsplit(list.files("C:/data/modelling/output/geodist", pattern=".RData$"), split="SDMeval_geodist_")), split=".RData"))
  species.list2 <- unlist(strsplit(unlist(strsplit(list.files("C:/data/modelling/output/baseline", pattern="*SDMeval.clipped*"), split="SDMeval.clipped_")), split=".RData"))
  species.list.s <- species.list1[species.list1 %in% species.list2]

  load("data/data.species.RData")
  data.sample <- data.species[data.species$set=="background",c("x", "y")]


#*#*#*#*#*#*#*#*#*#*#*#
# Model comparison ####
#*#*#*#*#*#*#*#*#*#*#*#

  #Compare models run with geodistance to baseline modells and baseline models clipped by FAO area

  s=1
  all.eval=NULL
  for(s in 1:length(species.list.s)){
    
    print(species.list.s[s])
    
  # Read model output
    load(paste("C:/data/modelling/output/baseline/SDMeval_fulld_", species.list.s[s], ".RData", sep=""))
    spec.model.baseline <- results
    load(paste("C:/data/modelling/output/baseline/SDMeval.clipped_", species.list.s[s], ".RData", sep=""))
    spec.model.clipped <- results
    load(paste("C:/data/modelling/output/geodist/SDMeval_geodist_", species.list.s[s], ".RData", sep=""))
    spec.model.geodist <- results

    
    s.eval <- c("BIO.bas.auc"=spec.model.baseline[[3]][[1]]@auc,
                "MAX.bas.auc"=spec.model.baseline[[3]][[2]]@auc,
                "BRT.bas.auc"=spec.model.baseline[[3]][[3]]@auc,
                "BIO.cli.auc"=spec.model.clipped[[3]][[1]]@auc,
                "MAX.cli.auc"=spec.model.clipped[[3]][[2]]@auc,
                "BRT.cli.auc"=spec.model.clipped[[3]][[3]]@auc,
                "BIO.geo.auc"=spec.model.geodist[[3]][[1]]@auc,
                "MAX.geo.auc"=spec.model.geodist[[3]][[2]]@auc,
                "BRT.geo.auc"=spec.model.geodist[[3]][[3]]@auc,
                "BIO.bas"=spec.model.baseline[[3]][[1]]@cor,
                "MAX.bas"=spec.model.baseline[[3]][[2]]@cor,
                "BRT.bas"=spec.model.baseline[[3]][[3]]@cor,
                "BIO.cli"=spec.model.clipped[[3]][[1]]@cor,
                "MAX.cli"=spec.model.clipped[[3]][[2]]@cor,
                "BRT.cli"=spec.model.clipped[[3]][[3]]@cor,
                "BIO.geo"=spec.model.geodist[[3]][[1]]@cor,
                "MAX.geo"=spec.model.geodist[[3]][[2]]@cor,
                "BRT.geo"=spec.model.geodist[[3]][[3]]@cor)
    all.eval <- rbind(all.eval, s.eval)
    rownames(all.eval)[s] <- species.list.s[s]
  }

  colnames(all.eval)
  write.table(all.eval, file="compare geodist model to normal and clipped maps - standard auc and cor.txt", sep="\t")  
  all.eval1 <- read.table(file="compare geodist model to normal and clipped maps - standard auc and cor.txt", sep="\t")  

all.eval1 <- all.eval
boxplot(all.eval1)
boxplot(all.eval1[,c(1,4,7,2,5,8,3,6,9)])
boxplot(all.eval1[,c(10,13,16,11,14,17,12,15,18)])
t.test(all.eval1[,5], all.eval1[,8], paired=T)
t.test(all.eval1[,6], all.eval1[,9], paired=T)


# #*#*#*#*#*#*#*#*#*#*#
# # IUCN RANGE MAP ####
# #*#*#*#*#*#*#*#*#*#*#
# 
# #Compare model outputs to IUCN and SAUP range maps
# s=1
# all.eval2=NULL
# for(s in 1:length(species.list.s)){
#   print(species.list.s[s])
#   spec.raster <- raster(paste("data/speciesraster/", species.list.s[s], ".grd", sep=""))
#   #create presence and absence points from expert range map for evaluation ----
#   species.sample <- cbind(data.sample, exp=extract(spec.raster, data.sample))
#   species.sample <- species.sample[complete.cases(species.sample),]
#   species.sample <- SpatialPointsDataFrame(species.sample[,1:2], species.sample)
#   p <- species.sample[species.sample@data[,3]>0,]
#   a <- species.sample[species.sample@data[,3]==0,]
# 
#   if(length(p)==0 | length(a)==0){
#     print("skipped")
#     next}
#   
#   # Read model output
#   load(paste("C:/data/modelling/output/baseline/SDMeval_fulld_", species.list.s[s], ".RData", sep=""))
#   spec.model.baseline <- results
#   load(paste("C:/data/modelling/output/baseline/SDMeval.clipped_", species.list.s[s], ".RData", sep=""))
#   spec.model.clipped <- results
#   load(paste("C:/data/modelling/output/geodist/SDMeval_geodist_", species.list.s[s], ".RData", sep=""))
#   spec.model.geodist <- results
#   
#   #Evaluation----
#   BIO.bas <- evaluate(extract(spec.model.baseline[[2]][[1]], p), extract(spec.model.baseline[[2]][[1]], a))
#   MAX.bas <- evaluate(extract(spec.model.baseline[[2]][[2]], p), extract(spec.model.baseline[[2]][[2]], a))
#   BRT.bas <- evaluate(extract(spec.model.baseline[[2]][[3]], p), extract(spec.model.baseline[[2]][[3]], a))
#   BIO.cli <- evaluate(extract(spec.model.clipped[[2]][[1]], p), extract(spec.model.clipped[[2]][[1]], a))
#   MAX.cli <- evaluate(extract(spec.model.clipped[[2]][[2]], p), extract(spec.model.clipped[[2]][[2]], a))
#   BRT.cli <- evaluate(extract(spec.model.clipped[[2]][[3]], p), extract(spec.model.clipped[[2]][[3]], a))
#   BIO.geo <- evaluate(extract(spec.model.geodist[[2]][[1]], p), extract(spec.model.geodist[[2]][[1]], a))
#   MAX.geo <- evaluate(extract(spec.model.geodist[[2]][[2]], p), extract(spec.model.geodist[[2]][[2]], a))
#   BRT.geo <- evaluate(extract(spec.model.geodist[[2]][[3]], p), extract(spec.model.geodist[[2]][[3]], a))
#   
#   s.eval <- c("BIO.bas.auc"=BIO.bas@auc,
#               "MAX.bas.auc"=MAX.bas@auc,
#               "BRT.bas.auc"=BRT.bas@auc,
#               "BIO.cli.auc"=BIO.cli@auc,
#               "MAX.cli.auc"=MAX.cli@auc,
#               "BRT.cli.auc"=BRT.cli@auc,
#               "BIO.geo.auc"=BIO.geo@auc,
#               "MAX.geo.auc"=MAX.geo@auc,
#               "BRT.geo.auc"=BRT.geo@auc,
#               "BIO.bas.cor"=BIO.bas@cor,
#               "MAX.bas.cor"=MAX.bas@cor,
#               "BRT.bas.cor"=BRT.bas@cor,
#               "BIO.cli.cor"=BIO.cli@cor,
#               "MAX.cli.cor"=MAX.cli@cor,
#               "BRT.cli.cor"=BRT.cli@cor,
#               "BIO.geo.cor"=BIO.geo@cor,
#               "MAX.geo.cor"=MAX.geo@cor,
#               "BRT.geo.cor"=BRT.geo@cor)
#   
#   all.eval2 <- rbind(all.eval2, s.eval)
#   rownames(all.eval2)[length(all.eval2[,1])] <- species.list.s[s]
# }
#   
# write.table(all.eval2, file="compare geodist model to normal and clipped maps - expert auc and cor.txt", sep="\t")  
# all.eval1 <-  read.table(file="compare geodist model to normal and clipped maps - standard auc and cor.txt", sep="\t")
# all.eval2 <- read.table(file="compare geodist model to normal and clipped maps - expert auc and cor.txt", sep="\t")
# 
# boxplot(all.eval1[,c(2,5,8,3,6,9)])
# boxplot(all.eval1[,c(11,14,17,12,15,18)])
# boxplot(all.eval1[,c(10,13,16,11,14,17,12,15,18)], names=c("normal", "clipped", "distance", "normal", "clipped", "distance", "normal", "clipped", "distance"))
# 
# t.test(all.eval1[,5], all.eval1[,8], paired=T)
# t.test(all.eval1[,6], all.eval1[,9], paired=T)
# t.test(all.eval1[,14], all.eval1[,17], paired=T)
# t.test(all.eval1[,15], all.eval1[,18], paired=T)
# 
# boxplot(all.eval2)
# boxplot(all.eval2[,c(2,5,8,3,6,9)])
# boxplot(all.eval2[,c(11,14,17,12,15,18)])
# t.test(all.eval2[,5], all.eval2[,8], paired=T)
# t.test(all.eval2[,6], all.eval2[,9], paired=T)
# 
# 

a <- as.data.frame(all.eval1)
#new table
bla <- rbind(data.frame("TaxonKey"=as.numeric(rownames(a)), "auc"=as.numeric(a$BIO.bas.auc), "cor"=as.numeric(a$BIO.bas.cor), "model"=as.factor("BIOCLIM"), "run"=as.factor("baseline")),
             data.frame("TaxonKey"=as.numeric(rownames(a)), "auc"=as.numeric(a$BIO.cli.auc), "cor"=as.numeric(a$BIO.cli.cor), "model"=as.factor("BIOCLIM"), "run"=as.factor("clipped")),
             data.frame("TaxonKey"=as.numeric(rownames(a)), "auc"=as.numeric(a$BIO.geo.auc), "cor"=as.numeric(a$BIO.geo.cor), "model"=as.factor("BIOCLIM"), "run"=as.factor("expert")),
             data.frame("TaxonKey"=as.numeric(rownames(a)), "auc"=as.numeric(a$MAX.bas.auc), "cor"=as.numeric(a$MAX.bas.cor), "model"=as.factor("MAXENT"), "run"=as.factor("baseline")),
             data.frame("TaxonKey"=as.numeric(rownames(a)), "auc"=as.numeric(a$MAX.cli.auc), "cor"=as.numeric(a$MAX.cli.cor), "model"=as.factor("MAXENT"), "run"=as.factor("clipped")),
             data.frame("TaxonKey"=as.numeric(rownames(a)), "auc"=as.numeric(a$MAX.geo.auc), "cor"=as.numeric(a$MAX.geo.cor), "model"=as.factor("MAXENT"), "run"=as.factor("expert")),
             data.frame("TaxonKey"=as.numeric(rownames(a)), "auc"=as.numeric(a$BRT.bas.auc), "cor"=as.numeric(a$BRT.bas.cor), "model"=as.factor("BRT"), "run"=as.factor("baseline")),
             data.frame("TaxonKey"=as.numeric(rownames(a)), "auc"=as.numeric(a$BRT.cli.auc), "cor"=as.numeric(a$BRT.cli.cor), "model"=as.factor("BRT"), "run"=as.factor("clipped")),
             data.frame("TaxonKey"=as.numeric(rownames(a)), "auc"=as.numeric(a$BRT.geo.auc), "cor"=as.numeric(a$BRT.geo.cor), "model"=as.factor("BRT"), "run"=as.factor("expert")))

library(ggplot2)
library(devtools)
library(easyGgplot2)

ggplot2.boxplot(data=bla, xName="run", yName="auc", groupName="model",
                faceting=TRUE, facetingVarNames="model", facetingDirection="horizontal", facetingFont=c(12, "bold", "black"),
                width=1, notch=T, outlier.shape=21, groupColors=c('#999999','#E69F00','#56B4E9'), showLegend=FALSE,
                #axisLine=c(0.5, "solid", "black"),  backgroundColor="white", removePanelGrid=TRUE,removePanelBorder=TRUE,
                xtitle="Modelrun", ytitle="AUC") + 
                theme(axis.text.x=element_text(angle=-45, hjust=0, color="gray40"), axis.text.y=element_text(color="gray40"), axis.title.x=element_text(vjust=-0.5, hjust=0.5), axis.title.y=element_text(vjust=1.5)) + expand_limits(y=c(0.6,1))

ggplot2.boxplot(data=bla, xName="run", yName="cor", groupName="model",
                faceting=TRUE, facetingVarNames="model", facetingDirection="horizontal", facetingFont=c(12, "bold", "black"),
                width=1, notch=T, outlier.shape=21, groupColors=c('#999999','#E69F00','#56B4E9'), showLegend=FALSE,
                #axisLine=c(0.5, "solid", "black"),  backgroundColor="white", removePanelGrid=TRUE,removePanelBorder=TRUE,
                xtitle="Modelrun", ytitle="COR") + 
  theme(axis.text.x=element_text(angle=-45, hjust=0, color="gray40"), axis.text.y=element_text(color="gray40"), axis.title.x=element_text(vjust=-0.5, hjust=0.5), axis.title.y=element_text(vjust=1.5)) + expand_limits(y=c(0.6,1))

