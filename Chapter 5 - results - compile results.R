# LIBRARIES ####
library(sp)
library(raster)
library(rgdal)
library(dismo)

# DATA ####
  setwd("C:/data/modelling")
  data.extent <- raster("C:/data/modelling/data/data.extent.grd")
  species.list <- read.csv("C:/data/modelling/data/species.list.csv")

path="C:/data/modelling/output/baseline"
filename="SDMeval_fulld_"

# Used to read results from R objects and combine into table
merge.results <- function(path, filename){
  #path should be full filepath of result files, e.g. "C:/data/modelling/output/baseline"
  #filename should be name structure for files given without TaxonKey and .RData, e.g. "SDMeval_fulld_"
  list.files <- unlist(strsplit(unlist(strsplit(list.files(path, pattern=paste("*",filename,"*",sep="")), split=filename)) , split=".RData"))
  eval.modelruns=NULL
  for(s in 1:length(list.files)){
    load(paste(path, "/", filename, list.files[s], ".RData", sep=""))
    s.model <- results
    s.eval <- data.frame("TaxonKey"=list.files[s],
                         "BIO.auc"=results[[3]][[1]]@auc,
                         "MAX.auc"=results[[3]][[2]]@auc,
                         "BRT.auc"=results[[3]][[3]]@auc,
                         "BIO.cor"=results[[3]][[1]]@cor,
                         "MAX.cor"=results[[3]][[2]]@cor,
                         "BRT.cor"=results[[3]][[3]]@cor)
    
    eval.modelruns <- rbind(eval.modelruns, s.eval)
    print(paste(list.files[s], "done"))}
  return(eval.modelruns)
  }

#data.plot <- read.csv("C:/Users/LaurensG/Desktop/completecomparison_SAUP.csv")


#baseline fulld
base.f <- merge.results("C:/data/modelling/output/baseline", "SDMeval_fulld_")
rownames(base.f) <- base.f$TaxonKey
colnames(base.f)[2:7] <- paste("base.f.", colnames(base.f)[2:7], sep="")

#bias fulld
bias.f <- merge.results("C:/data/modelling/output/bias1", "SDMeval_fulld_sb1_")
rownames(bias.f) <- bias.f$TaxonKey
colnames(bias.f)[2:7] <- paste("bias.f.", colnames(bias.f)[2:7], sep="")

#baseline clip
base.c <- merge.results("C:/data/modelling/output/baseline", "SDMeval.clipped_")
rownames(base.c) <- base.c$TaxonKey
colnames(base.c)[2:7] <- paste("base.c.", colnames(base.c)[2:7], sep="")

#geodist
geod.c <- merge.results("C:/data/modelling/output/geodist", "SDMeval_geodist_")
rownames(geod.c) <- geod.c$TaxonKey
colnames(geod.c)[2:7] <- paste("geod.c.", colnames(geod.c)[2:7], sep="")

# combine all in single table
all <- merge(base.f, bias.f, all=T)
all <- merge(all, base.c, all=T)
all <- merge(all, geod.c, all=T)

write.csv(all, file="C:/data/modelling/output/modelcomparison.csv")
all <- read.csv("C:/data/modelling/output/modelcomparison.csv")
all <- all[complete.cases(all),]

bla <- rbind(data.frame("model"=as.factor("BIOCLIM"), "run"=as.factor("baseline"), "metric"=as.factor("auc"), "value"=all[,"base.f.BIO.auc"]),
             data.frame("model"=as.factor("MAXENT"), "run"=as.factor("baseline"), "metric"=as.factor("auc"), "value"=all[,"base.f.MAX.auc"]),
             data.frame("model"=as.factor("BRT"), "run"=as.factor("baseline"), "metric"=as.factor("auc"), "value"=all[,"base.f.BRT.auc"]),
             data.frame("model"=as.factor("BIOCLIM"), "run"=as.factor("baseline"), "metric"=as.factor("cor"), "value"=all[,"base.f.BIO.cor"]),
             data.frame("model"=as.factor("MAXENT"), "run"=as.factor("baseline"), "metric"=as.factor("cor"), "value"=all[,"base.f.MAX.cor"]),
             data.frame("model"=as.factor("BRT"), "run"=as.factor("baseline"), "metric"=as.factor("cor"), "value"=all[,"base.f.BRT.cor"]),
             
             data.frame("model"=as.factor("BIOCLIM"), "run"=as.factor("bias"), "metric"=as.factor("auc"), "value"=all[,"bias.f.BIO.auc"]),
             data.frame("model"=as.factor("MAXENT"), "run"=as.factor("bias"), "metric"=as.factor("auc"), "value"=all[,"bias.f.MAX.auc"]),
             data.frame("model"=as.factor("BRT"), "run"=as.factor("bias"), "metric"=as.factor("auc"), "value"=all[,"bias.f.BRT.auc"]),
             data.frame("model"=as.factor("BIOCLIM"), "run"=as.factor("bias"), "metric"=as.factor("cor"), "value"=all[,"bias.f.BIO.cor"]),
             data.frame("model"=as.factor("MAXENT"), "run"=as.factor("bias"), "metric"=as.factor("cor"), "value"=all[,"bias.f.MAX.cor"]),
             data.frame("model"=as.factor("BRT"), "run"=as.factor("bias"), "metric"=as.factor("cor"), "value"=all[,"bias.f.BRT.cor"]),
             
             data.frame("model"=as.factor("BIOCLIM"), "run"=as.factor("clipped"), "metric"=as.factor("auc"), "value"=all[,"base.c.BIO.auc"]),
             data.frame("model"=as.factor("MAXENT"), "run"=as.factor("clipped"), "metric"=as.factor("auc"), "value"=all[,"base.c.MAX.auc"]),
             data.frame("model"=as.factor("BRT"), "run"=as.factor("clipped"), "metric"=as.factor("auc"), "value"=all[,"base.c.BRT.auc"]),
             data.frame("model"=as.factor("BIOCLIM"), "run"=as.factor("clipped"), "metric"=as.factor("cor"), "value"=all[,"base.c.BIO.cor"]),
             data.frame("model"=as.factor("MAXENT"), "run"=as.factor("clipped"), "metric"=as.factor("cor"), "value"=all[,"base.c.MAX.cor"]),
             data.frame("model"=as.factor("BRT"), "run"=as.factor("clipped"), "metric"=as.factor("cor"), "value"=all[,"base.c.BRT.cor"]),
             
             data.frame("model"=as.factor("BIOCLIM"), "run"=as.factor("expert"), "metric"=as.factor("auc"), "value"=all[,"geod.c.BIO.auc"]),
             data.frame("model"=as.factor("MAXENT"), "run"=as.factor("expert"), "metric"=as.factor("auc"), "value"=all[,"geod.c.MAX.auc"]),
             data.frame("model"=as.factor("BRT"), "run"=as.factor("expert"), "metric"=as.factor("auc"), "value"=all[,"geod.c.BRT.auc"]),
             data.frame("model"=as.factor("BIOCLIM"), "run"=as.factor("expert"), "metric"=as.factor("cor"), "value"=all[,"geod.c.BIO.cor"]),
             data.frame("model"=as.factor("MAXENT"), "run"=as.factor("expert"), "metric"=as.factor("cor"), "value"=all[,"geod.c.MAX.cor"]),
             data.frame("model"=as.factor("BRT"), "run"=as.factor("expert"), "metric"=as.factor("cor"), "value"=all[,"geod.c.BRT.cor"]))

##
##

bla <- rbind(data.frame("model"=as.factor("BIOCLIM"), "run"=as.factor("baseline"), "metric"=as.factor("auc"), "value"=all[,"all.auc.con.BIO"]),
             data.frame("model"=as.factor("MAXENT"), "run"=as.factor("baseline"), "metric"=as.factor("auc"), "value"=all[,"all.auc.con.MAX"]),
             data.frame("model"=as.factor("BRT"), "run"=as.factor("baseline"), "metric"=as.factor("auc"), "value"=all[,"all.auc.con.BRT"]),
             data.frame("model"=as.factor("BIOCLIM"), "run"=as.factor("baseline"), "metric"=as.factor("cor"), "value"=all[,"all.cor.con.BIO"]),
             data.frame("model"=as.factor("MAXENT"), "run"=as.factor("baseline"), "metric"=as.factor("cor"), "value"=all[,"all.cor.con.MAX"]),
             data.frame("model"=as.factor("BRT"), "run"=as.factor("baseline"), "metric"=as.factor("cor"), "value"=all[,"all.cor.con.BRT"]),
             
             data.frame("model"=as.factor("BIOCLIM"), "run"=as.factor("bias"), "metric"=as.factor("auc"), "value"=all[,"all.auc.con.BIO  bias.f.BIO.auc"]),
             data.frame("model"=as.factor("MAXENT"), "run"=as.factor("bias"), "metric"=as.factor("auc"), "value"=all[,"bias.f.MAX.auc"]),
             data.frame("model"=as.factor("BRT"), "run"=as.factor("bias"), "metric"=as.factor("auc"), "value"=all[,"bias.f.BRT.auc"]),
             data.frame("model"=as.factor("BIOCLIM"), "run"=as.factor("bias"), "metric"=as.factor("cor"), "value"=all[,"bias.f.BIO.cor"]),
             data.frame("model"=as.factor("MAXENT"), "run"=as.factor("bias"), "metric"=as.factor("cor"), "value"=all[,"bias.f.MAX.cor"]),
             data.frame("model"=as.factor("BRT"), "run"=as.factor("bias"), "metric"=as.factor("cor"), "value"=all[,"bias.f.BRT.cor"]),
             
             data.frame("model"=as.factor("BIOCLIM"), "run"=as.factor("clipped"), "metric"=as.factor("auc"), "value"=all[,"base.c.BIO.auc"]),
             data.frame("model"=as.factor("MAXENT"), "run"=as.factor("clipped"), "metric"=as.factor("auc"), "value"=all[,"base.c.MAX.auc"]),
             data.frame("model"=as.factor("BRT"), "run"=as.factor("clipped"), "metric"=as.factor("auc"), "value"=all[,"base.c.BRT.auc"]),
             data.frame("model"=as.factor("BIOCLIM"), "run"=as.factor("clipped"), "metric"=as.factor("cor"), "value"=all[,"base.c.BIO.cor"]),
             data.frame("model"=as.factor("MAXENT"), "run"=as.factor("clipped"), "metric"=as.factor("cor"), "value"=all[,"base.c.MAX.cor"]),
             data.frame("model"=as.factor("BRT"), "run"=as.factor("clipped"), "metric"=as.factor("cor"), "value"=all[,"base.c.BRT.cor"]),
             
             data.frame("model"=as.factor("BIOCLIM"), "run"=as.factor("expert"), "metric"=as.factor("auc"), "value"=all[,"geod.c.BIO.auc"]),
             data.frame("model"=as.factor("MAXENT"), "run"=as.factor("expert"), "metric"=as.factor("auc"), "value"=all[,"geod.c.MAX.auc"]),
             data.frame("model"=as.factor("BRT"), "run"=as.factor("expert"), "metric"=as.factor("auc"), "value"=all[,"geod.c.BRT.auc"]),
             data.frame("model"=as.factor("BIOCLIM"), "run"=as.factor("expert"), "metric"=as.factor("cor"), "value"=all[,"geod.c.BIO.cor"]),
             data.frame("model"=as.factor("MAXENT"), "run"=as.factor("expert"), "metric"=as.factor("cor"), "value"=all[,"geod.c.MAX.cor"]),
             data.frame("model"=as.factor("BRT"), "run"=as.factor("expert"), "metric"=as.factor("cor"), "value"=all[,"geod.c.BRT.cor"]))


# PLOT ----
# select what to plot
#plot <- bla[bla$run %in% c("baseline", "bias") & bla$metric=="cor",]

library(ggplot2)
#library(devtools)
library(easyGgplot2)

tiff("C:/Users/LaurensG/Desktop/plots/barplots_sb1_SAUP.tiff", width=16, height=12, unit="cm", res=600)
ggplot2.boxplot(data=data.plot, xName="run", yName="value", groupName="model",
                faceting=TRUE, facetingVarNames="model", facetingDirection="horizontal", facetingFont=c(12, "bold", "black"),
                width=1, notch=T, outlier.shape=21, groupColors=c('#999999','#E69F00','#56B4E9'), showLegend=FALSE,
                #axisLine=c(0.5, "solid", "black"),  backgroundColor="white", removePanelGrid=TRUE,removePanelBorder=TRUE,
                xtitle="Modelrun", ytitle="AUC") + 
                theme(axis.text.x=element_text(angle=-45, hjust=0, color="gray40"), axis.text.y=element_text(color="gray40"), 
                      axis.title.x=element_text(vjust=-0.5, hjust=0.5), axis.title.y=element_text(vjust=1.5),
                      plot.margin=unit(c(2,10,2,2),"mm")) + expand_limits(y=c(0.6,1))
dev.off()

#
#
#

ggplot2.boxplot(data=data.plot, xName="run", yName="cor", groupName="model",
                faceting=TRUE, facetingVarNames="model", facetingDirection="horizontal", facetingFont=c(12, "bold", "black"),
                width=1, notch=T, outlier.shape=21, groupColors=c('#999999','#E69F00','#56B4E9'), showLegend=FALSE,
                #axisLine=c(0.5, "solid", "black"),  backgroundColor="white", removePanelGrid=TRUE,removePanelBorder=TRUE,
                xtitle="Modelrun", ytitle="COR") + 
  theme(axis.text.x=element_text(angle=-45, hjust=0, color="gray40"), axis.text.y=element_text(color="gray40"), axis.title.x=element_text(vjust=-0.5, hjust=0.5), axis.title.y=element_text(vjust=1.5)) + expand_limits(y=c(0.6,1))




#*#*#*#*#*#*#*#*
# Add model results to table ----

speciestable <- read.csv("C:/Users/LaurensG/Desktop/speciestable.csv")

resIUCN <- read.csv("C:/Users/LaurensG/Desktop/completecomparison_SAUP.csv")
resSAUP <- read.csv("C:/Users/LaurensG/Desktop/completecomparison_IUCN.csv")
resBOTH <- read.csv("C:/Users/LaurensG/Desktop/completecomparison_IUCN_SAUP.csv")

stab <- merge(speciestable[,c(1:13)], resIUCN, by="TaxonKey", all.x=T)
stab <- merge(stab, resSAUP, by="TaxonKey", all.x=T)
stab <- merge(stab, resBOTH, by="TaxonKey", all.x=T)

stab[,c("occcount", "rsIUCN", "rsSAUP", "rsBOTH", "jaccard", "dempel", "habitat",)]
pairs(stab[,c("occcount", "rsIUCN", "rsSAUP", "rsBOTH", "jaccard", "dempel", "habitat")])

bla <- cbind(stab[,c("occcount", "rsIUCN", "rsSAUP", "rsBOTH", "jaccard", "dempel", "habitat")], 
      
      #normal vs bias, conventional evaluation
      "B.a.c.d.BIO"=stab$IUCN_SAUP.auc.con.BIO.b - stab$IUCN_SAUP.auc.con.BIO,
      "B.a.c.d.MAX"=stab$IUCN_SAUP.auc.con.MAX.b - stab$IUCN_SAUP.auc.con.MAX,
      "B.a.c.d.BRT"=stab$IUCN_SAUP.auc.con.BRT.b - stab$IUCN_SAUP.auc.con.BRT,
      
      "B.c.c.d.BIO"=stab$IUCN_SAUP.cor.con.BIO.b - stab$IUCN_SAUP.cor.con.BIO,
      "B.c.c.d.MAX"=stab$IUCN_SAUP.cor.con.MAX.b - stab$IUCN_SAUP.cor.con.MAX,
      "B.c.c.d.BRT"=stab$IUCN_SAUP.cor.con.BRT.b - stab$IUCN_SAUP.cor.con.BRT,
      
      #normal vs bias, expert evaluation
      "B.a.e.d.BIO"=stab$IUCN_SAUP.auc.exp.BIO.b - stab$IUCN_SAUP.auc.exp.BIO,
      "B.a.e.d.MAX"=stab$IUCN_SAUP.auc.exp.MAX.b - stab$IUCN_SAUP.auc.exp.MAX,
      "B.a.e.d.BRT"=stab$IUCN_SAUP.auc.exp.BRT.b - stab$IUCN_SAUP.auc.exp.BRT,
      
      "B.c.e.d.BIO"=stab$IUCN_SAUP.cor.exp.BIO.b - stab$IUCN_SAUP.cor.exp.BIO,
      "B.c.e.d.MAX"=stab$IUCN_SAUP.cor.exp.MAX.b - stab$IUCN_SAUP.cor.exp.MAX,
      "B.c.e.d.BRT"=stab$IUCN_SAUP.cor.exp.BRT.b - stab$IUCN_SAUP.cor.exp.BRT,

      #bias models, difference between conventional and expert evaluation
      "B.a.m.d.BIO"=stab$IUCN_SAUP.auc.exp.BIO.b - stab$IUCN_SAUP.auc.con.BIO.b,
      "B.a.m.d.MAX"=stab$IUCN_SAUP.auc.exp.MAX.b - stab$IUCN_SAUP.auc.con.MAX.b,
      "B.a.m.d.BRT"=stab$IUCN_SAUP.auc.exp.BRT.b - stab$IUCN_SAUP.auc.con.BRT.b,
      
      "B.c.m.d.BIO"=stab$IUCN_SAUP.cor.exp.BIO.b - stab$IUCN_SAUP.cor.con.BIO.b,
      "B.c.m.d.MAX"=stab$IUCN_SAUP.cor.exp.MAX.b - stab$IUCN_SAUP.cor.con.MAX.b,
      "B.c.m.d.BRT"=stab$IUCN_SAUP.cor.exp.BRT.b - stab$IUCN_SAUP.cor.con.BRT.b)


boxplot(B.a.c.d.MAX~dempel, data=bla)
boxplot(B.a.e.d.MAX~dempel, data=bla)
boxplot(B.a.m.d.MAX~dempel, data=bla)

boxplot(B.a.c.d.BRT~dempel, data=bla)
boxplot(B.a.e.d.BRT~dempel, data=bla)
boxplot(B.a.m.d.BRT~dempel, data=bla)


boxplot(B.c.c.d.MAX~dempel, data=bla)
boxplot(B.c.e.d.MAX~dempel, data=bla)
boxplot(B.c.m.d.MAX~dempel, data=bla)

boxplot(B.c.c.d.BRT~dempel, data=bla)
boxplot(B.c.e.d.BRT~dempel, data=bla)
boxplot(B.c.m.d.BRT~dempel, data=bla)

bla2 <- bla[,c(8:25)]
bla2 <- bla2[complete.cases(bla2),]
bla3 <- kmeans(bla2, centers=3)

library(cluster)
library(HSAUR)
data(pottery)
km    <- kmeans(pottery,3)
dissE <- daisy(pottery) 
dE2   <- dissE^2
sk2   <- silhouette(km$cl, dE2)
plot(sk2)

library(cluster)
library(fpc)

data(iris)
dat <- iris[, -5] # without known classification 
# Kmeans clustre analysis
clus <- kmeans(dat, centers=3)
# Fig 01
plotcluster(dat, clus$cluster)
##
##

with(bla2, pairs(bla2, col=c(1:3)[bla3$cluster])) 
