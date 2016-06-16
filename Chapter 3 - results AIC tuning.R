library(sp)
library(raster)
library(dismo)
library(rgdal)

Files <- list.files("G:/HPC/output/aic/", pattern="AIC_summary_")

ModelEvals=NULL
File=Files[1]
for(File in Files){
  ModelEval <- read.csv(paste("G:/HPC/output/aic/", File, sep=""))
  ModelEvals <- rbind(ModelEvals, ModelEval)
}

#Plot some graphs ----
boxplot(ModelEvals[,c("BIO.auc", "MAX.auc", "BRT.auc")])
boxplot(ModelEvals[,c("BIO.cor", "MAX.cor", "BRT.cor")])
colnames(ModelEvals) <- c("X", "TaxonKey", "beta", "npar", "MAX.AICc", "bag", "BRT.AICc", "logLik", "AUC.BIO.aic", "AUC.MAX.aic", "AUC.BRT.aic", "COR.BIO.aic", "COR.MAX.aic", "COR.BRT.aic")
ModelEvals <- ModelEvals[,2:length(ModelEvals)]

EvalTable <- merge(ModelEvalsBaseline, ModelEvals, by="TaxonKey")
attach(EvalTable)

#Change in AUC is slightly positive for MAXENT, inconclusive for BRT
summary(AUC.MAX - AUC.MAX.aic)
summary(AUC.BRT - AUC.BRT.aic)
t.test(AUC.MAX, AUC.MAX.aic, paired=T) #highly significant
t.test(AUC.BRT, AUC.BRT.aic, paired=T) #not significant

plot(log(NPresences), jitter(beta))
cor.test(log(NPresences), beta)
plot(log(NPresences), jitter(bag))
cor.test(log(NPresences), bag)

plot(EvalTable$AUC.BIO ~ EvalTable$AUC.BIO.aic)
plot(EvalTable$AUC.MAX ~ EvalTable$AUC.MAX.aic)
plot(EvalTable$AUC.BRT ~ EvalTable$AUC.BRT.aic)

plot(EvalTable$NPresences ~ jitter(EvalTable$beta))
plot(EvalTable$NPresences ~ jitter(EvalTable$bag))

cor(EvalTable[,c("NPresences", "beta")])
plot(EvalTable[,c("NPresences", "beta")])
cor.test(log(EvalTable$NPresences), EvalTable$beta)
plot(cbind(logN=log(EvalTable$NPresences),beta=jitter(EvalTable$beta, factor=3)))

t.test(EvalTable[,"NPresences"], EvalTable[,"beta"])
cor(EvalTable[,c("NPresences", "bag")])

plot(jitter(EvalTable$beta) ~ EvalTable$AUC.MAX)
plot(jitter(EvalTable$beta) ~ EvalTable$COR.MAX)
plot(jitter(EvalTable$bag) ~ EvalTable$AUC.BRT)
plot(jitter(EvalTable$bag) ~ EvalTable$COR.BRT)

plotdata <- rbind(data.frame("TaxonKey"=all[,"TaxonKey"], "model"=as.factor("BIOCLIM"), "run"=as.factor("baseline"), "metric"=as.factor("auc"), "value"=all[,"BIO.auc"]),
                  data.frame("TaxonKey"=all[,"TaxonKey"], "model"=as.factor("MAXENT"), "run"=as.factor("baseline"), "metric"=as.factor("auc"), "value"=all[,"MAX.auc"]),
                  data.frame("TaxonKey"=all[,"TaxonKey"], "model"=as.factor("BRT"), "run"=as.factor("baseline"), "metric"=as.factor("auc"), "value"=all[,"BRT.auc"]),
                  data.frame("TaxonKey"=all[,"TaxonKey"], "model"=as.factor("BIOCLIM"), "run"=as.factor("baseline"), "metric"=as.factor("cor"), "value"=all[,"BIO.cor"]),
                  data.frame("TaxonKey"=all[,"TaxonKey"], "model"=as.factor("MAXENT"), "run"=as.factor("baseline"), "metric"=as.factor("cor"), "value"=all[,"MAX.cor"]),
                  data.frame("TaxonKey"=all[,"TaxonKey"], "model"=as.factor("BRT"), "run"=as.factor("baseline"), "metric"=as.factor("cor"), "value"=all[,"BRT.cor"]))

library(ggplot2)
library(easyGgplot2)
test.geodist <- plotdata
# ~ Plot1 ----
# select what to plot
data.plot <- test.geodist[test.geodist$model %in% c("BIOCLIM", "MAXENT", "BRT") & test.geodist$run %in% c("baseline", "range map", "clipped") & test.geodist$metric=="auc",]
tiff("C:/Users/LaurensG/Desktop/plots/barplots_baseline_auc.tiff", width=16, height=12, unit="cm", res=600)
ggplot2.boxplot(data=data.plot, yName="value", groupName="model",
                faceting=TRUE, facetingVarNames="model", facetingDirection="horizontal", facetingFont=c(12, "bold", "black"),
                width=1, notch=T, outlier.shape=21, groupColors=c('#999999','#E69F00','#56B4E9'), showLegend=FALSE,
                #axisLine=c(0.5, "solid", "black"),  backgroundColor="white", removePanelGrid=TRUE,removePanelBorder=TRUE,
                xtitle="Model", ytitle="AUC") + #coord_cartesian(ylim = c(0.9,1)) +
  #scale_y_continuous(limits = c(0, 1)) +
  theme(axis.text.x=element_text(angle=-45, hjust=0, color="gray40"), axis.text.y=element_text(color="gray40"), 
        axis.title.x=element_text(vjust=-0.5, hjust=0.5), axis.title.y=element_text(vjust=1.5),
        plot.margin=unit(c(2,10,2,2),"mm")) + expand_limits(y=c(0.6,1))
dev.off()

# ~ Plot2 ----
# select what to plot
data.plot <- test.geodist[test.geodist$model %in% c("BIOCLIM", "MAXENT", "BRT") & test.geodist$run %in% c("baseline", "range map", "clipped") & test.geodist$metric=="cor",]
tiff("C:/Users/LaurensG/Desktop/plots/barplots_baseline_cor.tiff", width=16, height=12, unit="cm", res=600)
ggplot2.boxplot(data=data.plot, yName="value", groupName="model",
                faceting=TRUE, facetingVarNames="model", facetingDirection="horizontal", facetingFont=c(12, "bold", "black"),
                width=1, notch=T, outlier.shape=21, groupColors=c('#999999','#E69F00','#56B4E9'), showLegend=FALSE,
                #axisLine=c(0.5, "solid", "black"),  backgroundColor="white", removePanelGrid=TRUE,removePanelBorder=TRUE,
                xtitle="Model", ytitle="COR") + #coord_cartesian(ylim = c(0.9,1)) +
  #scale_y_continuous(limits = c(0, 1)) +
  theme(axis.text.x=element_text(angle=-45, hjust=0, color="gray40"), axis.text.y=element_text(color="gray40"), 
        axis.title.x=element_text(vjust=-0.5, hjust=0.5), axis.title.y=element_text(vjust=1.5),
        plot.margin=unit(c(2,10,2,2),"mm")) + expand_limits(y=c(0.6,1))
dev.off()
