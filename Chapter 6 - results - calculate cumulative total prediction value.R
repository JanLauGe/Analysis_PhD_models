setwd("C:/data/modelling")
data.extent <- raster("C:/data/modelling/data/data.extent.grd")
load("data/species/2_data.species.aggregated.RData")
load("data/data.species.RData")
a <- data.species[data.species$set=="background",c("x", "y")]

species.list1 <- unlist(strsplit(unlist(strsplit(list.files("C:/data/modelling/output/geodist2", pattern=".RData$"), split="SDMeval_geodist_")), split=".RData"))
species.list2 <- unlist(strsplit(unlist(strsplit(list.files("C:/data/modelling/output/baseline", pattern="*SDMeval.clipped*"), split="SDMeval.clipped_")), split=".RData"))

#normalise data
#x is the dataframe to be normalised, y is the reference scale data (can be x to rescale 0-1)
normalise <- function(x){
  x <- x[!is.na(x)]
  norm <- (x - min(x))/(max(x) - min(x))
  return(norm)}


all.preds <- list("baseline.BIO"=data.extent,
                  "baseline.MAX"=data.extent,
                  "baseline.BRT"=data.extent,
                  "bias.BIO"=data.extent,
                  "bias.MAX"=data.extent,
                  "bias.BRT"=data.extent,
                  "clipped.BIO"=data.extent,
                  "clipped.MAX"=data.extent,
                  "clipped.BRT"=data.extent,
                  "geodist.BIO"=data.extent,
                  "geodist.MAX"=data.extent,
                  "geodist.BRT"=data.extent)

s=1
species.list=species.list1[species.list1 %in% species.list2]
app=1
approach=c("baseline", "bias", "clipped", "geodistance")
paths=c("baseline/SDMeval_fulld_", "bias1/SDMeval_fulld_sb1_", "baseline/SDMeval.clipped_", "geodist2/SDMeval_geodist_")
m=1
models=c("BIOCLIM","MAXENT","BRT")
all.crhs=NULL
all.thres=NULL
all.cnpc=NULL
for(s in 1:length(species.list)){
  species=species.list[s]
  print(species)

  p <- data.species[data.species$TaxonKey==species.list[s] & data.species$set=="test",c("x","y")]
  pa <- data.species.aggregated[data.species.aggregated$TaxonKey==species.list[s] & data.species.aggregated$set=="test",c("x","y")]
  
  pred.alla=NULL
  for(app in 1:length(approach)){
    path=paths[app]
    load(paste("C:/data/modelling/output/", path, species.list[s], ".RData", sep=""))
    for(m in 1:length(models)){
      #create raster
      pred <- results[[2]][[m]]
      values(pred)[!is.na(values(pred))] <- normalise(values(pred))
      
      #thresholds
      eval <- evaluate(extract(pred, SpatialPoints(p)), extract(pred, SpatialPoints(a)))
      evalt <- threshold(eval)
      thres <- data.frame("TaxonKey"=species.list[s], "approach"=approach[app], "model"=models[m], "kappa"=evalt[[1]], "spec_sens"=evalt[[2]], "no_omission"=evalt[[3]], "prevalence"=evalt[[4]], "equal_sens_spec"=evalt[[5]], "sensitivity"=evalt[[6]])
      all.thres <- rbind(all.thres, thres)
      
      #get cumulative normalised RHS value
      crhs <- data.frame("TaxonKey"=species.list[s], "approach"=approach[app], "model"=models[m], "sum"=sum(values(pred)[!is.na(values(pred))]))
      all.crhs <- rbind(all.crhs, crhs)
      
      #get number of total predicted cells per threshold
      cnpc <- cbind(thres[,1:3], "kappa"=length(pred[pred >= thres[,4]]), "spec_sens"=length(pred[pred >= thres[,5]]), "no_omission"=length(pred[pred >= thres[,6]]), "prevalence"=length(pred[pred >= thres[,7]]), "equal_sens_spec"=length(pred[pred >= thres[,8]]), "sensitivity"=length(pred[pred >= thres[,9]]))
      all.cnpc <- rbind(all.cnpc, cnpc)
      
      #add raster to total cumulative rhs prediction
      all.preds[[m+3*(app-1)]] <- all.preds[[m+3*(app-1)]] + pred
    }
  }
}

save(all.preds, file="C:/Users/LaurensG/Desktop/all predictions raster list.RData")
write.csv(all.crhs, file="C:/Users/LaurensG/Desktop/all predictions cumulative rhs.csv")
write.csv(all.thres, file="C:/Users/LaurensG/Desktop/all predictions thresholds.csv")
write.csv(all.cnpc, file="C:/Users/LaurensG/Desktop/all predictions number of cells per threshold.csv")

plot(all.crhs[all.crhs$approach=="clipped","sum"], all.crhs[all.crhs$approach=="geodistance","sum"], col=all.crhs$model)
plot(log(all.crhs[all.crhs$approach=="baseline","sum"]), log(all.crhs[all.crhs$approach=="clipped","sum"]), col=all.crhs$model)

#threshold kappa clipped vs geodist
plot(all.thres[all.thres$approach=="clipped","kappa"], all.thres[all.thres$approach=="geodistance","kappa"], col=all.thres$model)
#bioclim the same, or lower, brt low and varying, maxent high and varying

#threshold kappa baseline vs bias
plot(all.thres[all.thres$approach=="baseline","kappa"], all.thres[all.thres$approach=="bias","kappa"], col=all.thres$model)
#bioclim the same, maxent high, brt low

pdata <- data.frame(cbind("baseline"=as.numeric(log(all.crhs[all.crhs$approach=="baseline","sum"])), "clipped"=as.numeric(log(all.crhs[all.crhs$approach=="clipped","sum"])), "model"=all.crhs[,"model"]))
qplot(x=baseline, y=clipped, data=pdata, colour=model)

#Number of predicted cells for each approach and model
sum(all.cnpc[all.cnpc$approach=="baseline" & all.cnpc$model=="BIOCLIM", "kappa"])
sum(all.cnpc[all.cnpc$approach=="clipped" & all.cnpc$model=="BIOCLIM", "kappa"])
sum(all.cnpc[all.cnpc$approach=="geodistance" & all.cnpc$model=="BIOCLIM", "kappa"])

sum(all.cnpc[all.cnpc$approach=="baseline" & all.cnpc$model=="MAXENT", "kappa"])
sum(all.cnpc[all.cnpc$approach=="clipped" & all.cnpc$model=="MAXENT", "kappa"])
sum(all.cnpc[all.cnpc$approach=="geodistance" & all.cnpc$model=="MAXENT", "kappa"])

sum(all.cnpc[all.cnpc$approach=="baseline" & all.cnpc$model=="BRT", "kappa"])
sum(all.cnpc[all.cnpc$approach=="clipped" & all.cnpc$model=="BRT", "kappa"])
sum(all.cnpc[all.cnpc$approach=="geodistance" & all.cnpc$model=="BRT", "kappa"])

###
###
ggplot(data=pdata, xName="run", yName="value", groupName="model") +
geom_point(aes(baseline, clipped, col=factor(model)))
+ geom_point(alpha = 1/10)
           
           
           , faceting=TRUE, facetingVarNames="model", facetingDirection="horizontal", facetingFont=c(12, "bold", "black"),
width=1, notch=T, outlier.shape=21, groupColors=c('#999999','#E69F00','#56B4E9'), showLegend=FALSE)
,
#axisLine=c(0.5, "solid", "black"),  backgroundColor="white", removePanelGrid=TRUE,removePanelBorder=TRUE,
xtitle="Modelrun", ytitle="AUC") + #coord_cartesian(ylim = c(0.9,1)) +
  #scale_y_continuous(limits = c(0, 1)) +
  theme(axis.text.x=element_text(angle=-45, hjust=0, color="gray40"), axis.text.y=element_text(color="gray40"), 
        axis.title.x=element_text(vjust=-0.5, hjust=0.5), axis.title.y=element_text(vjust=1.5),
        plot.margin=unit(c(2,10,2,2),"mm")) + expand_limits(y=c(0.6,1))




bla <- ggplot(data=pdata, aes(baseline, clipped)) +
  geom_point(colour=qdata$model, size=1)
geom_point(aes(baseline, clipped, col=model))
