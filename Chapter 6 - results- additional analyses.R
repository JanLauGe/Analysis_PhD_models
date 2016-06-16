library(raster)
library(rasterVis)
library(dismo)
library(plyr)
library(rgdal)
library(reshape)

data.extent <- raster("C:/data/modelling/data/data.extent.grd")
data.land <- readOGR("C:/data/modelling/data/envi/othershapefiles", layer="land")
data.specinfo <- read.csv("C:/data/modelling/data/species.list.csv")
data.files <- list.files("G:/HPC/output/finalmodel", pattern="SDMeval_clipped_finalmodel_", full.names=T) 

#Function to read all result maps for each of the species
GetResults <- function(file){
  load(file)
  result.maps.bio <- stack(lapply(results, FUN=function(x) x[[2]][[1]]))  
  result.maps.max <- stack(lapply(results, FUN=function(x) x[[2]][[2]]))
  
  if(class(results[[1]][[2]][[3]]) != "RasterLayer"){
    result.maps.brt <- NULL
  }else{
    result.maps.brt <- stack(lapply(results, FUN=function(x) x[[2]][[3]]))
  }  
  result.maps <- list(result.maps.bio, result.maps.max, result.maps.brt)
  return(result.maps)
}

SaveResults <- function(stack=RasterStack(NULL), specinfo){
  spec.pred.list = list(NULL)
  spec.vari.list = list(NULL)
  for(model in 1:3){
    if (model==1){
      modelname="BIOCLIM"
    }else if (model==2){
      modelname="MAXENT"
    }else if (model==3){
      modelname="BRT"
    }else{
      print("Error! Invalid model identifier")
    }
    if(class(stack[[model]]) == "NULL"){
      write.table("NULL", file=paste("C:/data/modelling/output/chapter6/grids/", specinfo[,"TaxonKey"], "_", model, "_mean", ".txt", sep=""))
      write.table("NULL", file=paste("C:/data/modelling/output/chapter6/grids/", specinfo[,"TaxonKey"], "_", model, "_vari", ".txt", sep=""))
      
      jpeg(paste("C:/data/modelling/output/chapter6/mapplots/", "Pred_", specinfo[,"TaxonKey"], "_", modelname, ".jpeg", sep=""), height=10, width=16, unit="cm", res=300)
      dev.off()
      jpeg(paste("C:/data/modelling/output/chapter6/mapplots/", "Vari_", specinfo[,"TaxonKey"], "_", modelname, ".jpeg", sep=""), height=10, width=16, unit="cm", res=300)
      dev.off()
      
    }else{
      spec.pred <- mean(stack[[model]])
      x = values(spec.pred); values(spec.pred) = (x-min(na.omit(x)))/(max(na.omit(x))-min(na.omit(x)))
      spec.vari <- calc(stack[[model]], fun=var)
    
      writeRaster(spec.pred, file=paste("C:/data/modelling/output/chapter6/grids/", specinfo[,"TaxonKey"], "_", model, "_mean", ".grd", sep=""), overwrite=T)
      writeRaster(spec.vari, file=paste("C:/data/modelling/output/chapter6/grids/", specinfo[,"TaxonKey"], "_", model, "_vari", ".grd", sep=""), overwrite=T)
      
      jpeg(paste("C:/data/modelling/output/chapter6/mapplots/", "Pred_", specinfo[,"TaxonKey"], "_", modelname, ".jpeg", sep=""), height=10, width=16, unit="cm", res=300)
        print(levelplot(spec.pred, margin=F, contour=F, colorkey=list(space="right"), pretty=T, par.settings=BuRdTheme, main=paste(specinfo[,"TaxonName"], " - ", modelname, "models")) + layer(sp.polygons(data.land, fill="black")))
      dev.off()
      jpeg(paste("C:/data/modelling/output/chapter6/mapplots/", "Vari_", specinfo[,"TaxonKey"], "_", modelname, ".jpeg", sep=""), height=10, width=16, unit="cm", res=300)
        print(levelplot(spec.vari, margin=F, contour=F, colorkey=list(space="right"), pretty=T, par.settings=BuRdTheme, main=paste(specinfo[,"TaxonName"], " - ", modelname, "models - variance")) + layer(sp.polygons(data.land, fill="black")))
      dev.off()
      
      spec.pred.list = append(spec.pred.list, spec.pred)
      spec.vari.list = append(spec.vari.list, spec.vari)
    }
  }
  return(list(spec.pred, spec.vari))
}

results.allpred.bioclim = data.extent
results.allpred.maxent = data.extent
results.allpred.brt = data.extent
results.allvari.bioclim = data.extent
results.allvari.maxent = data.extent
results.allvari.brt = data.extent
s = 1
#s=127
#s=462
for (s in 1:length(data.files)){
  TaxonKey = unlist(strsplit(unlist(strsplit(data.files[s], split="G:/HPC/output/finalmodel/SDMeval_clipped_finalmodel_")), split=".RData"))
  SpecInfo = data.specinfo[data.specinfo$TaxonKey==TaxonKey,]
  
  result.maps.individual <- GetResults(data.files[s])
  result.maps.aggregated <- SaveResults(stack=result.maps.individual, specinfo=SpecInfo)

  #Add up layers for mean overall prediction
  results.allpred.bioclim = result.maps.aggregated[[1]] + results.allpred.bioclim
  results.allpred.maxent = result.maps.aggregated[[1]] + results.allpred.maxent
  results.allpred.brt = result.maps.aggregated[[1]] + results.allpred.brt
  
  #Add up layers for mean overall variance
  results.allvari.bioclim = result.maps.aggregated[[1]] + results.allvari.bioclim
  results.allvari.maxent = result.maps.aggregated[[1]] + results.allvari.maxent
  results.allvari.brt = result.maps.aggregated[[1]] + results.allvari.brt
  
  print(paste("... species", SpecInfo[,"TaxonKey"], "done!"))
}

#Create joined predictions from grid ----
all.pred.bio <- data.extent
for (f in list.files("C:/data/modelling/output/chapter6/grids", pattern="1_mean.grd", full.names=T)){
  this.raster <- raster(f)
  this.raster[is.na(this.raster)] <- 0
  all.pred.bio <- all.pred.bio + this.raster
  print(paste("...done with", f))
}
writeRaster(all.pred.bio, file="C:/data/modelling/output/chapter6/aggregatedprediction_bio.grd")
all.pred.max <- data.extent
for (f in list.files("C:/data/modelling/output/chapter6/grids", pattern="2_mean.grd", full.names=T)){
  this.raster <- raster(f)
  this.raster[is.na(this.raster)] <- 0
  all.pred.max <- all.pred.max + this.raster
  print(paste("...done with", f))
}
writeRaster(all.pred.max, file="C:/data/modelling/output/chapter6/aggregatedprediction_max.grd", overwrite=T)
all.pred.brt <- data.extent
for (f in list.files("C:/data/modelling/output/chapter6/grids", pattern="3_mean.grd", full.names=T)){
  this.raster <- raster(f)
  this.raster[is.na(this.raster)] <- 0
  all.pred.brt <- all.pred.brt + this.raster
  print(paste("...done with", f))
}
writeRaster(all.pred.brt, file="C:/data/modelling/output/chapter6/aggregatedprediction_brt.grd", overwrite=T)

#Get the same maps for the clipped baseline
baseline.pred.bio <- data.extent
baseline.pred.max <- data.extent
baseline.pred.brt <- data.extent
for (f in list.files("G:/HPC/output/baseline", pattern="SDMeval_clipped_baseline_", full.names=T)){
  load(f)
  
  this.raster <- results[[2]][[1]]
  this.raster[is.na(this.raster)] <- 0
  baseline.pred.bio <- baseline.pred.bio + this.raster
  this.raster <- results[[2]][[2]]
  this.raster[is.na(this.raster)] <- 0
  baseline.pred.max <- baseline.pred.max + this.raster
  this.raster <- results[[2]][[3]]
  this.raster[is.na(this.raster)] <- 0
  baseline.pred.brt <- baseline.pred.brt + this.raster
  print(paste("...done with", f))
}
writeRaster(baseline.pred.bio, file="C:/data/modelling/output/chapter6/baseline_bio.asc", overwrite=T)
writeRaster(baseline.pred.max, file="C:/data/modelling/output/chapter6/baseline_max.asc", overwrite=T)
writeRaster(baseline.pred.brt, file="C:/data/modelling/output/chapter6/baseline_brt.asc", overwrite=T)

writeRaster(all.pred.bio, file="C:/data/modelling/output/chapter6/aggregatedprediction_bio.asc")
writeRaster(all.pred.max, file="C:/data/modelling/output/chapter6/aggregatedprediction_max.asc")
writeRaster(all.pred.brt, file="C:/data/modelling/output/chapter6/aggregatedprediction_brt.asc")

#Create maps of uncertainty from grid ----
vari.pred.bio <- data.extent
for (f in list.files("C:/data/modelling/output/chapter6/grids", pattern="1_vari.grd", full.names=T)){
  this.raster <- raster(f)
  this.raster[is.na(this.raster)] <- 0
  vari.pred.bio <- vari.pred.bio + this.raster
  print(paste("...done with", f))
}
vari.pred.max <- data.extent
for (f in list.files("C:/data/modelling/output/chapter6/grids", pattern="2_vari.grd", full.names=T)){
  this.raster <- raster(f)
  this.raster[is.na(this.raster)] <- 0
  vari.pred.max <- vari.pred.max + this.raster
  print(paste("...done with", f))
}
vari.pred.brt <- data.extent
for (f in list.files("C:/data/modelling/output/chapter6/grids", pattern="3_vari.grd", full.names=T)){
  this.raster <- raster(f)
  this.raster[is.na(this.raster)] <- 0
  vari.pred.brt <- vari.pred.brt + this.raster
  print(paste("...done with", f))
}
writeRaster(vari.pred.bio, file="C:/data/modelling/output/chapter6/aggregatedvariance_bio.asc")
writeRaster(vari.pred.max, file="C:/data/modelling/output/chapter6/aggregatedvariance_max.asc")
writeRaster(vari.pred.brt, file="C:/data/modelling/output/chapter6/aggregatedvariance_brt.asc")

writeRaster(vari.pred.bio / all.pred.bio, file="C:/data/modelling/output/chapter6/rel_aggregatedvariance_bio.asc")
writeRaster(vari.pred.max / all.pred.max, file="C:/data/modelling/output/chapter6/rel_aggregatedvariance_max.asc")
writeRaster(vari.pred.brt / all.pred.brt, file="C:/data/modelling/output/chapter6/rel_aggregatedvariance_brt.asc")


#writeRaster(spec.pred.list, file="C:/data/modelling/output/chapter6/aggregatedprediction.grd")
#writeRaster(spec.pred.list, file="C:/data/modelling/output/chapter6/aggregatedvariance.grd")

#Aggregated model eval plots -----
data.files2 <- list.files("G:/HPC/output/finalmodel", pattern="AICc_finalmodels_", full.names=T)
all.mean = data.frame("TaxonKey"=as.factor(NULL),
                      "rangemap"=as.logical(NULL),
                      "beta"=as.numeric(NULL),
                      "bag"=as.numeric(NULL),
                      "BIO.auc"=as.numeric(NULL),
                      "MAX.auc"=as.numeric(NULL),
                      "BRT.auc"=as.numeric(NULL),
                      "BIO.cor"=as.numeric(NULL),
                      "MAX.cor"=as.numeric(NULL),
                      "BRT.cor"=as.numeric(NULL))
all.vari = all.mean
for (s in 1:length(data.files2)){
  TaxonKey = unlist(strsplit(unlist(strsplit(data.files2[s], split="G:/HPC/output/finalmodel/AICc_finalmodels_")), split=".csv"))
  SpecInfo = data.specinfo[data.specinfo$TaxonKey==TaxonKey,]
  
  values <- read.csv(data.files2[s])
  values <- values[,c("TaxonKey", "rangemap", "beta", "bag", "BIO.auc", "MAX.auc", "BRT.auc", "BIO.cor", "MAX.cor", "BRT.cor")]
  values.mean <- apply(values[3:length(values)], FUN=mean, MARGIN=2)
  values.mean <- c(values[1,c("TaxonKey", "rangemap")], values.mean)
  values.vari <- apply(values[3:length(values)], FUN=var, MARGIN=2)
  values.vari <- c(values[1,c("TaxonKey", "rangemap")], values.vari)
  
  all.mean <- rbind(all.mean, values.mean)
  all.vari <- rbind(all.vari, values.vari)
}


values.auc <- melt(all.mean, id.vars=c("TaxonKey", "rangemap"), measure.vars=c("BIO.auc","MAX.auc","BRT.auc"))

data.plot <- values.auc
tiff("C:/Users/LaurensG/Desktop/plots/chapter6/boxplots_aic.tiff", width=16, height=12, unit="cm", res=600)
  ggplot(data = data.plot, aes(variable, value)) + geom_boxplot()
, showLegend=FALSE, width=1, notch=T, outlier.shape=21)

+ 
    xlab ("Modelrun") + ylab ("AUC") +
    theme(axis.text.x=element_text(angle=-45, hjust=0, color="gray40"), axis.text.y=element_text(color="gray40"), 
          axis.title.x=element_text(vjust=-0.5, hjust=0.5), axis.title.y=element_text(vjust=1.5))
dev.off()



# Cluster analysis ----
# Do the experiment for some point clouds.
#   cx, cy: vectors of of x and y of the centers of the point clouds;
#   sd: standard deviation of point cloud around each center;
#   nr_samples: number of points in each point cloud.
# Additional parameters will be passed to the XMeans function.
# Minimum number of clusters is L = 1, maximum H = 10.
# Returns a list with a plot and some info.
library(RWeka)
experiment <- function(cx, cy, sd, nr_samples = 25, ...) {
  df <- data.frame(
    x = rnorm(nr_samples * length(cx), cx, sd),
    y = rnorm(nr_samples * length(cy), cy, sd)
  )
  xmeans_df <- XMeans(df, Weka_control(I = 2, L = 1, H = 10, ...))
  nr_clusters <- xmeans_df$clusterer$numberOfClusters()
  df$cluster_ids <- factor(xmeans_df$class_ids)
  # Here is the result:
  list(
    plot =
      ggplot(df, aes(x = x, y = y, color = cluster_ids)) +
      geom_point(shape=1) +
      scale_colour_hue(l=50),
    info = xmeans_df$clusterer$toString(),
    nr_clusters = nr_clusters
  )
}

experiment_1 <- experiment(cx = c(1, 3, 2), cy = c(1, 2, 3), sd = 0.2)
experiment_1$plot


library(NbClust)

wssplot <- function(data, nc=15, seed=1234){
  wss <- (nrow(data)-1)*sum(apply(data,2,var))
  for (i in 2:nc){
    set.seed(seed)
    wss[i] <- sum(kmeans(data, centers=i)$withinss)}
  plot(1:nc, wss, type="b", xlab="Number of Clusters",
       ylab="Within groups sum of squares")}

cx=c(1, 3, 2)
cy=c(1, 2, 3)
df <- data.frame(
  x = rnorm(100 * length(cx), cx, 0.2),
  y = rnorm(100 * length(cy), cy, 0.2))

wssplot(df)

specinfo <- read.csv("C:/Users/LaurensG/Desktop/specinfo.csv", header=T, na.string="#N/A")
modelinfo <- specinfo[(complete.cases(specinfo[,63:114])),63:114]
wssplot(modelinfo)

clustdata <- modelinfo[,c("baseline.fulld.auc.BIO",
                          "baseline.fulld.cor.BIO",
                          "baseline.clipped.auc.BIO",
                          "baseline.clipped.cor.BIO",
                          "aic.fulld.auc.BIO",
                          "aic.fulld.cor.BIO",
                          "bias.fulld.auc.BIO",
                          "bias.fulld.cor.BIO",
                          "bias.clipped.auc.BIO",
                          "bias.clipped.cor.BIO",
                          "geodist.fulld.auc.BIO",
                          "geodist.fulld.cor.BIO",
                          "finalmodel.clipped.auc.BIO.mean",
                          "finalmodel.clipped.cor.BIO.mean")]

clustdata2 <- modelinfo[,c("baseline.fulld.auc.BRT",
                          "baseline.fulld.cor.BRT",
                          "baseline.clipped.auc.BRT",
                          "baseline.clipped.cor.BRT",
                          "aic.fulld.auc.BRT",
                          "aic.fulld.cor.BRT",
                          "bias.fulld.auc.BRT",
                          "bias.fulld.cor.BRT",
                          "bias.clipped.auc.BRT",
                          "bias.clipped.cor.BRT",
                          "geodist.fulld.auc.BRT",
                          "geodist.fulld.cor.BRT",
                          "finalmodel.clipped.auc.BRT.mean",
                          "finalmodel.clipped.cor.BRT.mean")]

clustdata[,5] = clustdata[,5] - clustdata[,1]
clustdata[,6] = clustdata[,6] - clustdata[,2]
clustdata[,7] = clustdata[,7] - clustdata[,1]
clustdata[,8] = clustdata[,8] - clustdata[,2]
clustdata[,9] = clustdata[,9] - clustdata[,1]
clustdata[,10] = clustdata[,10] - clustdata[,2]
clustdata[,11] = clustdata[,11] - clustdata[,1]
clustdata[,12] = clustdata[,12] - clustdata[,2]
clustdata[,13] = clustdata[,13] - clustdata[,1]
clustdata[,14] = clustdata[,14] - clustdata[,2]

clustdata <- clustdata[,5:14]

clustdata <- modelinfo[,c("baseline.fulld.auc.BIO",
                          "baseline.fulld.cor.BIO",
                          "baseline.clipped.auc.BIO",
                          "baseline.clipped.cor.BIO",
                          "aic.fulld.auc.BIO",
                          "aic.fulld.cor.BIO",
                          "bias.fulld.auc.BIO",
                          "bias.fulld.cor.BIO",
                          "bias.clipped.auc.BIO",
                          "bias.clipped.cor.BIO",
                          "geodist.fulld.auc.BIO",
                          "geodist.fulld.cor.BIO",
                          "finalmodel.clipped.auc.BIO.mean",
                          "finalmodel.clipped.cor.BIO.mean",
                          
                          "baseline.fulld.auc.MAX",
                          "baseline.fulld.cor.MAX",
                          "baseline.clipped.auc.MAX",
                          "baseline.clipped.cor.MAX",
                          "aic.fulld.auc.MAX",
                          "aic.fulld.cor.MAX",
                          "bias.fulld.auc.MAX",
                          "bias.fulld.cor.MAX",
                          "bias.clipped.auc.MAX",
                          "bias.clipped.cor.MAX",
                          "geodist.fulld.auc.MAX",
                          "geodist.fulld.cor.MAX",
                          "finalmodel.clipped.auc.MAX.mean",
                          "finalmodel.clipped.cor.MAX.mean",
                          
                          "baseline.fulld.auc.BRT",
                          "baseline.fulld.cor.BRT",
                          "baseline.clipped.auc.BRT",
                          "baseline.clipped.cor.BRT",
                          "aic.fulld.auc.BRT",
                          "aic.fulld.cor.BRT",
                          "bias.fulld.auc.BRT",
                          "bias.fulld.cor.BRT",
                          "bias.clipped.auc.BRT",
                          "bias.clipped.cor.BRT",
                          "geodist.fulld.auc.BRT",
                          "geodist.fulld.cor.BRT",
                          "finalmodel.clipped.auc.BRT.mean",
                          "finalmodel.clipped.cor.BRT.mean")]

clustdata[,"aic.fulld.auc.BIO"] = clustdata[,"aic.fulld.auc.BIO"] - clustdata[,"baseline.fulld.auc.BIO"]
clustdata[,"aic.fulld.cor.BIO"] = clustdata[,"aic.fulld.cor.BIO"] - clustdata[,"baseline.fulld.cor.BIO"]
clustdata[,"bias.fulld.auc.BIO"] = clustdata[,"bias.fulld.auc.BIO"] - clustdata[,"baseline.fulld.auc.BIO"]
clustdata[,"bias.fulld.cor.BIO"] = clustdata[,"bias.fulld.cor.BIO"] - clustdata[,"baseline.fulld.cor.BIO"]
clustdata[,"bias.clipped.auc.BIO"] = clustdata[,"bias.clipped.auc.BIO"] - clustdata[,"baseline.clipped.auc.BIO"]
clustdata[,"bias.clipped.cor.BIO"] = clustdata[,"bias.clipped.cor.BIO"] - clustdata[,"baseline.clipped.cor.BIO"]
clustdata[,"geodist.fulld.auc.BIO"] = clustdata[,"geodist.fulld.auc.BIO"] - clustdata[,"baseline.clipped.auc.BIO"]
clustdata[,"geodist.fulld.cor.BIO"] = clustdata[,"geodist.fulld.cor.BIO"] - clustdata[,"baseline.clipped.cor.BIO"]
clustdata[,"finalmodel.clipped.auc.BIO.mean"] = clustdata[,"finalmodel.clipped.auc.BIO.mean"] - clustdata[,"baseline.clipped.auc.BIO"]
clustdata[,"finalmodel.clipped.cor.BIO.mean"] = clustdata[,"finalmodel.clipped.cor.BIO.mean"] - clustdata[,"baseline.clipped.cor.BIO"]

clustdata[,"aic.fulld.auc.MAX"] = clustdata[,"aic.fulld.auc.MAX"] - clustdata[,"baseline.fulld.auc.MAX"]
clustdata[,"aic.fulld.cor.MAX"] = clustdata[,"aic.fulld.cor.MAX"] - clustdata[,"baseline.fulld.cor.MAX"]
clustdata[,"bias.fulld.auc.MAX"] = clustdata[,"bias.fulld.auc.MAX"] - clustdata[,"baseline.fulld.auc.MAX"]
clustdata[,"bias.fulld.cor.MAX"] = clustdata[,"bias.fulld.cor.MAX"] - clustdata[,"baseline.fulld.cor.MAX"]
clustdata[,"bias.clipped.auc.MAX"] = clustdata[,"bias.clipped.auc.MAX"] - clustdata[,"baseline.clipped.auc.MAX"]
clustdata[,"bias.clipped.cor.MAX"] = clustdata[,"bias.clipped.cor.MAX"] - clustdata[,"baseline.clipped.cor.MAX"]
clustdata[,"geodist.fulld.auc.MAX"] = clustdata[,"geodist.fulld.auc.MAX"] - clustdata[,"baseline.clipped.auc.MAX"]
clustdata[,"geodist.fulld.cor.MAX"] = clustdata[,"geodist.fulld.cor.MAX"] - clustdata[,"baseline.clipped.cor.MAX"]
clustdata[,"finalmodel.clipped.auc.MAX.mean"] = clustdata[,"finalmodel.clipped.auc.MAX.mean"] - clustdata[,"baseline.clipped.auc.MAX"]
clustdata[,"finalmodel.clipped.cor.MAX.mean"] = clustdata[,"finalmodel.clipped.cor.MAX.mean"] - clustdata[,"baseline.clipped.cor.MAX"]

clustdata[,"aic.fulld.auc.BRT"] = clustdata[,"aic.fulld.auc.BRT"] - clustdata[,"baseline.fulld.auc.BRT"]
clustdata[,"aic.fulld.cor.BRT"] = clustdata[,"aic.fulld.cor.BRT"] - clustdata[,"baseline.fulld.cor.BRT"]
clustdata[,"bias.fulld.auc.BRT"] = clustdata[,"bias.fulld.auc.BRT"] - clustdata[,"baseline.fulld.auc.BRT"]
clustdata[,"bias.fulld.cor.BRT"] = clustdata[,"bias.fulld.cor.BRT"] - clustdata[,"baseline.fulld.cor.BRT"]
clustdata[,"bias.clipped.auc.BRT"] = clustdata[,"bias.clipped.auc.BRT"] - clustdata[,"baseline.clipped.auc.BRT"]
clustdata[,"bias.clipped.cor.BRT"] = clustdata[,"bias.clipped.cor.BRT"] - clustdata[,"baseline.clipped.cor.BRT"]
clustdata[,"geodist.fulld.auc.BRT"] = clustdata[,"geodist.fulld.auc.BRT"] - clustdata[,"baseline.clipped.auc.BRT"]
clustdata[,"geodist.fulld.cor.BRT"] = clustdata[,"geodist.fulld.cor.BRT"] - clustdata[,"baseline.clipped.cor.BRT"]
clustdata[,"finalmodel.clipped.auc.BRT.mean"] = clustdata[,"finalmodel.clipped.auc.BRT.mean"] - clustdata[,"baseline.clipped.auc.BRT"]
clustdata[,"finalmodel.clipped.cor.BRT.mean"] = clustdata[,"finalmodel.clipped.cor.BRT.mean"] - clustdata[,"baseline.clipped.cor.BRT"]

clustdata <- clustdata[,c("aic.fulld.auc.BIO",
                          "aic.fulld.cor.BIO",
                          "bias.fulld.auc.BIO",
                          "bias.fulld.cor.BIO",
                          "bias.clipped.auc.BIO",
                          "bias.clipped.cor.BIO",
                          "geodist.fulld.auc.BIO",
                          "geodist.fulld.cor.BIO",
                          "finalmodel.clipped.auc.BIO.mean",
                          "finalmodel.clipped.cor.BIO.mean",
                          
                          "aic.fulld.auc.MAX",
                          "aic.fulld.cor.MAX",
                          "bias.fulld.auc.MAX",
                          "bias.fulld.cor.MAX",
                          "bias.clipped.auc.MAX",
                          "bias.clipped.cor.MAX",
                          "geodist.fulld.auc.MAX",
                          "geodist.fulld.cor.MAX",
                          "finalmodel.clipped.auc.MAX.mean",
                          "finalmodel.clipped.cor.MAX.mean",
                          
                          "aic.fulld.auc.BRT",
                          "aic.fulld.cor.BRT",
                          "bias.fulld.auc.BRT",
                          "bias.fulld.cor.BRT",
                          "bias.clipped.auc.BRT",
                          "bias.clipped.cor.BRT",
                          "geodist.fulld.auc.BRT",
                          "geodist.fulld.cor.BRT",
                          "finalmodel.clipped.auc.BRT.mean",
                          "finalmodel.clipped.cor.BRT.mean")]

#nc <- NbClust(clustdata2, min.nc=2, max.nc=15, method="kmeans")
#table(nc$Best.n[1,])
#barplot(table(nc$Best.n[1,]),
#        xlab="Numer of Clusters", ylab="Number of Criteria",
#        main="Number of Clusters Chosen by 26 Criteria")
#fit.km$size


#library(HSAUR)
#km <- kmeans(clustdata,3)
#dissE <- daisy(clustdata) 
#dE2   <- dissE^2
#sk2   <- silhouette(km$cl, dE2)
#plot(sk2)
dev.off()

library(cluster)
library(fpc)
fit.km <- kmeans(clustdata, 6, nstart=50)
plotcluster(clustdata, fit.km$cluster)
cluster <- cbind(clustdata, fit.km$cluster)

dev.off()
par(mfrow=c(3,5))
plot(cluster[,c("aic.fulld.auc.BIO", "aic.fulld.cor.BIO")], col=cluster[,"fit.km$cluster"])
plot(cluster[,c("bias.fulld.auc.BIO", "bias.fulld.cor.BIO")], col=cluster[,"fit.km$cluster"])
plot(cluster[,c("bias.clipped.auc.BIO", "bias.clipped.cor.BIO")], col=cluster[,"fit.km$cluster"])
plot(cluster[,c("geodist.fulld.auc.BIO", "geodist.fulld.cor.BIO")], col=cluster[,"fit.km$cluster"])
plot(cluster[,c("finalmodel.clipped.auc.BIO.mean", "finalmodel.clipped.cor.BIO.mean")], col=cluster[,"fit.km$cluster"])

plot(cluster[,c("aic.fulld.auc.MAX", "aic.fulld.cor.MAX")], col=cluster[,"fit.km$cluster"])
plot(cluster[,c("bias.fulld.auc.MAX", "bias.fulld.cor.MAX")], col=cluster[,"fit.km$cluster"])
plot(cluster[,c("bias.clipped.auc.MAX", "bias.clipped.cor.MAX")], col=cluster[,"fit.km$cluster"])
plot(cluster[,c("geodist.fulld.auc.MAX", "geodist.fulld.cor.MAX")], col=cluster[,"fit.km$cluster"])
plot(cluster[,c("finalmodel.clipped.auc.MAX.mean", "finalmodel.clipped.cor.MAX.mean")], col=cluster[,"fit.km$cluster"])

plot(cluster[,c("aic.fulld.auc.BRT", "aic.fulld.cor.BRT")], col=cluster[,"fit.km$cluster"])
plot(cluster[,c("bias.fulld.auc.BRT", "bias.fulld.cor.BRT")], col=cluster[,"fit.km$cluster"])
plot(cluster[,c("bias.clipped.auc.BRT", "bias.clipped.cor.BRT")], col=cluster[,"fit.km$cluster"])
plot(cluster[,c("geodist.fulld.auc.BRT", "geodist.fulld.cor.BRT")], col=cluster[,"fit.km$cluster"])
plot(cluster[,c("finalmodel.clipped.auc.BRT.mean", "finalmodel.clipped.cor.BRT.mean")], col=cluster[,"fit.km$cluster"])





plot(cluster[,3:4], col=cluster[,"fit.km$cluster"])
plot(cluster[,5:6], col=cluster[,11])
plot(cluster[,7:8], col=cluster[,11])
plot(cluster[,9:10], col=cluster[,11])
#plot(cluster[,11:12], col=cluster[,15])
#plot(cluster[,13:14], col=cluster[,15])

#Cell by cell comparison of maps----

#read in all layers
GetPredictions <- function(TaxonKey){
  if (file.exists(paste("G:/HPC/output/baseline/SDMeval_clipped_baseline_", TaxonKey, ".RData", sep=""))){
    load(paste("G:/HPC/output/baseline/SDMeval_clipped_baseline_", TaxonKey, ".RData", sep=""))
    bl.bio <- results[[2]][[1]]; bl.max <- results[[2]][[2]]; bl.brt <- results[[2]][[3]]
  }else{
    bl.bio <- data.extent; bl.max <- data.extent; bl.brt <- data.extent
  }
  if (file.exists(paste("G:/HPC/output/aic/SDMeval_fulld_aic_", TaxonKey, ".RData", sep=""))){
    load(paste("G:/HPC/output/aic/SDMeval_fulld_aic_", TaxonKey, ".RData", sep=""))
    aic.bio <- results[[2]][[1]]; aic.max <- results[[2]][[2]]; aic.brt <- results[[2]][[3]]
  }else{
    aic.bio <- data.extent; aic.max <- data.extent; aic.brt <- data.extent
  }
  if (file.exists(paste("G:/HPC/output/bias/SDMeval_clipped_sb_", TaxonKey, ".RData", sep=""))){
    load(paste("G:/HPC/output/bias/SDMeval_clipped_sb_", TaxonKey, ".RData", sep=""))
    sb.bio <- results[[2]][[1]]; sb.max <- results[[2]][[2]]; sb.brt <- results[[2]][[3]]
  }else{
    sb.bio <- data.extent; sb.max <- data.extent; sb.brt <- data.extent
  }
  if (file.exists(paste("G:/HPC/output/geodist/SDMeval_fulld_geodist_", TaxonKey, ".RData", sep=""))){
    load(paste("G:/HPC/output/geodist/SDMeval_fulld_geodist_", TaxonKey, ".RData", sep=""))
    geo.bio <- results[[2]][[1]]; geo.max <- results[[2]][[2]]; geo.brt <- results[[2]][[3]]
  }else{
    geo.bio <- data.extent; geo.max <- data.extent; geo.brt <- data.extent
  }

  if (file.exists(paste("C:/data/modelling/output/chapter6/grids/", TaxonKey, "_1_mean.grd", sep=""))){
    fm.bio <- raster(paste("C:/data/modelling/output/chapter6/grids/", TaxonKey, "_1_mean.grd", sep=""))
  }else{
    fm.bio <- data.extent
  }
  if (file.exists(paste("C:/data/modelling/output/chapter6/grids/", TaxonKey, "_2_mean.grd", sep=""))){
    fm.max <- raster(paste("C:/data/modelling/output/chapter6/grids/", TaxonKey, "_2_mean.grd", sep=""))
  }
  else{
    fm.max <- data.extent
  }
  if (file.exists(paste("C:/data/modelling/output/chapter6/grids/", TaxonKey, "_3_mean.grd", sep=""))){
    fm.brt <- raster(paste("C:/data/modelling/output/chapter6/grids/", TaxonKey, "_3_mean.grd", sep=""))
  }
  else{
    fm.brt <- data.extent
  }
  allmaps <- stack(bl.bio, bl.max, bl.brt, aic.bio, aic.max, aic.brt, sb.bio, sb.max, sb.brt, geo.bio, geo.max, geo.brt, fm.bio, fm.max, fm.brt)
  names(allmaps) <- c("bl.bio", "bl.max", "bl.brt", "aic.bio", "aic.max", "aic.brt", " sb.bio", "sb.max", "sb.brt", "geo.bio", "geo.max", "geo.brt", "fm.bio", "fm.max", "fm.brt")
  return(allmaps)
}

#comparison methodology
GetMapComparisonStats <- function(baseline.map=NULL, predicted.map=NULL){
  x <- values(baseline.map) ; y <- values(predicted.map)
  complete.difference <- (y - x)
  valid.difference <- complete.difference[!is.na(complete.difference)]
  relevant.difference <- valid.difference[valid.difference > 0.1]
  if (length(relevant.difference) == 0){
    dmean <- 0
    dvar <- 0
    dnum <- 0
    dtotal <- 0
    dacc <- 1
    dmax <- 0
    dmep <- 0
  }else{
    dmean <- mean(relevant.difference)
    dvar <- var(relevant.difference)
    dnum <- length(relevant.difference)
    dtotal <- length(valid.difference)
    dacc <- (dtotal - dnum) / dtotal
    dmax <- max(relevant.difference) 
    dmep <- dmean / dmax
  }
  value.pairs <- cbind(x,y)
  dcor <- cor(value.pairs[complete.cases(value.pairs),1], value.pairs[complete.cases(value.pairs),2])  
  return(c("dmean"=dmean, "dvar"=dvar, "dnum"=dnum, "dtotal"=dtotal, "dacc"=dacc, "dmax"=dmax, "dmep"=dmep, "dcor"=dcor))}

#Get all comparisions
all.comparisons <- array(data=c(0), dim=c(1086,8,12), dimnames=c("species", "metric", "modelrun"))
baseline.map = all.predictions[[model]]
predicted.map = all.predictions[[model+(run*3)]]
s = 1
species.list = c(600015)
for (s in 1:length(data.specinfo[,"TaxonKey"])){
  TaxonKey = data.specinfo[s,"TaxonKey"]
  #600015
  print(paste("now on species: ", TaxonKey))
  all.predictions <- GetPredictions(TaxonKey)
  for (model in 1:3){
    #print(paste("model:", model))
    for (run in 1:4){
      #print(run)
      this.run <- GetMapComparisonStats(all.predictions[[model]], all.predictions[[model+(run*3)]])
      #print(this.run)
      all.comparisons[s,,model+(run*3)-3] <- this.run
    }
  }
}
all.comparisons
save(all.comparisons, file="C:/data/modelling/output/chapter6/all.comparisons.array.RData")

#Merge all the data
spec.info = read.csv("C:/Users/LaurensG/Desktop/specinfo.csv", na.string="#N/A")
spec.info = spec.info[,1:3]
#all names: c("c.bl.bio", "c.bl.max", "c.bl.brt", "c.aic.bio", "c.aic.max", "c.aic.brt", "c.sb.bio", "c.sb.max", "c.sb.brt", "c.geo.bio", "c.geo.max", "c.geo.brt", "c.fm.bio", "c.fm.max", "c.fm.brt")
names = c("c.aic.bio", "c.aic.max", "c.aic.brt", "c.sb.bio", "c.sb.max", "c.sb.brt", "c.geo.bio", "c.geo.max", "c.geo.brt", "c.fm.bio", "c.fm.max", "c.fm.brt")
for(run in 1:dim(all.comparisons)[3]){
  all.comparisons[,,run]
  cluster.data <- all.comparisons[,,run]
  colnames(cluster.data) <- c(paste(names[run], "dmean", sep="."), paste(names[run], "dvar", sep="."), paste(names[run], "dnum", sep="."), paste(names[run], "dtotal", sep="."), paste(names[run], "dacc", sep="."), paste(names[run], "dmax", sep="."), paste(names[run], "dmep", sep="."), paste(names[run], "dcor", sep="."))
  
  spec.info <- cbind(spec.info, cluster.data)
}
write.csv(spec.info, file="C:/data/modelling/output/chapter6/spec.info.comparison.csv")
save(spec.info, file="C:/data/modelling/output/chapter6/spec.info.comparison.RData")

#Calculate clusters based on comparison
library(cluster)
library(NbClust)
library(fpc)
#get spec.info
load(file="C:/data/modelling/output/chapter6/spec.info.comparison.RData")
cluster.info <- spec.info[,1:114]
cluster.data <- spec.info[,115:210]
#Normalize data and select columns
normalize <- function(values){
  return((values- min(values, na.rm=TRUE))/(max(values, na.rm=TRUE) - min(values, na.rm=TRUE)))
}
for (col in 1:length(cluster.data)){
  cluster.data[,col] <- normalize(cluster.data[,col])  
}
cd <- cluster.data[,c("c.fm.bio.dmean", "c.fm.bio.dvar", "c.fm.bio.dacc", "c.fm.bio.dcor", "c.fm.max.dmean", "c.fm.max.dvar", "c.fm.max.dacc", "c.fm.max.dcor", "c.fm.brt.dmean", "c.fm.brt.dvar", "c.fm.brt.dacc", "c.fm.brt.dcor")]

  "c.aic.max.dmean", "c.aic.max.dvar", "c.aic.max.dacc", "c.aic.max.dcor",
                      "c.sb.max.dmean", "c.sb.max.dvar", "c.sb.max.dacc", "c.sb.max.dcor",
                      "c.geo.max.dmean", "c.geo.max.dvar", "c.geo.max.dacc", "c.geo.max.dcor",
                      "c.fm.max.dmean", "c.fm.max.dvar", "c.fm.max.dacc", "c.fm.max.dcor",
                      "c.aic.brt.dmean", "c.aic.brt.dvar", "c.aic.brt.dacc", "c.aic.brt.dcor",
                      "c.sb.brt.dmean", "c.sb.brt.dvar", "c.sb.brt.dacc", "c.sb.brt.dcor",
                      "c.geo.brt.dmean", "c.geo.brt.dvar", "c.geo.brt.dacc", "c.geo.brt.dcor",
                      "c.fm.brt.dmean", "c.fm.brt.dvar", "c.fm.brt.dacc", "c.fm.brt.dcor",
                      "c.geo.bio.dmean", "c.geo.bio.dvar", "c.geo.bio.dacc", "c.geo.bio.dcor",
                      "c.fm.bio.dmean", "c.fm.bio.dvar", "c.fm.bio.dacc", "c.fm.bio.dcor",)]
cd <- cbind(cd, cluster.info)
cd <- cd[complete.cases(cd[,1:12]),]

nc <- NbClust(cd[,1:12], min.nc=3, max.nc=15, method="kmeans")
km <- kmeans(cd[,1:12], 3, nstart=100)
plotcluster(cd[,1:12], km$cluster)
cd <- cbind(cd, km$cluster)

plot(table(cd[,"km$cluster"], cd$dempel))
plot(table(cd[,"km$cluster"], cd$habitat))
plot(table(cd$habitat, cd[,"km$cluster"]), col=cd[,"km$cluster"])

plot(as.numeric(cd$fb_length), as.numeric(cd$fb_weight), col=cd[,"km$cluster"])
plot(cd$meanlat, cd$deep, col=cd[,"km$cluster"])

cor(cd[])

plot(log(cd$fb_length), log(cd$fb_weight), col=cd[,"km$cluster"])
plot(table(cd$habitat, cd[,"km$cluster"]), col=cd[,"km$cluster"])

41] "TaxonKey"                        "TaxonName"                      
[43] "group1"                          "group2"                         
[45] "group3"                          "dempel"                         
[47] "fb_length"                       "fb_weight"                      
[49] "fb_age"                          "shallow"                        
[51] "deep"                            "usual.shallow"                  
[53] "usual.deep"                      "habitat"                        
[55] "fb_lifestyle"                    "fb_latlon"                      
[57] "meanlat"                         "fb_climate"                     
[59] "iucnmap"                         "rsIUCN"                         
[61] "saupmap"                         "rsSAUP"                         
[63] "jaccard"                         "rsBOTH"  

cor(cd[complete.cases(cd[,c("c.fm.bio.dmean", "c.fm.bio.dvar", "c.fm.bio.dacc", "c.fm.bio.dcor", "c.fm.max.dmean", "c.fm.max.dvar", "c.fm.max.dacc", "c.fm.max.dcor", "c.fm.brt.dmean", "c.fm.brt.dvar", "c.fm.brt.dacc", "c.fm.brt.dcor",
                            "fb_length", "fb_weight", "fb_age", "shallow", "deep", "usual.shallow", "usual.deep", "meanlat", "rsIUCN", "rsSAUP")]),c("c.fm.bio.dmean", "c.fm.bio.dvar", "c.fm.bio.dacc", "c.fm.bio.dcor", "fb_length", "fb_weight", "fb_age", "shallow", "deep", "usual.shallow", "usual.deep", "meanlat", "rsIUCN", "rsSAUP")])

