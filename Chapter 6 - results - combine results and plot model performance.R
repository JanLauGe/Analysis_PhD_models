library(foreign)
species.list <- read.csv("C:/Users/LaurensG/Desktop/species.list cross check v22 20151007.csv")

load("C:/data/modelling/data/species/data.obis.RData")
load("C:/data/modelling/data/species/data.gbif.RData")
load("C:/data/modelling/data/species/data.species.RData")


# SPECIES LIST ----
#Counting number of records and joining results to species.list
nobis <- aggregate(data.obis, by=list(data.obis$TaxonKey), FUN=length)[,c("Group.1", "TaxonKey")]
colnames(nobis) <- c("TaxonKey", "nrecs-obis")
ngbif <- aggregate(data.gbif, by=list(data.gbif$TaxonKey), FUN=length)[,c("Group.1", "TaxonKey")]
colnames(ngbif) <- c("TaxonKey", "nrecs-gbif")
nrecords <- aggregate(data.species, by=list(data.species$TaxonKey), FUN=length)[,c("Group.1", "TaxonKey")]
colnames(nrecords) <- c("TaxonKey", "nrecs")

mergetospec <- function(species.list, adtab){
  return(merge(species.list, adtab, by="TaxonKey", all.x=T, all.y=F))}

species.list <- mergetospec(species.list, nobis)
species.list <- mergetospec(species.list, ngbif)
species.list <- mergetospec(species.list, nrecords)


# COMPILE RUN RESULTS ----
runs = c("baseline", "aic", "bias", "geodist")
modes = c("fulld", "clipped")
metrics = c("auc", "cor")
models = c("BIO", "MAX", "BRT")

run <- runs[1]
mode <- modes[1]
metric <- metrics[1]
model <- models[1]

across.runs = data.frame("TaxonKey" = 600005)
for(run in runs){
  path = paste("G:/HPC/output/", run, sep="")
  across.modes = data.frame("TaxonKey" = 600005)
  for(mode in modes){
    across.models = data.frame("TaxonKey" = 600005)
    for(model in models){
      across.species = NULL
      files <- list.files(path, pattern=mode)
      file <- files[1]
      if(length(files) > 0){
        for(file in files){
          load(paste(path, file, sep="/"))
          spec.value <- data.frame("TaxonKey" = gsub("[^0-9]","",file))
          spec.value[,paste(run, mode, "auc", model, sep=".")] <- results[[3]][[which(models %in% model)]]@auc
          spec.value[,paste(run, mode, "cor", model, sep=".")] <- results[[3]][[which(models %in% model)]]@cor
          
          #join results of this species to the results of the other species
          across.species <- rbind(across.species, spec.value)
          print(paste(spec.value[[1]], "done!"))
        }
      }else{
        print(paste(">>> no files for '", mode, "' in: E:/HPC/output/", run, sep=""))
        across.species = data.frame("TaxonKey" = 600005)
      }
      print(head(across.species))
      across.models <- merge(across.models, across.species, all.x=T, all.y=T, by="TaxonKey")
      print(paste("##", model, "done!!"))
    }
    across.modes <- merge(across.modes, across.models, all.x=T, all.y=T, by="TaxonKey")
    print(paste("####", mode, "done!!!"))
  }
  across.runs <- merge(across.runs, across.modes, all.x=T, all.y=T, by="TaxonKey")
  writeLines(paste("####################\n##", run, "done!!!!"))
}


write.csv(across.runs, "C:/Users/LaurensG/Desktop/allruns.csv")

# finalmodel ----
#reading model results and adding them to table
files <- list.files("G:/HPC/output/finalmodel", pattern="*.csv")
file = files[1]
results = data.frame("TaxonKey"=as.character(NULL),
                     "finalmodel.clipped.auc.BIO.mean"=as.numeric(NULL), "finalmodel.clipped.auc.MAX.mean"=as.numeric(NULL), "finalmodel.clipped.auc.BRT.mean"=as.numeric(NULL),
                     "finalmodel.clipped.auc.BIO.var"=as.numeric(NULL), "finalmodel.clipped.auc.MAX.var"=as.numeric(NULL), "finalmodel.clipped.auc.BRT.var"=as.numeric(NULL),
                     "finalmodel.clipped.beta.mean"=as.numeric(NULL), "finalmodel.clipped.bag.mean"=as.numeric(NULL), "finalmodel.clipped.beta.var"=as.numeric(NULL), "finalmodel.clipped.bag.var"=as.numeric(NULL),
                     "finalmodel.clipped.cor.BIO.mean"=as.numeric(NULL), "finalmodel.clipped.cor.MAX.mean"=as.numeric(NULL), "finalmodel.clipped.cor.BRT.mean"=as.numeric(NULL),  
                     "finalmodel.clipped.cor.BIO.var"=as.numeric(NULL), "finalmodel.clipped.cor.MAX.var"=as.numeric(NULL), "finalmodel.clipped.cor.BRT.var"=as.numeric(NULL)) 

for(file in files){
  load(paste("SDMeval_clipped_finalmodel_", species.list$TaxonKey[1], ".RData", sep=""))
  #result <- read.csv(paste("G:/HPC/output/finalmodel/", file, sep=""))
  result <- results
  result <- data.frame("TaxonKey"=result$TaxonKey[1],
                        "finalmodel.clipped.auc.BIO.mean"=mean(result$BIO.auc), "finalmodel.clipped.auc.MAX.mean"=mean(result$MAX.auc), "finalmodel.clipped.auc.BRT.mean"=mean(result$BRT.auc),
                        "finalmodel.clipped.auc.BIO.var"=var(result$BIO.auc), "finalmodel.clipped.auc.MAX.var"=var(result$MAX.auc), "finalmodel.clipped.auc.BRT.var"=var(result$BRT.auc),
                        "finalmodel.clipped.beta.mean"=mean(result$beta), "finalmodel.clipped.bag.mean"=mean(result$bag), "finalmodel.clipped.beta.var"=var(result$beta), "finalmodel.clipped.bag.var"=var(result$bag),
                        "finalmodel.clipped.cor.BIO.mean"=mean(result$BIO.cor), "finalmodel.clipped.cor.MAX.mean"=mean(result$MAX.cor), "finalmodel.clipped.cor.BRT.mean"=mean(result$BRT.cor),  
                        "finalmodel.clipped.cor.BIO.var"=var(result$BIO.cor), "finalmodel.clipped.cor.MAX.var"=var(result$MAX.cor), "finalmodel.clipped.cor.BRT.var"=var(result$BRT.cor)) 
  results <- rbind(results, result)
  print(paste(result$TaxonKey[1],"done"))
}

#MERGE ----
all.results <- merge(across.runs, results, all.x=T, all.y=T, by="TaxonKey") 
other.data <- read.csv("C:/Users/LaurensG/Desktop/full_results_table.csv")[,2:51]
all.info <- merge(other.data, all.results, all.x=T, all.y=T, by="TaxonKey") 

head(all.info)
write.csv(all.info, "C:/Users/LaurensG/Desktop/all.info.csv")

all.info2 <- all.info
all.info <- read.csv("C:/Users/LaurensG/Desktop/specinfo.csv")

#some plots
opar <- par()
boxplot(across.runs[complete.cases(across.runs),c(2,20,32,4,22,34,6,24,36)])
tiff("C:/Users/LaurensG/Desktop/plotall.tiff", height=100, width=100, unit="cm", res=100)
plot(all.info)
dev.off()


# Plots of uncertainty ----
path = paste("G:/HPC/output/finalmodel", sep="")
files <- list.files(path, pattern="*.RData")
file = files[1]
allvar <- data.extent
allrhs <- data.extent
for(file in files){
  load(paste(path, file, sep="/"))
  allvalues <- values(stack(results[[1]][[2]][[2]],results[[2]][[2]][[2]],results[[3]][[2]][[2]],results[[4]][[2]][[2]],results[[5]][[2]][[2]],results[[6]][[2]][[2]],results[[7]][[2]][[2]],results[[8]][[2]][[2]],results[[9]][[2]][[2]],results[[10]][[2]][[2]]))
  
  #mean for species richness
  meanraster <- data.extent
  values(meanraster) <- apply(X=allvalues, MARGIN=1, FUN=mean)
  allrhs <- allrhs + meanraster
  
  #var for spatial variance  
  varraster <- data.extent
  values(varraster) <- apply(X=allvalues, MARGIN=1, FUN=var)
  allvar <- allvar + varraster
  
  print(paste("done with ", file))
}

writeRaster(allvar, "C:/data/modelling/output/finalmodel_spatial_var.asc", format="ascii", overwrite=T)
writeRaster(allvar, "C:/data/modelling/output/finalmodel_spatial_var.grd", overwrite=T)
writeRaster(allrhs, "C:/data/modelling/output/finalmodel_richness.asc", format="ascii")
writeRaster(allrhs, "C:/data/modelling/output/finalmodel_richness.grd")
writeRaster(allvar/allrhs, "C:/data/modelling/output/finalmodel_spatial_var_corrected.asc", format="ascii")
writeRaster(allvar/allrhs, "C:/data/modelling/output/finalmodel_spatial_var_corrected.grd")

###

load("G:/HPC/output/baseline/SDMeval_fulld_baseline_600143.RData")
baseline_BIOCLIM_600143_Thu_alb <- results[[2]][[1]]
baseline_MAXENT_600143_Thu_alb <- results[[2]][[2]]
baseline_BRT_600143_Thu_alb <- results[[2]][[3]]
load("G:/HPC/output/baseline/SDMeval_fulld_baseline_601342.RData")
baseline_BIOCLIM_601342_Ple_pla <- results[[2]][[1]]
baseline_MAXENT_601342_Ple_pla <- results[[2]][[2]]
baseline_BRT_601342_Ple_pla <- results[[2]][[3]]
load("G:/HPC/output/baseline/SDMeval_fulld_baseline_605367.RData")
baseline_BIOCLIM_605367_Epi_are <- results[[2]][[1]]
baseline_MAXENT_605367_Epi_are <- results[[2]][[2]]
baseline_BRT_605367_Epi_are <- results[[2]][[3]]
load("G:/HPC/output/baseline/SDMeval_fulld_baseline_600656.RData")
baseline_BIOCLIM_600656_Cen_fab <- results[[2]][[1]]
baseline_MAXENT_600656_Cen_fab <- results[[2]][[2]]
baseline_BRT_600656_Cen_fab <- results[[2]][[3]]
load("G:/HPC/output/baseline/SDMeval_fulld_baseline_600138.RData")
baseline_BIOCLIM_600138_Som_mic <- results[[2]][[1]]
baseline_MAXENT_600138_Som_mic <- results[[2]][[2]]
baseline_BRT_600138_Som_mic <- results[[2]][[3]]

writeRaster(baseline_BIOCLIM_600143_Thu_alb,"C:/data/modelling/figures/arcGIS maps/mapdata/baseline_BIOCLIM_600143_Thu_alb.asc", format="ascii")
writeRaster(baseline_MAXENT_600143_Thu_alb,"C:/data/modelling/figures/arcGIS maps/mapdata/baseline_MAXENT_600143_Thu_alb.asc", format="ascii")
writeRaster(baseline_BRT_600143_Thu_alb,"C:/data/modelling/figures/arcGIS maps/mapdata/baseline_BRT_600143_Thu_alb.asc", format="ascii")
writeRaster(baseline_BIOCLIM_601342_Ple_pla, "C:/data/modelling/figures/arcGIS maps/mapdata/baseline_BIOCLIM_601342_Ple_pla.asc", format="ascii")
writeRaster(baseline_MAXENT_601342_Ple_pla, "C:/data/modelling/figures/arcGIS maps/mapdata/baseline_MAXENT_601342_Ple_pla.asc", format="ascii")
writeRaster(baseline_BRT_601342_Ple_pla, "C:/data/modelling/figures/arcGIS maps/mapdata/baseline_BRT_601342_Ple_pla.asc", format="ascii")
writeRaster(baseline_BIOCLIM_605367_Epi_are, "C:/data/modelling/figures/arcGIS maps/mapdata/baseline_BIOCLIM_605367_Epi_are.asc", format="ascii")
writeRaster(baseline_MAXENT_605367_Epi_are, "C:/data/modelling/figures/arcGIS maps/mapdata/baseline_MAXENT_605367_Epi_are.asc", format="ascii")
writeRaster(baseline_BRT_605367_Epi_are, "C:/data/modelling/figures/arcGIS maps/mapdata/baseline_BRT_6605367_Epi_are.asc", format="ascii")
writeRaster(baseline_BIOCLIM_600656_Cen_fab, "C:/data/modelling/figures/arcGIS maps/mapdata/baseline_BIOCLIM_600656_Cen_fab.asc", format="ascii")
writeRaster(baseline_MAXENT_600656_Cen_fab, "C:/data/modelling/figures/arcGIS maps/mapdata/baseline_MAXENT_600656_Cen_fab.asc", format="ascii")
writeRaster(baseline_BRT_600656_Cen_fab, "C:/data/modelling/figures/arcGIS maps/mapdata/baseline_BRT_600656_Cen_fab.asc", format="ascii")
writeRaster(baseline_BIOCLIM_600138_Som_mic, "C:/data/modelling/figures/arcGIS maps/mapdata/baseline_BIOCLIM_600138_Som_mic.asc", format="ascii")
writeRaster(baseline_MAXENT_600138_Som_mic, "C:/data/modelling/figures/arcGIS maps/mapdata/baseline_MAXENT_600138_Som_mic.asc", format="ascii")
writeRaster(baseline_BRT_600138_Som_mic, "C:/data/modelling/figures/arcGIS maps/mapdata/baseline_BRT_600138_Som_mic.asc", format="ascii")

600143_Thu_alb, 
601342_Ple_pla, 
605367_Epi_are, 
600656_Cen_fab, 
600138_Som_mic, 


600143 Thunnus albacares
601342 Pleuronectes platessa
605367 Epinephelus areolatus
600656 Centroscyllium fabricii
600138 Somniosus microcephalus


#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#
#exploratory scatterplots----
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#

specinfo <- read.csv("C:/Users/LaurensG/Desktop/specinfo.csv", na.string="#N/A")

#add some transformed variables
specinfo$polarity <- abs(specinfo$meanlat)
specinfo$log.nrecs.obis <- log(specinfo$nrecs.obis)
specinfo$log.nrecs.gbif <- log(specinfo$nrecs.gbif)
specinfo$log.nrecs <- log(specinfo$nrecs)

taxinfo <- specinfo[,c("TaxonKey", "fb_length", "fb_weight", "fb_age", "nrecs.obis", "nrecs.gbif", "nrecs")]
                       # "baseline.fulld.auc.BIO", "baseline.fulld.auc.MAX", "baseline.fulld.auc.BRT", "baseline.fulld.cor.BIO", "baseline.fulld.cor.MAX", "baseline.fulld.cor.BRT")]
ecoinfo <- specinfo[,c("TaxonKey", "shallow", "deep", "usual.shallow", "usual.deep", "nrecs.obis", "nrecs.gbif", "nrecs")]
                      #, "nrecs.obis", "nrecs.gbif", "nrecs", "baseline.fulld.auc.BIO", "baseline.fulld.auc.MAX", "baseline.fulld.auc.BRT", "baseline.fulld.cor.BIO", "baseline.fulld.cor.MAX", "baseline.fulld.cor.BRT")]
geoinfo <- specinfo[,c("TaxonKey", "meanlat", "polarity", "rsIUCN", "rsSAUP", "rsBOTH", "nrecs.obis", "nrecs.gbif", "nrecs")]
                       #, "nrecs.obis", "nrecs.gbif", "nrecs", "baseline.fulld.auc.BIO", "baseline.fulld.auc.MAX", "baseline.fulld.auc.BRT", "baseline.fulld.cor.BIO", "baseline.fulld.cor.MAX", "baseline.fulld.cor.BRT")]
useinfo <- specinfo[,c("TaxonKey", "Importance", "MainCatchingMethod", "nrecs.obis", "nrecs.gbif", "nrecs")]
                       #, "nrecs.obis", "nrecs.gbif", "nrecs", "baseline.fulld.auc.BIO", "baseline.fulld.auc.MAX", "baseline.fulld.auc.BRT", "baseline.fulld.cor.BIO", "baseline.fulld.cor.MAX", "baseline.fulld.cor.BRT")]


library(GGally)
tiff("C:/Users/LaurensG/Desktop/plots/chapter2/taxinfo.tiff", res=600, width=40, height=35, units="cm")
ggpairs(taxinfo[,2:7], title="Exploration for Taxonomic Bias",  columnLabels=c("Length", "Weight", "Age", "OBIS Records", "GBIF Records", "Total Unique Records"))
dev.off()

tiff("C:/Users/LaurensG/Desktop/plots/chapter2/ecoinfo.tiff", res=300, width=30, height=30, units="cm")
ggpairs(ecoinfo[,2:8], title="Exploration for Taxonomic Bias",  columnLabels=c("Shallow", "Deep", "Usual Shallow", "Usual Deep", "OBIS Records", "GBIF Records", "Total Unique Records"))
dev.off()

tiff("C:/Users/LaurensG/Desktop/plots/chapter2/geoinfo.tiff", res=300, width=30, height=30, units="cm")
ggpairs(geoinfo[,2:9], title="Exploration for Spatial Bias",  columnLabels=c("Mean Latitude", "Distance from Equator", "Range Size IUCN", "Range Size SAUP", "Range Size Combined", "OBIS Records", "GBIF Records", "Total Unique Records"))
dev.off()

ggpairs(useinfo[,c(2,3,6)])
#ggpairs(useinfo[,2:7])




tiff("C:/Users/LaurensG/Desktop/taxinfo.tiff", res=300, width=30, height=30, units="cm")
pairs(taxinfo)
dev.off()
tiff("C:/Users/LaurensG/Desktop/ecoinfo.tiff", res=300, width=30, height=30, units="cm")
pairs(ecoinfo)
dev.off()
tiff("C:/Users/LaurensG/Desktop/geoinfo.tiff", res=300, width=30, height=30, units="cm")
pairs(geoinfo)
dev.off()
tiff("C:/Users/LaurensG/Desktop/useinfo.tiff", res=300, width=30, height=30, units="cm")
pairs(useinfo)
dev.off()

plot(specinfo$nrecs, specinfo$baseline.fulld.auc.MAX)

###
# DEFINE FUNCTIONS ----

library(ggplot2)
library(grid)

singlescatterplot <- function(data=data.frame(), col1=string(), col2=string(), lab1=as.character(NULL), lab2=as.character(NULL), pointcolour=NULL){
  corvalout <- cor(data[,c(col1, col2)], use="complete.obs")
  corsigout <- cor.test(x=data[,which(colnames(data)==col1)], y=data[,which(colnames(data)==col2)])$p.value
  
  if(corsigout > 0.05){
    signifout <- c("ns")
  }else if(corsigout <= 0.05 & corsigout > 0.01){
    signifout <- c("*")    
  }else if(corsigout <= 0.01 & corsigout > 0.001){
    signifout <- c("**")    
  }else if(corsigout <= 0.001){
    signifout <- c("***")
  }else{
    print("oops")
  }
  
  textadout <- paste("n: ", "\t", dim(na.omit(data[,c(col1,col2)]))[1], "\n",
            "correlation: ", "\t", round(corvalout[2], 3), "\n",
            "p-value: ", "\t", round(corsigout, 3), "  ", signifout, "\n  ", sep="")
  
  if(missing(lab1)){
    lab1=element_blank()
  }
  if(missing(lab2)){
    lab2=element_blank()
  }
  
  if(!missing(pointcolour)){
    ggplotout <- 
      ggplot(data, aes_string(col1, col2, col=pointcolour))+ geom_point(size=0.5) + labs(x=lab1, y=lab2) +
      ggtitle(textadout) + theme(plot.title=element_text(size=8, hjust=1), axis.title=element_text(size=8), axis.text=element_blank()) 
  }else if(missing(pointcolour)){
    ggplotout <- 
      ggplot(data, aes_string(col1, col2)) + geom_point(size=0.5) + labs(x=lab1, y=lab2) +
      ggtitle(textadout) + theme(plot.title=element_text(size=8, hjust=1), axis.title=element_text(size=8), axis.text=element_blank()) 
  }
  return(ggplotout)
}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

# MAKE MULTIPLOT SCATTERPLOTS -----

#~Manual plots----

plot11 <- singlescatterplot(specinfo, "fb_length", "nrecs.obis", lab2="OBIS records")
plot12 <- singlescatterplot(specinfo, "fb_weight", "nrecs.obis")
plot13 <- singlescatterplot(specinfo, "fb_age", "nrecs.obis")
plot14 <- singlescatterplot(specinfo, "shallow", "nrecs.obis")
plot15 <- singlescatterplot(specinfo, "deep", "nrecs.obis")
plot16 <- singlescatterplot(specinfo, "meanlat", "nrecs.obis")
plot17 <- singlescatterplot(specinfo, "polarity", "nrecs.obis")
#
plot21 <- singlescatterplot(specinfo, "fb_length", "nrecs.gbif", lab2="GBIF records")
plot22 <- singlescatterplot(specinfo, "fb_weight", "nrecs.gbif")
plot23 <- singlescatterplot(specinfo, "fb_age", "nrecs.gbif")
plot24 <- singlescatterplot(specinfo, "shallow", "nrecs.gbif")
plot25 <- singlescatterplot(specinfo, "deep", "nrecs.gbif")
plot26 <- singlescatterplot(specinfo, "meanlat", "nrecs.gbif")
plot27 <- singlescatterplot(specinfo, "polarity", "nrecs.gbif")
#
plot31 <- singlescatterplot(specinfo, "fb_length", "nrecs", "Length", "Unique records")
plot32 <- singlescatterplot(specinfo, "fb_weight", "nrecs", "Weight")
plot33 <- singlescatterplot(specinfo, "fb_age", "nrecs", "Age")
plot34 <- singlescatterplot(specinfo, "shallow", "nrecs", "Upper depth limit")
plot35 <- singlescatterplot(specinfo, "deep", "nrecs", "Lower depth limit")
plot36 <- singlescatterplot(specinfo, "meanlat", "nrecs", "Mean latitude")
plot37 <- singlescatterplot(specinfo, "polarity", "nrecs", "Distance from equator")

tiff("C:/Users/LaurensG/Desktop/PP cor plots test.tiff", height=14, width=24, units="cm", res=600)
multiplot(plot11, plot21, plot31, plot12, plot22, plot32, plot13, plot23, plot33, 
          plot14, plot24, plot34, plot15, plot25, plot35, plot16, plot26, plot36,
          plot17, plot27, plot37, cols=7)
dev.off()
#
#

#
#
#
#
#
for(row in 1:19){
  assign(paste("plot", row, sep=""), singlescatterplot(specinfo, as.character(plotinfo[row,1]), as.character(plotinfo[row,2]), as.character(plotinfo[row,3]), as.character(plotinfo[row,4])))
  #cat ("Press [enter] to continue")
  #line <- readline()
}
multiplot(plot1, plot2, plot3, plot4,
          plot5, plot6, plot7, plot8,
          plot9, plot10, plot11, plot12,
          plot13, plot14, plot15, plot16,
          plot17, plot18, plot19, cols=3)
dev.off()


#' categorical variables to plot:
#' "group1", "group2", "group3", "dempel", "habitat", "fb_lifestyle", "fb_climate", "Importance", "MainCatchingMethod"
#' "Taxonomic group", "Taxonomic group", "Taxonomic group", "Zone", "Habitat", "Lifestyle", "Climate", "Economic importance", "Main catching method"
#' 
#' continuous variables to plot:
#' "fb_length", "fb_weight", "fb_age", "shallow", "deep", "usual.shallow", "usual.deep", "meanlat", "polarity", "rsIUCN", "rsSAUP", "jaccard", "rsBOTH", "nrecs.obis", "nrecs.gbif", "nrecs", "log.nrecs.obis", "log.nrecs.gbif", "log.nrecs"
#' "Length", "Weight", "Age", "upper depth limit", "lower depth limit", "upper depth interval", "lower depth interval", "Mean latitude", "Distance from equator", "IUCN range size", "SAUP range size", "Jaccard index", "Combined range size", "Number of OBIS records", "Number of GBIF records", "Number of records", "Log number of OBIS records", "Log number of GBIF records", "Log number of records"

plotinfo <- data.frame(
  colname1=c("group1", "group2", "group3", "dempel", "habitat", "fb_lifestyle", "fb_climate", "Importance", "MainCatchingMethod"),
  colname2=rep("baseline.fulld.auc.MAX", 9),
  varname1=c("Taxonomic group", "Taxonomic group", "Taxonomic group", "Zone", "Habitat", "Lifestyle", "Climate", "Economic importance", "Main catching method"),
  varname2=rep("AUC", 9))

plotinfo <- data.frame(
  colname1=c("fb_length", "fb_weight", "fb_age", "shallow", "deep", "usual.shallow", "usual.deep", "meanlat", "polarity", "rsIUCN", "rsSAUP", "jaccard", "rsBOTH", "nrecs.obis", "nrecs.gbif", "nrecs", "log.nrecs.obis", "log.nrecs.gbif", "log.nrecs"),
  colname2=rep("nrecs", 19),
  varname1=c("Length", "Weight", "Age", "Upper depth limit", "Lower depth limit", "Upper depth interval", "Lower depth interval", "Mean latitude", "Distance from equator", "IUCN range size", "SAUP range size", "Jaccard index", "Combined range size", "Number of OBIS records", "Number of GBIF records", "Number of records", "Log number of OBIS records", "Log number of GBIF records", "Log number of records"),
  varname2=rep("records", 19))

tiff("C:/Users/LaurensG/Desktop/PP cor.tiff", height=40, width=40, units="cm", res=300)
for(row in 1:19){
  assign(paste("plot", row, sep=""), singlescatterplot(specinfo, as.character(plotinfo[row,1]), as.character(plotinfo[row,2]), as.character(plotinfo[row,3]), as.character(plotinfo[row,4]), "dempel"))
  #cat ("Press [enter] to continue")
  #line <- readline()
}
multiplot(plot1, plot2, plot3, plot4,
          plot5, plot6, plot7, plot8,
          plot9, plot10, plot11, plot12,
          plot13, plot14, plot15, plot16,
          plot17, plot18, plot19, cols=4)
dev.off()


#~Correlation of predictive performance with ecological traits----
plotinfo <- data.frame(
  colname1=c("group1", "group2", "group3", "dempel", "habitat", "fb_lifestyle", "fb_climate", "Importance", "MainCatchingMethod"),
  colname2=rep("baseline.fulld.auc.MAX", 9),
  varname1=c("Taxonomic group", "Taxonomic group", "Taxonomic group", "Zone", "Habitat", "Lifestyle", "Climate", "Economic importance", "Main catching method"),
  varname2=rep("AUC", 9))

plotinfo <- data.frame(
  colname1=c("fb_length", "fb_weight", "fb_age", "shallow", "deep", "usual.shallow", "usual.deep", "meanlat", "polarity", "rsIUCN", "rsSAUP", "jaccard", "rsBOTH", "nrecs.obis", "nrecs.gbif", "nrecs", "log.nrecs.obis", "log.nrecs.gbif", "log.nrecs"),
  colname2=rep("baseline.fulld.auc.MAX", 19),
  varname1=c("Length", "Weight", "Age", "Upper depth limit", "Lower depth limit", "Upper depth interval", "Lower depth interval", "Mean latitude", "Distance from equator", "IUCN range size", "SAUP range size", "Jaccard index", "Combined range size", "Number of OBIS records", "Number of GBIF records", "Number of records", "Log number of OBIS records", "Log number of GBIF records", "Log number of records"),
  varname2=rep("AUC", 19))

#opar <- par()
#par(mfrow=c(5,5))

tiff("C:/Users/LaurensG/Desktop/PP cor.tiff", height=40, width=40, units="cm", res=300)
for(row in 1:19){
  assign(paste("plot", row, sep=""), singlescatterplot(specinfo, as.character(plotinfo[row,1]), as.character(plotinfo[row,2]), as.character(plotinfo[row,3]), as.character(plotinfo[row,4]), "dempel"))
  #cat ("Press [enter] to continue")
  #line <- readline()
}
multiplot(plot1, plot2, plot3, plot4,
          plot5, plot6, plot7, plot8,
          plot9, plot10, plot11, plot12,
          plot13, plot14, plot15, plot16,
          plot17, plot18, plot19, cols=4)
dev.off()

##
##
##



row=1
singlescatterplot(specinfo, as.character(plotinfo[row,1]), as.character(plotinfo[row,2]), as.character(plotinfo[row,3]), as.character(plotinfo[row,4]))
row=2
singlescatterplot(specinfo, as.character(plotinfo[row,1]), as.character(plotinfo[row,2]), as.character(plotinfo[row,3]), as.character(plotinfo[row,4]))
row=3
singlescatterplot(specinfo, as.character(plotinfo[row,1]), as.character(plotinfo[row,2]), as.character(plotinfo[row,3]), as.character(plotinfo[row,4]))
row=4
singlescatterplot(specinfo, as.character(plotinfo[row,1]), as.character(plotinfo[row,2]), as.character(plotinfo[row,3]), as.character(plotinfo[row,4]))
row=5
singlescatterplot(specinfo, as.character(plotinfo[row,1]), as.character(plotinfo[row,2]), as.character(plotinfo[row,3]), as.character(plotinfo[row,4]))
row=6
singlescatterplot(specinfo, as.character(plotinfo[row,1]), as.character(plotinfo[row,2]), as.character(plotinfo[row,3]), as.character(plotinfo[row,4]))
row=7
singlescatterplot(specinfo, as.character(plotinfo[row,1]), as.character(plotinfo[row,2]), as.character(plotinfo[row,3]), as.character(plotinfo[row,4]))
row=8
singlescatterplot(specinfo, as.character(plotinfo[row,1]), as.character(plotinfo[row,2]), as.character(plotinfo[row,3]), as.character(plotinfo[row,4]))
row=9
singlescatterplot(specinfo, as.character(plotinfo[row,1]), as.character(plotinfo[row,2]), as.character(plotinfo[row,3]), as.character(plotinfo[row,4]))
row=10
singlescatterplot(specinfo, as.character(plotinfo[row,1]), as.character(plotinfo[row,2]), as.character(plotinfo[row,3]), as.character(plotinfo[row,4]))
row=11
singlescatterplot(specinfo, as.character(plotinfo[row,1]), as.character(plotinfo[row,2]), as.character(plotinfo[row,3]), as.character(plotinfo[row,4]))
row=12
singlescatterplot(specinfo, as.character(plotinfo[row,1]), as.character(plotinfo[row,2]), as.character(plotinfo[row,3]), as.character(plotinfo[row,4]))
row=13
singlescatterplot(specinfo, as.character(plotinfo[row,1]), as.character(plotinfo[row,2]), as.character(plotinfo[row,3]), as.character(plotinfo[row,4]))
row=14
singlescatterplot(specinfo, as.character(plotinfo[row,1]), as.character(plotinfo[row,2]), as.character(plotinfo[row,3]), as.character(plotinfo[row,4]))
row=15
singlescatterplot(specinfo, as.character(plotinfo[row,1]), as.character(plotinfo[row,2]), as.character(plotinfo[row,3]), as.character(plotinfo[row,4]))





#species records
ggplot(specinfo, aes(x=nrecs, y=baseline.fulld.auc.MAX)) + geom_point()
cor(specinfo[,c("nrecs", "baseline.fulld.auc.MAX")], use="complete.obs")
cor.test(x=specinfo$nrecs, y=specinfo$baseline.fulld.auc.MAX)

ggplot(specinfo, aes(x=nrecs, y=baseline.fulld.cor.MAX)) + geom_point()
cor(specinfo[,c("nrecs", "baseline.fulld.cor.MAX")], use="complete.obs")
cor.test(x=specinfo$nrecs, y=specinfo$baseline.fulld.auc.MAX)

#length
ggplot(specinfo, aes(x=fb_length, y=baseline.fulld.cor.MAX)) + geom_point()
cor(specinfo[,c("nrecs", "baseline.fulld.cor.MAX")], use="complete.obs")
cor.test(x=specinfo$nrecs, y=specinfo$baseline.fulld.auc.MAX)

#weight
ggplot(specinfo, aes(x=log(fb_weight), y=baseline.fulld.cor.MAX)) + geom_point()
cor(specinfo[,c("fb_weight", "baseline.fulld.cor.MAX")], use="complete.obs")
cor.test(x=specinfo$fb_weight, y=specinfo$baseline.fulld.auc.MAX)


###
###
par(mfrow=c(10,10))

tiff("C:/Users/LaurensG/Desktop/testplot2.tiff", res=600, width=30, height=30, units="cm")
plot(specinfo[,c("group1", "group2", "group3", "dempel", "fb_length", "fb_weight", "fb_age", "shallow", "deep", "habitat", "fb_lifestyle", "meanlat", "fb_climate", "baseline.fulld.auc.BIO", "baseline.fulld.auc.MAX", "baseline.fulld.auc.BRT", "baseline.fulld.cor.BIO", "baseline.fulld.cor.MAX", "baseline.fulld.cor.BRT")], )
dev.off()

specinfo2 <- specinfo[complete.cases(specinfo[,c("group3", "baseline.fulld.auc.MAX")]),c("group3", "baseline.fulld.auc.MAX")]
specinfo2$group3 <- reorder(specinfo2$group3, specinfo2$baseline.fulld.auc.MAX, mean)
par(mar=c(5.1,10,4.1,2.1))
boxplot(specinfo2$baseline.fulld.auc.MAX ~ specinfo2$group3, horizontal=T, las=2)

#dempel
boxplot(specinfo$baseline.fulld.auc.MAX ~ specinfo$dempel)

tiff("C:/Users/LaurensG/Desktop/chapter 6 plot.tiff", res=300, width=15, height=15, units="cm")
par(mfrow = c(2, 2))
boxplot(results[,2:4], names=c("BIO", "MAX", "BRT"), main="mean AUC")
boxplot(results[,5:7], names=c("BIO", "MAX", "BRT"), main="var AUC")
boxplot(results[,12:14], names=c("BIO", "MAX", "BRT"), main="mean COR")
boxplot(results[,15:17], names=c("BIO", "MAX", "BRT"), main="var COR")
dev.off()


species.list3 <- mergetospec(species.list2, results)
boxplot(species.list3[,c(51, 53, 55, 63, 65, 67, 57, 59, 61, 69, 70, 71)])

species.list4 <- species.list3[,c(51, 53, 55, 63, 65, 67, 57, 59, 61, 69, 70, 71)]
species.list4 <- species.list4[complete.cases(species.list4),]

boxplot(species.list4)
#save output


write.csv(species.list, file="C:/Users/LaurensG/Desktop/full_results_table.csv")
