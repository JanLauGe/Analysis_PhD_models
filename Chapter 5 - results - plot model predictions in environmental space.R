require(dismo)
require(ggplot2)
require(gridExtra)
require(gtable)

#*#*#*#*#*#*#*#*#
# Data prep -----

# load data
#data.envi <- stack(data.extent, data.envi[[1:7]])
species.list <- read.csv("C:/Users/LaurensG/Desktop/species.list.csv", header=T)
data.extent <- raster("C:/data/modelling/data/data.extent.mod.grd")
data.envi <- mask(stack("C:/data/modelling/data/envi/data.envi.grd"), data.extent)
load("C:/data/modelling/data/data.species.RData")

# set species
#species="Pleuronectes platessa"
#species="Thunnus albacares"

for(s in c("Pleuronectes platessa", "Thunnus albacares", "Epinephelus areolatus", "Centroscyllium fabricii", "Somniosus microcephalus")){ 
  species <- s

p <- data.species[data.species$TaxonName==species,c(3,4)]
a <- data.species[data.species$set=="background",3:4]
x <- as.data.frame(na.omit(extract(data.envi[[1:2]],p)))
#y <- as.data.frame(na.omit(extract(data.envi[[1:2]],a)))
y <- as.data.frame(values(data.envi[[1:2]]))[complete.cases(as.data.frame(values(data.envi[[1:2]]))),]

#normalise envidata
normalise <- function(x, y){
  normx=NULL
  for(col in 1:length(colnames(x))){
    mini = min(min(x[,col]), min(y[,col]))
    maxi = max(max(x[,col]), max(y[,col]))
    normcol <- (x[,col]-mini)/(maxi-mini)
    normx <- cbind(normx, normcol)}
  colnames(normx) <- colnames(x)
  return(normx)}

x.norm <- as.data.frame(normalise(x,y))
y.norm <- as.data.frame(normalise(y,y))

#*#*#*#*#*#*#*#*#*#*#
# Envidata plot -----

# calibrate
  enviplot <- function(x, y, x.norm, y.norm, species){
    grid3 <- ggplot(y, aes(depthmean, sstanmean)) + 
    geom_point(colour="lightgrey", size=1) + geom_point(data=x, color="black", size=1) +
    labs(x="Depth [m]", y="Sea Surface Temperature [°C]") +
    #scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1), labels=c("0", "0.25", "0.5", "0.75", "1"), expand=c(0,0)) +
    #scale_x_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1), labels=c("0", "0.25", "0.5", "0.75", "1"), expand=c(0,0)) +
    theme(plot.margin=unit(c(2,2,0,0),"mm"), 
            legend.position="none", 
            panel.grid.major=element_blank(), 
            panel.grid.minor=element_blank(),
            panel.background=element_blank(),
            panel.border=element_blank(), axis.line=element_line(color="grey"))

  grid1 <- ggplot(y.norm, aes(depthmean)) + 
    geom_histogram(data=y.norm, aes(depthmean, y=..density..), colour="darkgrey", fill="lightgrey", binwidth=0.025) + 
    geom_histogram(data=x.norm, aes(depthmean, y=..density..*0.25), colour="darkgrey", fill="black", binwidth=0.025) +
    scale_y_continuous(expand=c(0,0), labels=NULL) +
    scale_x_continuous(expand=c(0.055,-0.0125), labels=NULL) +
    labs(x=NULL, y="proportion of records") +
    theme(plot.margin=unit(c(2,2,2,0),"mm"), 
          legend.position="none", 
          axis.ticks=element_line(colour="white"), 
          axis.text=element_text(colour="white"), 
          axis.title.x=element_blank(), 
          panel.grid.major=element_blank(), 
          panel.grid.minor=element_blank(),
          panel.background=element_blank(),
          panel.border=element_blank(), axis.line=element_blank())#element_line(color="grey"))

  grid4 <- ggplot(y.norm, aes(sstanmean)) + 
    geom_histogram(data=y.norm, aes(sstanmean, y=..density..), colour="darkgrey", fill="grey", binwidth=0.025) + 
    geom_histogram(data=x.norm, aes(sstanmean, y=..density..*0.25), colour="darkgrey", fill="black", binwidth=0.025) +
    scale_y_continuous(expand=c(0,0), NULL) +
    scale_x_continuous(expand=c(-0.0125,-0.0125), labels=NULL) +
    labs(x=NULL, y="proportion of records") + coord_flip() + 
    theme(plot.margin=unit(c(2,2,0,2),"mm"), 
          legend.position="none", 
          axis.ticks=element_line(colour="white"), 
          axis.text=element_text(colour="white"), 
          axis.title.y=element_blank(), 
          panel.grid.major=element_blank(), 
          panel.grid.minor=element_blank(),
          panel.background=element_blank(),
          panel.border=element_blank(), axis.line=element_line(color="grey"))

  grid2 <- grid.text(paste(unlist(strsplit(species, split=" ")),"\n",unlist(strsplit(species, split=" "))[2], "\n\n", length(x[,1]), " records", sep=""),
                       gp=gpar(fontsize=12), just="left", draw=F, x=0.1, y=0.35)

  #Arrange the four charts
  return(grid.arrange(grid1, grid2, grid3, grid4, ncol=2, nrow=2, widths=c(4, 1.4), heights=c(1.4, 4)))}

# plot graph
  #tiff(paste("C:/Users/LaurensG/Desktop/plots/ggplot_", species, "_1_envispace.tif", sep=""), width=16, height=16, unit="cm", res=600)
  #enviplot(x, y, x.norm, y.norm, species)
  #dev.off()


#*#*#*#*#*#*#*#*#*#
# Model plots ----

# calibrate
  modelplot <- function(pred, species, x, x.norm, y, y.norm, model){
  # prepare data
    #z <- as.data.frame(cbind(y, "z"=na.omit(extract(pred,a))))
    z <- as.data.frame(cbind(y, "z"=values(pred)[complete.cases(values(pred))]))
    z <- z[order(z$z),]
    z.norm <- as.data.frame(normalise(z,z))
    z.norm <- z.norm[order(z.norm$z),]
    z.norm <- cbind(z.norm, "depthbins"=as.factor(cut(z.norm$depthmean, breaks=seq(0, 1, 0.025), labels=seq(0.025, 1, 0.025))), "sstbins"=as.factor(cut(z.norm$sstanmean, breaks=seq(0, 1, 0.025), labels=seq(0.025, 1, 0.025))))
    depthbin <- aggregate(z.norm$z, by=list(z.norm$depthbins), FUN="max")
    depthbin$Group.1 <- as.character(depthbin$Group.1)
    depthbin <- rbind(depthbin, cbind("Group.1"="0.925", "x"=mean(max(depthbin[depthbin$Group.1=="0.9",2]), max(depthbin[depthbin$Group.1=="0.95",2]))))
    depthbin <- depthbin[order(depthbin$Group.1),]
    depthbin <- data.frame("bin"=as.integer(as.factor(depthbin$Group.1)), "max"=as.numeric(depthbin$x))
    sstbin <- aggregate(z.norm$z, by=list(z.norm$sstbin), FUN="max")
    sstbin <- data.frame("bin"=as.integer(sstbin$Group.1), "max"=sstbin$x)
  
    grid3 <- ggplot(z, aes(depthmean, sstanmean)) + 
      geom_point(aes(depthmean, sstanmean, col=z), size=1) +
      scale_colour_gradientn(colours = c("#0000FFFF", "#3399FFFF", "#CCFFFFFF", "#CCFF00FF", "#FFFF00FF", "#FF0000FF")) +
      labs(x="Depth [m]", y="Sea Surface Temperature [°C]") +
      #scale_colour_gradientn(colours=rev(rainbow(7))) +
      #scale_colour_gradientn(low="lightgrey", high="black", space="Lab") +
      #scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1), labels=c("0", "0.25", "0.5", "0.75", "1"), expand=c(0,0)) +
      #scale_x_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1), labels=c("0", "0.25", "0.5", "0.75", "1"), expand=c(0,0)) +
      theme(plot.margin=unit(c(2,2,0,0),"mm"), 
            legend.position="none", 
            panel.grid.major=element_blank(), 
            panel.grid.minor=element_blank(),
            panel.background=element_blank(),
            panel.border=element_blank(), axis.line=element_line(color="grey"))
    
    grid1 <- ggplot(depthbin, aes(x=bin, y=max, col=max)) + geom_line(size=1) + 
      scale_colour_gradientn(colours = c("#0000FFFF", "#3399FFFF", "#CCFFFFFF", "#CCFF00FF", "#FFFF00FF", "#FF0000FF")) +
      #scale_colour_gradient(low="white", high="black", space="Lab", limits=c(0,1)) +
      scale_y_continuous(limits=c(0,1), breaks=c(0, 0.5, 1), labels=c("0", "0.5", "1"), expand=c(0,0)) +
      scale_x_continuous(expand=c(-0,-0), labels=NULL) +
      labs(x=NULL, y="RHS")+
      theme(plot.margin=unit(c(2,2,2,2),"mm"), 
            legend.position="none", 
            axis.ticks.x=element_blank(), 
            axis.text.x=element_blank(), 
            axis.title.x=element_blank(), 
            panel.grid.major=element_blank(), 
            panel.grid.minor=element_blank(),
            panel.background=element_blank(),
            panel.border=element_blank(), axis.line=element_line(color="grey"))
  
    grid4 <- ggplot(sstbin, aes(x=bin, y=max, col=max)) + geom_line(size=1) +
      scale_y_continuous(limits=c(0,1), breaks=c(0, 0.5, 1), labels=c("0", "0.5", "1"), expand=c(0,0)) +
      scale_x_continuous(expand=c(-0,-0), labels=NULL) +
      scale_colour_gradientn(colours = c("#0000FFFF", "#3399FFFF", "#CCFFFFFF", "#CCFF00FF", "#FFFF00FF", "#FF0000FF")) +
      #scale_colour_gradient(low="white", high="black", space="Lab", limits=c(0,1)) +
      labs(x=NULL, y="RHS") + coord_flip() + 
      theme(plot.margin=unit(c(2,2,0,2),"mm"), 
            legend.position="none", 
            axis.ticks.y=element_blank(), 
            axis.text.y=element_blank(), 
            axis.title.y=element_blank(), 
            panel.grid.major=element_blank(), 
            panel.grid.minor=element_blank(),
            panel.background=element_blank(),
            panel.border=element_blank(), axis.line=element_line(color="grey"))
  
    grid2 <- grid.text(paste(model, "\n\n", unlist(strsplit(species, split=" ")),"\n",unlist(strsplit(species, split=" "))[2], "\n\n", length(x[,1]), " records", sep=""),
                       gp=gpar(fontsize=12), just="left", draw=F, x=0.1, y=0.5)
    
    #Arrange the four charts
    return(
      grid.arrange(grid1, grid2, grid3, grid4, ncol=2, nrow=2, widths=c(4, 1.4), heights=c(1.4, 4))
    )}

# Get model results ----
  #load(paste("C:/data/modelling/output/baseline/SDMeval_fulld_", species.list[species.list$TaxonName==species, "TaxonKey"], ".RData", sep=""))
  #load(paste("C:/data/modelling/output/bias1/SDMeval_fulld_sb1_", species.list[species.list$TaxonName==species, "TaxonKey"], ".RData", sep=""))
  load(paste("C:/data/modelling/output/geodist2/SDMeval_geodist_", species.list[species.list$TaxonName==species, "TaxonKey"], ".RData", sep=""))

# # BIOCLIM
#   # run bioclim model
#     bc <- bioclim(x)
#     pred.bc <- predict(bc, x=data.envi[[1:2]])
#   # plot graph
#     tiff(paste("C:/Users/LaurensG/Desktop/plots/ggplot_", species, "_2_bioclim.tif", sep=""), width=16, height=16, unit="cm", res=600)
#     modelplot(pred.bc, species=species, x, x.norm, y, y.norm, model="BIOCLIM")
#     dev.off()
# 
# # MAXENT
#   # run maxent model
#     f.MAXENTrun <- function(envidata, p, a, p.output.dir){
#       m <- maxent(x=envidata, p=p, a=a,
#                   removeDuplicates=T, overwrite=T, path=paste(p.output.dir, "/maxent", sep=""),
#                   #See Maxent.jar programme -> help for a list and summary of these commands
#                   args=c("removeduplicates=true", "linear=TRUE", "quadratic=TRUE", "product=TRUE", "threshold=TRUE", "hinge=TRUE",
#                          "responsecurves=true", "pictures=true", "plots=true", "writeplotdata=true",
#                          "askoverwrite=false", "skipifexists=false", "randomseed=false", "randomtestpoints=0", "jackknife=false",
#                          "writeclampgrid=false", "writemess=false", "warnings=true", "allowpartialdata=false",
#                          "betamultiplier=1.0"))    
#       return(m)}
#     mx <- f.MAXENTrun(data.envi[[1:2]], p, a, "C:/data/modelling/output/bin")
#     pred.mx <- predict(mx, x=data.envi[[1:2]])
#   
#   # plot graph
#     tiff(paste("C:/Users/LaurensG/Desktop/plots/ggplot_", species, "_3_maxent.tif", sep=""), width=16, height=16, unit="cm", res=600)
#     modelplot(pred.mx, species=species, x, x.norm, y, y.norm, model="MAXENT")
#     dev.off()
# 
# 
#   # Run a Boosted Regression Tree
#     f.BRTrun <- function(p, a, envidata, data.extent, learning.rate, bag.fraction){
#     r <- rbind(cbind(p, "occ"=1), cbind(a, "occ"=0))
#     brtdata <- cbind(r, extract(envidata, r[,1:2]))
#     nichemodel <- gbm.step(data=brtdata, gbm.x=4:length(colnames(brtdata)), gbm.y=3, family="bernoulli", tree.complexity=3, learning.rate=learning.rate, bag.fraction=bag.fraction)
#     prediction.values <- predict.gbm(nichemodel, as.data.frame(values(envidata)), n.trees=nichemodel$gbm.call$best.trees, type="response")
#     prediction <- data.extent
#     values(prediction) <- prediction.values
#     prediction <- mask(prediction, data.extent)
#     #return(list(nichemodel, prediction))}
#     return(prediction)}
#     pred.brt <- f.BRTrun(p, a, data.envi[[1:2]], data.extent, 0.01, 0.5)
#   # plot graph
#     tiff(paste("C:/Users/LaurensG/Desktop/plots/ggplot_", species, "_4_brt.tif", sep=""), width=16, height=16, unit="cm", res=600)
#     modelplot(pred.brt, species=species, x, x.norm, y, y.norm, model="BRT")
#     dev.off()

model="BIOCLIM"
bla <- mask(results[[2]][[1]], data.extent)
tiff(paste("C:/Users/LaurensG/Desktop/plots/ggplot_geodist_", species, "_2_bioclim.tif", sep=""), width=16, height=16, unit="cm", res=600)
modelplot(bla, species=species, x, x.norm, y, y.norm, model=model)
dev.off()

model="MAXENT"
bla <- mask(results[[2]][[2]], data.extent)
tiff(paste("C:/Users/LaurensG/Desktop/plots/ggplot_geodist_", species, "_3_maxent.tif", sep=""), width=16, height=16, unit="cm", res=600)
modelplot(bla, species=species, x, x.norm, y, y.norm, model=model)
dev.off()

model="BRT"
bla <- mask(results[[2]][[3]], data.extent)
tiff(paste("C:/Users/LaurensG/Desktop/plots/ggplot_geodist_", species, "_4_brt.tif", sep=""), width=16, height=16, unit="cm", res=600)
modelplot(bla, species=species, x, x.norm, y, y.norm, model=model)
dev.off()


}
###

###


####colour gradient
#scale_colour_gradientn(colours=rev(rainbow(4)))
