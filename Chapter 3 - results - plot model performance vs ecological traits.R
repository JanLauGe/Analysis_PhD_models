
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#
#exploratory scatterplots----
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#

specinfo <- read.csv("C:/Users/LaurensG/Desktop/specinfo.csv", na.string="#N/A")
results.aic2 <- read.csv("C:/data/modelling/output/results_aic_tuning2.csv")

#add some transformed variables
specinfo$polarity <- abs(specinfo$meanlat)
specinfo$log.nrecs.obis <- log(specinfo$nrecs.obis)
specinfo$log.nrecs.gbif <- log(specinfo$nrecs.gbif)
specinfo$log.nrecs <- log(specinfo$nrecs)

info <- merge
info <- merge(specinfo, results.aic2, by="TaxonKey", all.x=F, all.y=T)



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
      ggplot(data, aes_string(col1, col2)) + geom_point(size=0.5, colour=pointcolour) + labs(x=lab1, y=lab2) +
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
ylabs = c("BIOCLIM AUC", "MAXENT AUC", "BRT AUC")
ypars = c("baseline.fulld.auc.BIO", "baseline.fulld.auc.MAX", "baseline.fulld.auc.BRT")

plot11 <- singlescatterplot(specinfo, "fb_length", ypars[1], lab1=NULL, lab2=ylabs[1], pointcolour="#999999")
plot12 <- singlescatterplot(specinfo, "fb_weight", ypars[1], lab1=NULL, lab2=NULL, pointcolour="#999999")
plot13 <- singlescatterplot(specinfo, "fb_age", ypars[1], lab1=NULL, lab2=NULL, pointcolour="#999999")
plot14 <- singlescatterplot(specinfo, "shallow", ypars[1], lab1=NULL, lab2=NULL, pointcolour="#999999")
plot15 <- singlescatterplot(specinfo, "deep", ypars[1], lab1=NULL, lab2=NULL, pointcolour="#999999")
plot16 <- singlescatterplot(specinfo, "meanlat", ypars[1], lab1=NULL, lab2=NULL, pointcolour="#999999")
plot17 <- singlescatterplot(specinfo, "polarity", ypars[1], lab1=NULL, lab2=NULL, pointcolour="#999999")
#
plot21 <- singlescatterplot(specinfo, "fb_length", ypars[2], lab1=NULL, lab2=ylabs[2], pointcolour="#E69F00")
plot22 <- singlescatterplot(specinfo, "fb_weight", ypars[2], lab1=NULL, lab2=NULL, pointcolour="#E69F00")
plot23 <- singlescatterplot(specinfo, "fb_age", ypars[2], lab1=NULL, lab2=NULL, pointcolour="#E69F00")
plot24 <- singlescatterplot(specinfo, "shallow", ypars[2], lab1=NULL, lab2=NULL, pointcolour="#E69F00")
plot25 <- singlescatterplot(specinfo, "deep", ypars[2], lab1=NULL, lab2=NULL, pointcolour="#E69F00")
plot26 <- singlescatterplot(specinfo, "meanlat", ypars[2], lab1=NULL, lab2=NULL, pointcolour="#E69F00")
plot27 <- singlescatterplot(specinfo, "polarity", ypars[2], lab1=NULL, lab2=NULL, pointcolour="#E69F00")
#
plot31 <- singlescatterplot(specinfo, "fb_length", ypars[3], lab1="Length", lab2=ylabs[3], pointcolour="#56B4E9")
plot32 <- singlescatterplot(specinfo, "fb_weight", ypars[3], lab1="Weight", lab2=NULL, pointcolour="#56B4E9")
plot33 <- singlescatterplot(specinfo, "fb_age", ypars[3], lab1="Age", lab2=NULL, pointcolour="#56B4E9")
plot34 <- singlescatterplot(specinfo, "shallow", ypars[3], lab1="Upper depth limit", lab2=NULL, pointcolour="#56B4E9")
plot35 <- singlescatterplot(specinfo, "deep", ypars[3], lab1="Lower depth limit", lab2=NULL, pointcolour="#56B4E9")
plot36 <- singlescatterplot(specinfo, "meanlat", ypars[3], lab1="Mean latitude", lab2=NULL, pointcolour="#56B4E9")
plot37 <- singlescatterplot(specinfo, "polarity", ypars[3], lab1="Distance from equator", lab2=NULL, pointcolour="#56B4E9")

tiff("C:/Users/LaurensG/Desktop/plots/baseline PP ecoinfo.tiff", height=14, width=24, units="cm", res=600)
multiplot(plot11, plot21, plot31, plot12, plot22, plot32, plot13, plot23, plot33, 
          plot14, plot24, plot34, plot15, plot25, plot35, plot16, plot26, plot36,
          plot17, plot27, plot37, cols=7)
dev.off()



#~Boxplots of predictive performance with ecological trait categories----
library(cowplot)
library(reshape)

info[]

yvar = c("baseline.fulld.auc.BIO", "baseline.fulld.auc.MAX", "baseline.fulld.auc.BRT")
ylab = "BIOCLIM - AUC"

#Drop unwanted factor classes
data41 <- info[info$group1 != "other",c("group1", yvar)]; data41$group1 <- factor(data41$group1, levels=c("vertebrates","molluscs", "arthropods")); data41 <- melt(data41)
data42 <- info[,c("group2", yvar)]; data42$group2 <- factor(data42$group2, levels=c("ray-finned fish", "sharks and rays", "octopods and squids", "clams snails and slugs", "crabs and lobsters", "prawns and shrimps", "other")); data42 <- melt(data42)
data51 <- info[,c("dempel", yvar)]; data51$dempel <- factor(data51$dempel, levels=c("Demersal", "Pelagic")); data51 <- melt(data51)
data52 <- info[,c("fb_lifestyle", yvar)]; data52$fb_lifestyle <- factor(data52$fb_lifestyle, levels=c("catadromous", "amphidromous", "anadromous", "oceanodromous", "non-migratory")); data52 <- melt(data52)
data61 <- info[,c("habitat", yvar)]; data61$habitat <- factor(data61$habitat, levels=c("bathydemersal", "bathypelagic", "benthopelagic", "demersal", "pelagic-neritic", "pelagic-oceanic", "reef-associated")); data61 <- melt(data61)
data62 <- info[,c("fb_climate", yvar)]; data62$fb_climate <- factor(data62$fb_climate, levels=c("Polars", "Boreal", "Temperate", "Subtropical", "Tropical", "Deep-water")); data62 <- melt(data62)
data71 <- info[,c("Importance", yvar)]; data71$Importance <- factor(data71$Importance, levels=c("highly commercial", "commercial", "minor commercial", "subsistence fisheries", "of no interest")); data71 <- melt(data71)
data72 <- info[,c("MainCatchingMethod", yvar)]; data72$MainCatchingMethod <- factor(data72$MainCatchingMethod, levels=c("trawls", "seines", "hooks and lines", "gillnets", "traps", "various gears", "other")); data72 <- melt(data72)

plot41 <- ggplot(data41, aes(x=group1, y=value, fill=variable)) +  geom_boxplot(outlier.shape=1, outlier.size=1) + scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) + ggtitle("Taxonomic Group") + ylab(ylab) + theme(plot.margin=unit(c(0,0,0,0), "lines"), legend.position="none", axis.title.y=element_text(size=10), axis.text.y = element_text(size=8), plot.title=element_text(face="bold", size=12), axis.title.x = element_blank(), axis.text.x = element_text(angle=-30, hjust=0, vjust=1, size=8))
plot42 <- ggplot(data42, aes(x=group2, y=value, fill=variable)) +  geom_boxplot(outlier.shape=1, outlier.size=1) + scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) + ggtitle("Taxonomic Subgroup") + ylab("") + theme(plot.margin=unit(c(0,0,0,0), "lines"), legend.position="none", axis.title.y=element_text(size=10), axis.text.y = element_blank(), plot.title=element_text(face="bold", size=12), axis.title.x = element_blank(), axis.text.x = element_text(angle=-30, hjust=0, vjust=1, size=8))
plot51 <- ggplot(data51, aes(x=dempel, y=value, fill=variable)) +  geom_boxplot(outlier.shape=1, outlier.size=1) + scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) + ggtitle("Marine Zone") + ylab(ylab) + theme(plot.margin=unit(c(0,0,0,0), "lines"), legend.position="none", axis.title.y=element_text(size=10), axis.text.y = element_text(size=8), plot.title=element_text(face="bold", size=12), axis.title.x = element_blank(), axis.text.x = element_text(angle=-30, hjust=0, vjust=1, size=8))
plot52 <- ggplot(data52, aes(x=fb_lifestyle, y=value, fill=variable)) +  geom_boxplot(outlier.shape=1, outlier.size=1) + scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) + ggtitle("Lifestyle") + ylab("") + theme(plot.margin=unit(c(0,0,0,0), "lines"), legend.position="none", axis.title.y=element_text(size=10), axis.text.y = element_blank(), plot.title=element_text(face="bold", size=12), axis.title.x = element_blank(), axis.text.x = element_text(angle=-30, hjust=0, vjust=1, size=8))
plot61 <- ggplot(data61, aes(x=habitat, y=value, fill=variable)) +  geom_boxplot(outlier.shape=1, outlier.size=1) + scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) + ggtitle("Habitat") + ylab(ylab) + theme(plot.margin=unit(c(0,0,0,0), "lines"), legend.position="none", axis.title.y=element_text(size=10), axis.text.y = element_text(size=8), plot.title=element_text(face="bold", size=12), axis.title.x = element_blank(), axis.text.x = element_text(angle=-30, hjust=0, vjust=1, size=8))
plot62 <- ggplot(data62, aes(x=fb_climate, y=value, fill=variable)) +  geom_boxplot(outlier.shape=1, outlier.size=1) + scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) + ggtitle("Climate") + ylab("") + theme(plot.margin=unit(c(0,0,0,0), "lines"), legend.position="none", axis.title.y=element_text(size=10), axis.text.y = element_blank(), plot.title=element_text(face="bold", size=12), axis.title.x = element_blank(), axis.text.x = element_text(angle=-30, hjust=0, vjust=1, size=8))
plot71 <- ggplot(data71, aes(x=Importance, y=value, fill=variable)) +  geom_boxplot(outlier.shape=1, outlier.size=1) + scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) + ggtitle("Economic Importance") + ylab(ylab) + theme(plot.margin=unit(c(0,0,0,0), "lines"), legend.position="none", axis.title.y=element_text(size=10), axis.text.y = element_text(size=8), plot.title=element_text(face="bold", size=12), axis.title.x = element_blank(), axis.text.x = element_text(angle=-30, hjust=0, vjust=1, size=8))
plot72 <- ggplot(data72, aes(x=MainCatchingMethod, y=value, fill=variable)) +  geom_boxplot(outlier.shape=1, outlier.size=1) + scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) + ggtitle("Main Catching Method") + ylab("") + theme(plot.margin=unit(c(0,0,0,0), "lines"), legend.position="none", axis.title.y=element_text(size=10), axis.text.y = element_blank(), plot.title=element_text(face="bold", size=12), axis.title.x = element_blank(), axis.text.x = element_text(angle=-30, hjust=0, vjust=1, size=8))

tiff("C:/Users/LaurensG/Desktop/plots/baseline PP ecocat.tiff", height=23, width=16, units="cm", res=600)
plot_grid(plot41, plot42, plot51, plot52, plot61, plot62, plot71, plot72, labels=c("a)", "b)", "c)", "d)", "e)", "f)", "g)", "h)"), ncol = 2, nrow = 4, rel_widths = 1, rel_heights = 1,  align = "h")
dev.off()

