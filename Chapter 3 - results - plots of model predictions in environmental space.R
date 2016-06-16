
require(dismo)
require(ggplot2)
require(gtable)


#data.extent <- raster("C:/data/modelling/data/envi/data.extent.grd")
#data.envi <- stack(data.extent, data.envi[[1:7]])
data.envi <- stack("C:/data/modelling/data/envi/data.envi.grd")
load("C:/data/modelling/data/data.species.RData")

#####
species="Pleuronectes platessa"

p <- data.species[data.species$TaxonName==species,c(3,4)]
a <- data.species[data.species$set=="background",3:4]
x <- as.data.frame(na.omit(extract(data.envi[[1:2]],p)))
y <- as.data.frame(na.omit(extract(data.envi[[1:2]],a)))
bc <- bioclim(x)
pred <- predict(bc, x=data.envi[[1:2]])

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

x <- as.data.frame(normalise(x,y))
y <- as.data.frame(normalise(y,y))

#######
tiff("C:/Users/LaurensG/Desktop/ggplot grob.tif", width=16, height=16, unit="cm", res=300)
  scatter <- ggplot(y, aes(depthmean, sstanmean)) + 
    geom_point(colour="grey") + 
    geom_point(data=x, color="black") +
    theme(plot.margin=unit(c(0,0,0,0),"mm")) +
    xlim(-0.025,1) + ylim(-0.025,1) + xlab("Depth") + ylab("SST")
  histtop <- ggplot(y, aes(depthmean)) + 
    geom_histogram(colour="darkgrey", fill="grey", binwidth=0.025) + 
    geom_histogram(data=x, colour="darkgrey", fill="black", binwidth=0.025) +
    labs(x=NULL, y=NULL) + xlim(0,1) +
    theme(plot.margin=unit(c(0,0,0,0),"mm"), 
          legend.position="none", 
          axis.ticks=element_blank(), 
          axis.text=element_blank(), 
          axis.title=element_blank(), 
          plot.margin=unit(c(0,0,0,0),"mm"), 
          panel.grid.major=element_blank(), 
          panel.grid.minor=element_blank(),
          panel.background=element_blank())
  histright <- ggplot(y, aes(sstanmean)) + 
    geom_histogram(colour="darkgrey", fill="grey", binwidth=0.025) + 
    geom_histogram(data=x, colour="darkgrey", fill="black", binwidth=0.025) + 
    coord_flip() + labs(x=NULL, y=NULL) + xlim(0,1) +
    theme(plot.margin=unit(c(0,0,0,0),"mm"), 
          legend.position="none", 
          axis.ticks=element_blank(), 
          axis.text=element_blank(), 
          axis.title=element_blank(), 
          plot.margin=unit(c(0,0,0,0),"mm"), 
          panel.grid.major=element_blank(), 
          panel.grid.minor=element_blank(),
          panel.background=element_blank())

g <- ggplotGrob(scatter)
g <- gtable_add_cols(g, unit(3.77,"cm"))
g <- gtable_add_grob(g, ggplotGrob(histright), t=1, l=ncol(g), b=3, r=ncol(g))
g <- gtable_add_rows(g, unit(3.77,"cm"), 0)
g <- gtable_add_grob(g, ggplotGrob(histtop), t=1, l=4, b=1, r=4)

grid.newpage()
grid.draw(g)
grid.text(paste(unlist(strsplit(species, split=" ")),"\n",unlist(strsplit(species, split=" "))[2], 
                "\n\n", length(x[,1]), " records", sep=""),
          gp=gpar(fontsize=12), just="left", draw=T, x=0.79, y=0.87)

dev.off()


