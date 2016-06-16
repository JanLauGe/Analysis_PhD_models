library(GGally)

specinfo <- read.csv("C:/Users/LaurensG/Desktop/specinfo.csv")
specinfo <- apply(specinfo, MARGIN=2, function(x) {as.numeric(x)})
ecoinfo <- specinfo[,c("fb_length", "fb_weight", "fb_age", "usual.shallow", "usual.deep", "fb_latlon", "meanlat", "rsIUCN", "rsSAUP", "nrecs.obis", "nrecs.gbif", "nrecs")]
ecoinfo1 <- specinfo[,c("fb_length", "fb_weight", "fb_age", "nrecs.obis", "nrecs.gbif", "nrecs")]
ecoinfo2 <- specinfo[,c("usual.shallow", "usual.deep", "fb_latlon", "meanlat", "rsIUCN", "rsSAUP", "nrecs.obis", "nrecs.gbif", "nrecs")]

tiff("C:/Users/LaurensG/Desktop/ecoinfo plot 1.tiff", units="cm", height=32, width=32, res=300)
ggpairs(ecoinfo1)
dev.off()
tiff("C:/Users/LaurensG/Desktop/ecoinfo plot 2.tiff", units="cm", height=32, width=32, res=300)
ggpairs(ecoinfo2)
dev.off()
