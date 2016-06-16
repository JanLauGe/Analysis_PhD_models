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

bla <- rbind(data.frame("model"=as.factor("BIOCLIM"), "run"=as.factor("baseline"), "metric"=as.factor("auc"), "value"=all[,"BIO.auc"]),
             data.frame("model"=as.factor("MAXENT"), "run"=as.factor("baseline"), "metric"=as.factor("auc"), "value"=all[,"MAX.auc"]),
             data.frame("model"=as.factor("BRT"), "run"=as.factor("baseline"), "metric"=as.factor("auc"), "value"=all[,"BRT.auc"]),
             data.frame("model"=as.factor("BIOCLIM"), "run"=as.factor("baseline"), "metric"=as.factor("cor"), "value"=all[,"BIO.cor"]),
             data.frame("model"=as.factor("MAXENT"), "run"=as.factor("baseline"), "metric"=as.factor("cor"), "value"=all[,"MAX.cor"]),
             data.frame("model"=as.factor("BRT"), "run"=as.factor("baseline"), "metric"=as.factor("cor"), "value"=all[,"BRT.cor"]))
             

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
#library(easyGgplot2)

data.plot1 <- bla[bla$metric=="auc",]
data.plot2 <- bla[bla$metric=="cor",]

ggplot(data.plot1, aes(x=model, y=value, fill=model, facet=run)) + geom_boxplot()

tiff("C:/Users/LaurensG/Desktop/plots/barplots_aic.tiff", width=16, height=12, unit="cm", res=600)
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

ggplot2.boxplot(data=data.plot, xName="run", yName="value", groupName="model",
                faceting=TRUE, facetingVarNames="model", facetingDirection="horizontal", facetingFont=c(12, "bold", "black"),
                width=1, notch=T, outlier.shape=21, groupColors=c('#999999','#E69F00','#56B4E9'), showLegend=FALSE,
                #axisLine=c(0.5, "solid", "black"),  backgroundColor="white", removePanelGrid=TRUE,removePanelBorder=TRUE,
                xtitle="Modelrun", ytitle="COR") + 
                theme(axis.text.x=element_text(angle=-45, hjust=0, color="gray40"), axis.text.y=element_text(color="gray40"), 
                      axis.title.x=element_text(vjust=-0.5, hjust=0.5), axis.title.y=element_text(vjust=1.5),
                      plot.margin=unit(c(2,10,2,2),"mm")) + expand_limits(y=c(0.6,1))
dev.off()



results <- merge(base.f, species.list, by="TaxonKey")
plot(results[,c("BIO.auc", "MAX.auc", "BRT.auc", "BIO.cor", "MAX.cor", "BRT.cor")])

library("GGally")
ggpairs(data=results, columns=2:7, axisLabels="none")
          
pairsplotting <- function (data, columns = 1:ncol(data), title = "", upper = list(), 
                           lower = list(), diag = list(), params = NULL, ..., axisLabels = "show", 
                           columnLabels = colnames(data[, columns]), legends = FALSE, 
                           verbose = FALSE) 
{
  printInfo <- FALSE
  verbose = verbose || printInfo
  axisLabelChoices <- c("show", "internal", "none")
  axisLabelChoice <- pmatch(axisLabels, axisLabelChoices)
  if (is.na(axisLabelChoice)) {
    warning("axisLabels not in c('show', 'internal', 'none').  Reverting to 'show'")
    axisLabelChoice <- 1
  }
  axisLabels <- axisLabelChoices[axisLabelChoice]
  if (any(columns > ncol(data))) {
    stop(paste("Make sure your 'columns' values are less than ", 
               ncol(data), ".\n\tcolumns = c(", paste(columns, collapse = ", "), 
               ")", sep = ""))
  }
  if (any(columns < 1)) {
    stop(paste("Make sure your 'columns' values are positive.", 
               "\n\tcolumns = c(", paste(columns, collapse = ", "), 
               ")", sep = ""))
  }
  if (any((columns%%1) != 0)) {
    stop(paste("Make sure your 'columns' values are integers.", 
               "\n\tcolumns = c(", paste(columns, collapse = ", "), 
               ")", sep = ""))
  }
  if (length(columnLabels) != length(columns)) {
    stop("The length of the 'columnLabels' does not match the length of the 'columns' being used.")
  }
  if (!is.list(upper) && upper == "blank") {
    upper <- list()
    upper$continuous = "blank"
    upper$combo = "blank"
    upper$discrete = "blank"
  }
  if (!is.list(lower) && lower == "blank") {
    lower <- list()
    lower$continuous = "blank"
    lower$combo = "blank"
    lower$discrete = "blank"
  }
  if (!is.list(diag) && diag == "blank") {
    diag <- list()
    diag$continuous = "blank"
    diag$discrete = "blank"
  }
  if (!is.list(upper)) 
    stop("upper is not a list")
  if (is.null(upper$continuous)) {
    upper$continuous <- "cor"
  }
  if (is.null(upper$combo)) {
    upper$combo <- "box"
  }
  if (is.null(upper$discrete)) {
    upper$discrete <- "facetbar"
  }
  if (!is.list(lower)) 
    stop("lower is not a list")
  if (is.null(lower$continuous)) {
    lower$continuous <- "points"
  }
  if (is.null(lower$combo)) {
    lower$combo <- "facethist"
  }
  if (is.null(lower$discrete)) {
    lower$discrete <- "facetbar"
  }
  if (is.null(diag$continuous)) {
    diag$continuous <- "density"
  }
  if (is.null(diag$discrete)) {
    diag$discrete <- "bar"
  }
  data <- as.data.frame(data)
  for (i in 1:dim(data)[2]) {
    if (is.character(data[, i])) {
      data[, i] <- as.factor(data[, i])
    }
  }
  numCol <- length(columns)
  if (printInfo) 
    cat("data col: ", numCol, "\n")
  ggpairsPlots <- list()
  grid <- rev(expand.grid(y = 1:ncol(data[columns]), x = 1:ncol(data[columns])))
  all <- do.call("rbind", lapply(1:nrow(grid), function(i) {
    xcol <- grid[i, "x"]
    ycol <- grid[i, "y"]
    data.frame(xvar = names(data[columns])[ycol], yvar = names(data[columns])[xcol])
  }))
  if (printInfo) {
    cat("\n\n\nALL\n")
    print(all)
  }
  dataTypes <- plot_types(data[columns])
  if (printInfo) {
    cat("\n\n\nDATA TYPES\n")
    print(dataTypes)
  }
  if (identical(axisLabels, "internal")) {
    dataTypes$Type <- as.character(dataTypes$Type)
    dataTypes$Type[dataTypes$posx == dataTypes$posy] <- "label"
    dataTypes$Type <- as.factor(dataTypes$Type)
  }
  for (i in 1:nrow(dataTypes)) {
    p <- "blank"
    type <- dataTypes[i, "Type"]
    posX <- as.numeric(as.character(dataTypes[i, "posx"]))
    posY <- as.numeric(as.character(dataTypes[i, "posy"]))
    xColName <- as.character(dataTypes[i, "xvar"])
    yColName <- as.character(dataTypes[i, "yvar"])
    up <- posX > posY
    if (printInfo) 
      cat("Pos #", i, "\t(", posX, ",", posY, ")\t type: ")
    section_aes <- section_params <- NULL
    if (type == "scatterplot") {
      if (printInfo) 
        cat("scatterplot\n")
      subType <- "points"
      if (up) {
        subType <- upper$continuous
        section_aes <- upper$aes_string
        section_params <- upper$params
      }
      else {
        subType <- lower$continuous
        section_aes <- lower$aes_string
        section_params <- lower$params
      }
      combo_aes <- addAndOverwriteAes(aes_string(x = xColName, 
                                                 y = yColName, ...), section_aes)
      if (subType == "density") {
        combo_aes <- addAndOverwriteAes(combo_aes, aes_string(group = combo_aes$colour))
        combo_aes
      }
      combo_params <- addAndOverwriteAes(params, section_params)
      p <- make_ggpair_text(subType, combo_aes, combo_params, 
                            printInfo)
    }
    else if (type == "box-hori" || type == "box-vert") {
      if (printInfo) 
        cat("box-hori-vert\n")
      subType <- "box"
      section_aes <- NULL
      if (up) {
        subType <- upper$combo
        section_aes <- upper$aes_string
        section_params <- upper$params
      }
      else {
        subType <- lower$combo
        section_aes <- lower$aes_string
        section_params <- lower$params
      }
      combo_aes <- addAndOverwriteAes(aes_string(x = xColName, 
                                                 y = yColName, ...), section_aes)
      if (subType != "dot") 
        combo_aes <- mapping_color_fill(combo_aes)
      combo_params <- addAndOverwriteAes(params, section_params)
      p <- make_ggpair_text(subType, combo_aes, combo_params, 
                            printInfo)
    }
    else if (type == "mosaic") {
      if (printInfo) 
        cat("mosaic\n")
      subType <- "facetbar"
      section_aes <- NULL
      if (up) {
        subType <- upper$discrete
        section_aes <- upper$aes_string
        section_params <- upper$params
      }
      else {
        subType <- lower$discrete
        section_aes <- lower$aes_string
        section_params <- lower$params
      }
      combo_aes <- addAndOverwriteAes(aes_string(x = xColName, 
                                                 y = yColName, ...), section_aes)
      combo_params <- addAndOverwriteAes(params, section_params)
      if (subType == "ratio") {
        p <- ggally_ratio(data[, c(yColName, xColName)])
      }
      else if (subType == "facetbar") {
        if (!is.null(combo_aes$colour)) {
          combo_aes <- addAndOverwriteAes(combo_aes, 
                                          aes_string(fill = combo_aes$colour))
        }
        p <- make_ggpair_text(subType, combo_aes, combo_params, 
                              printInfo)
      }
      else if (subType == "blank") {
        p <- "ggally_blank('blank')"
      }
      else {
        p <- ggally_text("Incorrect\nPlot", size = 6)
      }
    }
    else if (type == "stat_bin-num") {
      if (printInfo) 
        cat("stat_bin-num\n")
      subType <- diag$continuous
      combo_aes <- addAndOverwriteAes(aes_string(x = xColName, 
                                                 ...), diag$aes_string)
      if (subType != "density") 
        combo_aes <- mapping_color_fill(combo_aes)
      combo_params <- addAndOverwriteAes(params, diag$params)
      if (subType != "blank") {
        p <- make_ggpair_text(paste(subType, "Diag", 
                                    sep = "", collapse = ""), combo_aes, combo_params, 
                              printInfo)
      }
      else {
        p <- "blank"
      }
    }
    else if (type == "stat_bin-cat") {
      if (printInfo) 
        cat("stat_bin-cat\n")
      subType <- diag$discrete
      combo_aes <- addAndOverwriteAes(aes_string(x = xColName, 
                                                 ...), diag$aes_string)
      combo_aes <- mapping_color_fill(combo_aes)
      combo_params <- addAndOverwriteAes(params, diag$params)
      p <- make_ggpair_text(paste(subType, "Diag", sep = "", 
                                  collapse = ""), combo_aes, combo_params, printInfo)
    }
    else if (type == "label") {
      combo_aes <- addAndOverwriteAes(aes_string(x = xColName, 
                                                 ...), diag$aes_string)
      combo_params <- addAndOverwriteAes(params, diag$params)
      combo_params <- addAndOverwriteAes(combo_params, 
                                         c(label = columnLabels[posX]))
      p <- make_ggpair_text("diagAxis", combo_aes, combo_params, 
                            printInfo)
    }
    ggpairsPlots[[length(ggpairsPlots) + 1]] <- p
  }
  plotMatrix <- list(data = data, columns = columns, plots = ggpairsPlots, 
                     title = title, verbose = verbose, printInfo = printInfo, 
                     axisLabels = axisLabels, columnLabels = columnLabels, 
                     legends = legends, gg = NULL)
  attributes(plotMatrix)$class <- c("gg", "ggpairs")
  plotMatrix
}



assignInNamespace("ggally_cor", ggally_cor, "GGally")

tiff("C:/Users/LaurensG/Desktop/plots/pairs_plots_baseline.tiff", width=32, height=32, unit="cm", res=600)
pairsplot <- ggpairs(p, 
              lower=list(continuous="smooth", params=c(colour="blue")), 
              #diag=list(continuous="bar", params=c(colour="blue")), 
              upper=list(params=list(corSize=6)), axisLabels='show',
              theme(legend.position = "none", 
                    panel.grid.major = element_blank(), 
                    axis.ticks = element_blank(), 
                    panel.border = element_rect(linetype = "dashed", colour = "black", fill = NA)))

dev.off()

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

