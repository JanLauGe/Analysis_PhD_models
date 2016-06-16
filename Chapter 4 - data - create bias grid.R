library(sp)
library(raster)
library(rgdal)

setwd("C:/Users/LaurensG/Desktop/sampling bias/gridcreation")
data.extent <- raster("data.extent.grd")

load("results.1.RData")
load("results.2.RData")
load("results.3.RData")
load("results.4.RData")
load("results.5.RData")
load("results.6.RData")
load("results.7.RData")
load("results.8.RData")
load("results.9.RData")
load("results.10.RData")
load("results.11.RData")
load("results.12.RData")
load("results.13.RData")
load("results.14.RData")
load("results.15.RData")
load("results.16.RData")

gridvalues <- rbind(results.1, results.2, results.3, results.4, results.5, results.6, results.7, results.8, results.9, results.10, results.11, results.12, results.13, results.14, results.15, results.16)
#gridvalues <- c(results.1$SC, results.2$SC, results.3$SC, results.4$SC, results.5$SC, rep(10000, 4050) , results.7$SC, results.8$SC, results.9$SC, results.10$SC, results.11$SC, results.12$SC, results.13$SC, results.14$SC, results.15$SC, results.16$SC)

par(mfrow=c(3,2))

sb.ns <- data.extent
values(sb.ns) <- gridvalues$Species
plot(sb.ns, main="n species")
sb.nr <- data.extent
values(sb.nr) <- gridvalues$nrecords
plot(sb.nr, main="number of records")


sb.sc <- data.extent
values(sb.sc) <- gridvalues$SC
plot(sb.sc, main="iNEXT sample completeness")

sb.ch <- data.extent
values(sb.ch) <- gridvalues$chao
plot(sb.ch, main="chao")

sb.j1 <- data.extent
values(sb.j1) <- gridvalues$jack1
plot(sb.j1, main="jack1")

#sb.j2 <- data.extent
#values(sb.j2) <- gridvalues$jack2
#plot(sb.j2, main="jack2")

sb.bo <- data.extent
values(sb.bo) <- gridvalues$boot
plot(sb.bo, main="boot")

##################
sb.sc <- data.extent
values(sb.sc) <- gridvalues$SC
plot(sb.sc)
plot(mask(sb.sc, data.extent))

sb.all <- data.extent
values(sb.all) <- 4*gridvalues$Species / gridvalues$chao + gridvalues$jack1 + gridvalues$jack2 + gridvalues$boot
plot(sb.all)
plot(mask(sb.all, data.extent))

sb.cha <- data.extent
values(sb.cha) <- gridvalues$Species / gridvalues$chao
plot(sb.cha)
plot(mask(sb.cha, data.extent))
