#install.packages("devtools")
library(devtools)
#devtools::install_github("ropensci/rfishbase@rfishbase2.0")
devtools::install_github("ropensci/rfishbase")
library("rfishbase")

species.list <- read.csv("C:/data/modelling/data/species/speciesinfo.new.csv", header=T)

s=1
fao=NULL
notfound=NULL
#Find commercial species on fishbase
for(s in 1:length(species.list$TaxonName)){
  try(s.fao <- faoareas(species_list=as.vector(species.list$TaxonName[s])))
  if(exists("s.fao")){
    print(species.list$TaxonKey[s])
    s.fao <- cbind(TaxonKey=as.numeric(species.list$TaxonKey[s]), TaxonName=as.character(species.list$TaxonName[s]), s.fao)
    fao <- c(fao, list(s.fao))
  }else{
    print(c(species.list$TaxonKey[s], "not found"))
    notfound <- c(notfound, species.list$TaxonKey[s])
    s.fao <- cbind(TaxonKey=as.numeric(species.list$TaxonKey[s]), TaxonName=as.character(species.list$TaxonName[s]), Error="not found")
    fao <- c(fao, list(s.fao))
  }
  rm(s.fao)
}

?ecology
result <- ecology(species_list=as.vector(species.list$TaxonName[s]))
