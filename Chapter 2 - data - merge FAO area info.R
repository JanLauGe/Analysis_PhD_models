load("C:/data/modelling/data/species/data.FAOareas")
species.list <- read.csv("C:/data/modelling/data/species/speciesinfo.csv", header=T)

s=1
allfao = NULL
for(s in 1:length(species.list$TaxonKey)){
  if(data.faoareas[[s]][3]=="not found"){
    allfao <- rbind(allfao, cbind(species.list$TaxonKey[s], as.character(species.list$TaxonName[s]), "NA"))
  }else{
    allfao <- rbind(allfao, cbind(species.list$TaxonKey[s], as.character(species.list$TaxonName[s]), as.character(data.faoareas[[s]][3])))
  }}

bla <- as.data.frame(allfao)
write.csv(bla, file="C:/data/modelling/data/species/faoareas fishbase.csv")
