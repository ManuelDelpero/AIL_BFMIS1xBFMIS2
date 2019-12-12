# KASP for AIL BFMI S1xS2
#
# copyright (c) 2018-2021 - Brockmann group - HU Berlin Manuel Delpero & Danny Arends
# 
# first written december, 2019

setwd("C:/Users/Manuel/Desktop/AIL_S1xS2/RAWDATA")

map <- read.table("map.cleaned.txt", sep="\t")
lods <- read.table("lods.txt", sep = "\t",  header = TRUE)
chromosomes <- c(1:19, "X", "Y")
map <- map[,c(1,2)]

res <- c()

for (x in phenotypes){
  Lodscores <- lods[x,]
  ord <- sort(Lodscores, decreasing = TRUE)
  nqtl <- 0
  while ((any(ord > 4)) && (nqtl < 3)){
    info <- cbind(map[names(which.max(ord)),], rownames(Lodscores), max(Lodscores) )   
    markers <- rownames(map[which(map[, "chr"] == info[, "chr"]),])    
    Loscores <- Lodscores[-which(names(Lodscores) %in% markers)]
    nqtl <- nqtl + 1
    res <- rbind(res, info)
  }
}
colnames(res) <- c("Chr", "Pos", "Trait", "Lod")
