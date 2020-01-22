# Peak detect 2.0
#
# copyright (c) 2018-2021 - Brockmann group - HU Berlin Manuel Delpero 
# 
# first written december, 2019

setwd("C:/Users/Manuel/Desktop/AIL_S1xS2/RAWDATA")

lods <- read.table("lodmatrixADDDOM_nosum.txt", sep = "\t",  header = TRUE)
map <- read.csv("map.cleaned.txt", header=TRUE, sep="\t", check.names=FALSE)
map <- map[, c(1,2)]
lods <- t(lods)
phenotypes <- rownames(lods)

chromosomes <- c(1:19, "X", "Y")

annotation <- c()
for(chr in chromosomes){
  annotation <- rbind(annotation, map[map[,"chr"] == chr,])
}


# Get the main QTLs for each phenotype
res <- c()
for (x in phenotypes){
  Lodscores <- lods[x,]
  ord <- sort(Lodscores, decreasing = TRUE)
  nqtl <- 0
  while ((any(ord > 4)) && (nqtl < 8)){
    info <- cbind(rownames(map[names(ord[1]),]) , x, map[names(ord[1]),], ord[1][[1]])   
    markers <- rownames(map[which(map[, "chr"] == info[, "chr"]),])    
    ord <- ord[-which(names(ord) %in% markers)]
    ord <- sort(ord, decreasing = TRUE)
    nqtl <- nqtl + 1
    res <- rbind(res, info)
  }
}
colnames(res) <- c("TopMarker", "Trait", "Chr", "PosTopMarker", "LodTopMarker")
rownames(res) <- NULL
res <- res[-which(is.na(res[, "Chr"] )),]

# Get regions with a lod drop of 1.5 from the top marker
rights <- c()
lefts <- c()
for (top in 1:nrow(res)){
  trait <- as.character(res[top, "Trait"])
  toplod <- as.character(res[top,1])
  topmark <- map[names(lods[trait,which((lods[trait, ] > lods[trait, toplod] - 1.5) & (lods[trait, ] < toplod))]), ]
  topmark <- topmark[which(topmark[, 1] == res[top, 3]),] 
  topmark <- topmark[,2]
  left <- (topmark[order(topmark, decreasing = FALSE)])[1]
  right <- tail((topmark[order(topmark, decreasing = FALSE)]), n = 1)
  lefts <- c(lefts, left)
  rights <- c(rights, right)
  }
 
 # to do (overlap the regions and reduce the dimensions)
regions <- cbind(res, "LeftPos" = lefts, "RightPos" = rights)  
QTLregions <- cbind(as.character(regions[, "Trait"]), as.numeric(as.character(regions[, "Chr"])), regions[, "LeftPos"], regions[, "PosTopMarker"], regions[, "RightPos"], as.numeric(as.character(regions[,"LodTopMarker"])),  as.character(regions[, "TopMarker"]))
colnames(QTLregions) <- c("Phenotype",	"Chr",	"StartPos",	"TopPos",	"StopPos",	"LOD", "TopMarker")
write.table(QTLregions, file = "QTLregions2212020.txt", quote = FALSE, sep = "\t", row.names = FALSE)