# Get peaks using the DOM + ADD model
#
# copyright (c) 2018-2021 - Brockmann group - HU Berlin Manuel Delpero & Danny Arends
# 
# first written december, 2019

setwd("C:/Users/Manuel/Desktop/AIL_S1xS2/RAWDATA")

genotypes <- read.csv("genotypes.cleaned.txt", header = TRUE, check.names = FALSE, sep="\t", colClasses="character")
markerannot <- read.csv("map.cleaned.txt", header=TRUE, sep="\t", check.names=FALSE)
mprofiles <- read.table("lodmatrixADDDOM_nosum.txt", sep = "\t",  header = TRUE)
markerannot <- markerannot[, c(1,2)]
markerannot$Index <- seq.int(nrow(markerannot))
colnames(markerannot) <- c("Chromosome", "Position", "Index")
markerannot <- markerannot[rownames(genotypes),]
markerannot <- markerannot[-which(is.na(markerannot[,2])),]

results <- NULL
for(x in 1:ncol(mprofiles)){
  peaks <- peak.detect(mprofiles[,x], markerannot, loddrop = 1.5)
  if(!is.null(peaks)){
    for(p in 1:nrow(peaks)){
      leftPos <- markerannot[peaks[p,1], "Position"]
      topPos <- markerannot[peaks[p,2], "Position"]
      rightPos <- markerannot[peaks[p,3], "Position"]
      topChr <- as.character(markerannot[peaks[p,2], "Chromosome"])
      results <- rbind(results, c(colnames(mprofiles)[x], topChr, leftPos, topPos, rightPos, round(mprofiles[peaks[p,2],x],2), rownames(markerannot)[peaks[p,]]))
    }
  }
}

colnames(results) <- c("Phenotype", "Chr", "StartPos", "TopPos", "StopPos", "LOD", "flankLeft", "TopMarker", "FlankRight")
results <- data.frame(results)

write.table(results, "QTL_regions11320.txt", row.names=FALSE, quote=FALSE, sep='\t')