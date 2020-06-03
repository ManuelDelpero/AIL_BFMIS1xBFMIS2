# Analysis for master student (Max)
#
# copyright (c) 2018-2021 - Brockmann group - HU Berlin Manuel Delpero 
# 
# first written december, 2019

setwd("C:/Users/Manuel/Desktop/AIL_S1xS2/RAWDATA")

phenotypes <- read.csv("allPhenotypes.txt", header = TRUE, check.names = FALSE, sep = "\t", colClasses = "character")
genotypes <- read.csv("genotypesComplete.txt", header = TRUE, check.names = FALSE, sep="\t", colClasses="character")
lodmatrixDOM <- read.table("lodmatrixDOMComplete.txt", header = TRUE, sep = "\t", check.names = FALSE)
lodmatrixADD <- read.table("lodmatrixADDComplete.txt", header = TRUE, sep = "\t", check.names = FALSE)
mprofiles <- read.table("lodmatrixADDDOMComplete.txt", header = TRUE, sep = "\t", check.names = FALSE)
markerannot <- read.csv("map.cleaned.txt", header=TRUE, sep="\t", check.names=FALSE)
KASPregions <- read.csv("KASP_regions.txt", sep = "\t", check.names = FALSE, header = TRUE)
markerannot <- markerannot[, c(1,2)]
markerannot$Index <- seq.int(nrow(markerannot))
colnames(markerannot) <- c("Chromosome", "Position", "Index")
markerannot <- markerannot[-which(is.na(markerannot[,2])),]
colnames(genotypes) <- gsub("AIL", "", colnames(genotypes))
phenotypes <- phenotypes[colnames(genotypes),]
chromosomes <- c(1:19, "X", "Y")

annotation <- c()
for(chr in chromosomes){
  annotation <- rbind(annotation, markerannot[markerannot[,"Chromosome"] == chr,])
}

lodmatrixDOM <- lodmatrixDOM[rownames(annotation),]
lodmatrixADD <- lodmatrixADD[rownames(annotation),]
lodmatrixADDDOM <- mprofiles[rownames(annotation),]


# Define LOD score values before and after KASP assay
KASPgenotypes <- c("UNC28010943", "UNC24184030", "UNCHS041907", "JAX00063853", "UNCHS043909", "UNC5812781", "UNC25805470", "UNCHS030444", "UNCHS041714", "UNCHS019508", "UNC27568354")
KASPregions <- KASPregions[which(KASPregions[, "Marker"] %in% KASPgenotypes),]
write.table(KASPregions, file = "KASPregionsLast.txt", sep = "\t", quote = FALSE)
