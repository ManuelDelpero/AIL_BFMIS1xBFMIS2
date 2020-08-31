# Basic statistics for manuscript tables
#
# copyright (c) 2018-2021 - Brockmann group - HU Berlin, Manuel Delpero
# last modified September, 2019
# first written August, 2019

setwd("C:/Users/Manuel/Desktop/AIL_S1xS2/RAWDATA")

genotypes <- read.csv("genotypesComplete.txt", header = TRUE, check.names = FALSE, sep="\t", colClasses="character")
phenotypes <- read.csv("PhenotypesComplete.txt", header = TRUE, check.names = FALSE, sep="\t", row.names=1)
markerannot <- read.csv("map.cleaned.txt", header=TRUE, sep="\t", check.names=FALSE)
lodmatrix <- read.csv("lodmatrixADDComplete.txt", header=TRUE, sep="\t", check.names=FALSE)
markerannot <- markerannot[rownames(genotypes),]
markerannot <- markerannot[order(markerannot[,"bp_mm10"]),]

chromosomes <- c(1:19, "X", "Y")

annotation <- c()
for(chr in chromosomes){
  annotation <- rbind(annotation, markerannot[markerannot[,"chr"] == chr,])
}

lodannotmatrix <- cbind(annotation[rownames(lodmatrix), ], lodmatrix)

# Calculate the original tissue weight
tissues <- colnames(phenotypes[, c(46:55)])
for (x in tissues){
  phenotypes[, x] <- phenotypes[, x] * phenotypes[, "Gewicht"]
}

# Make sure that the ordering between phenotypes and genotypes matches !!!!!
# Also sort the markers by their natural chromosome ordering
genotypes <- genotypes[, rownames(phenotypes)]

## Figure out means for each genotypes

# Gonadal fat chr3, 12, 15 and 17
stat <- function(lodannotmatrix, pheno = "Gon", chr = 3){  
  lodannotmatrix <- lodannotmatrix[, c("chr", "bp_mm10", pheno)]
  chr3 <- lodannotmatrix[which(lodannotmatrix[,"chr"] == chr ),]
  topmarkerID <- rownames(chr3[which.max(chr3[,pheno]),])
  topmarker <- t(genotypes[topmarkerID,])
  groupsSize <- apply(genotypes,1,  table)[[topmarkerID]]
  genopheno <- cbind(topmarker, phenotypes[,pheno])
  colnames(genopheno) <- c("Genotype", pheno)
  mean1 <- mean(as.numeric(genopheno[which(genopheno[, "Genotype"] == "A"), 2]), na.rm = TRUE)
  mean2 <- mean(as.numeric(genopheno[which(genopheno[, "Genotype"] == "B"), 2]), na.rm = TRUE)
  meanHET <- mean(as.numeric(genopheno[which(genopheno[, "Genotype"] == "H"), 2]), na.rm = TRUE)
  if (mean1 > mean2){
    meanS1 <- mean1
    meanS2 <- mean2
  }else{
    meanS1 <- mean2
    meanS2 <- mean1
  }
  return(c(pheno, meanS1, "S1", meanS2, "S2", meanHET, "HET"))
}
stat(lodannotmatrix, pheno = "Gon", chr = 3) 
stat(lodannotmatrix, pheno = "Gon", chr = 12) 
stat(lodannotmatrix, pheno = "Gon", chr = 15)
stat(lodannotmatrix, pheno = "Gon", chr = 17)

stat(lodannotmatrix, pheno = "Leber", chr = 17)

stat(lodannotmatrix, pheno = "Triglycerides", chr = 17)

stat(lodannotmatrix, pheno = "Gluc172", chr = 3)
stat(lodannotmatrix, pheno = "Gluc172", chr = 17)

stat(lodannotmatrix, pheno = "Triglycerides", chr = 17)
stat(lodannotmatrix, pheno = "Triglycerides", chr = 7)

stat(lodannotmatrix, pheno = "Triglycerides", chr = 17)