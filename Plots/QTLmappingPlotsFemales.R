# Plots to represent the results of the QTL mapping in the females for AIL BFMI S1xS2
#
# copyright (c) 2018-2021 - Brockmann group - HU Berlin Manuel Delpero 
# 
# first written maj, 2020

setwd("C:/Users/Manuel/Desktop/AIL_S1xS2/RAWDATA")

lodmatrixADD <- read.csv("lodmatrixADDF.txt", header = TRUE, sep = "\t", check.names = FALSE)
lodmatrixDOM <- read.csv("lodmatrixDOM_nosumF.txt", header = TRUE, sep = "\t", check.names = FALSE)
lodmatrixADDOM <- read.csv("lodmatrixADDDOM_nosumF.txt", header = TRUE, sep = "\t", check.names = FALSE)
phenotypes <- read.csv("allPhenotypes.txt", header = TRUE, check.names = FALSE, sep = "\t", colClasses = "character")
genotypes <- read.csv("genotypesComplete.txt", header = TRUE, check.names = FALSE, sep="\t", colClasses="character")
phenotypes <- phenotypes[which(phenotypes[, "Sex"] == "f"),]
genotypes <- genotypes[, rownames(phenotypes)]

# In this case we just have LOD score values for 11 markers so we can just make so effect plots and no Manhattan plots
Sign<- -log10(0.005/11)
Hsign <- -log10(0.001/11)
apply(lodmatrixADD, 2, function(x){ x[which(x > Sign)]})
apply(lodmatrixDOM, 2, function(x){ x[which(x > Sign)]})
apply(lodmatrixADDOM, 2, function(x){ x[which(x > Sign)]})

# liver and d98 show QTLs but not the other phenotypes. One QTL on Chr 7 for liver weight and one QTL on chr 16 for bodyweight at day 98
# Effect plot top marker for liver weight
UNCHS019508 <- cbind(phenotypes[, "Leber"], t(genotypes["UNCHS019508",]))
boxplot(as.numeric(UNCHS019508[which(UNCHS019508[,2] == "A"),1]), as.numeric(UNCHS019508[which(UNCHS019508[,2] == "H"),1]), as.numeric(UNCHS019508[which(UNCHS019508[,2] == "B"),1]), col = c("lightskyblue1", "cyan3", "dodgerblue4"), main = "Liver weight [Chr 7]", ylab = "Weight [Gr.]", xlab = "Genotypes [UNCHS019508]" , las = 2, t = "n", xaxt = "n",  ylim = c(0, 4))
  axis(1, at = 1:3 , c("TT", "TC", "CC"))
  legend("topright", bg="gray",
  legend = c( "BFMI-S1", "HET", "BFMI-S2"),
    col = c("lightskyblue1", "cyan3", "dodgerblue4"),
    pch = 15,
    pt.cex = 1.7,
    pt.bg = "lightsteelblue1",
    cex = 1,
	bty = "n",
    text.col = "black")

# Effect plot top marker for D98
UNCHS041907 <- cbind(phenotypes[, "D98"], t(genotypes["UNCHS041907",]))
boxplot(as.numeric(UNCHS041907[which(UNCHS041907[,2] == "A"),1]), as.numeric(UNCHS041907[which(UNCHS041907[,2] == "H"),1]), as.numeric(UNCHS041907[which(UNCHS041907[,2] == "B"),1]), col = c("lightskyblue1", "cyan3", "dodgerblue4"), main = "Body weight day 98 [Chr 16]", ylab = "Weight [Gr.]", xlab = "Genotypes [UNCHS041907]" , las = 2, t = "n", xaxt = "n",  ylim = c(20, 45))
  axis(1, at = 1:3 , c("TT", "TC", "CC"))
  legend("topright", bg="gray",
  legend = c( "BFMI-S1", "HET", "BFMI-S2"),
    col = c("lightskyblue1", "cyan3", "dodgerblue4"),
    pch = 15,
    pt.cex = 1.7,
    pt.bg = "lightsteelblue1",
    cex = 1,
	bty = "n",
    text.col = "black")

