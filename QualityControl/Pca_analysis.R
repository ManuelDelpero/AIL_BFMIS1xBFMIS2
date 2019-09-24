# AIL_S1xS2 Analysis on Phenotypes
#
# copyright (c) - Manuel Delpero
# first written April, 2019
#
# Script for the analysis on the selected animals to genotype

setwd("C:/Users/Manuel/Desktop/AIL_S1xS2/RAWDATA")
mriLEAN <- read.csv("MRIlean.txt", sep = "\t", header=TRUE, check.names=FALSE)
mriFAT <- read.csv("MRIfat.txt", sep = "\t", header=TRUE, check.names=FALSE)
allPhenotypes <- read.csv("allPhenotypes.txt",sep = "\t", header=TRUE, check.names=FALSE, na.strings=c("NA", "", "-", "    ", "x"), colClasses="character")
extremes <- read.csv("groups.txt" ,sep = "\t", header=TRUE, check.names=FALSE)

# PCA with the extremes for Gon Fat and Triglycerides
GonTrig <- allPhenotypes[, c("Gon", "Triglycerides")]
GonTrig <- GonTrig[rownames(GonTrig) %in% extremes[,1],]
out <- which(is.na(GonTrig[,1]))
GonTrigEx <- GonTrig[-out,]
pcaGonTrigEx <- prcomp(GonTrigEx)
plot(pcaGonTrigeEx$x[,1], pcaGonTrigEx$x[,2])

subsetPCA <- allPhenotypes[which(apply(apply(allPhenotypes,1,is.na),2,sum) == 0), -c(1,2)]
subsetPCAnum <- apply(subsetPCA, 1, function(x){as.numeric(as.character(x))})
pcares <- prcomp(t(subsetPCAnum))
rownames(pcares$x) <- colnames(subsetPCAnum)
selected <- rownames(GonTrigEx)

pdf ("PCA-plot_allPhenotypes.pdf")
mcol <- as.numeric(rownames(pcares$x) %in% selected) + 1
plot(pcares$x[,1], pcares$x[,2], col = c("black", "orange")[mcol], type='p', pch=20)
dev.off()

# PCA with the extremes for all the samples
GonTrig <- allPhenotypes[, c("Gon", "Triglycerides")]
out1 <- which(is.na(GonTrig[,"Gon"]))
out2 <- which(is.na(GonTrig[,"Triglycerides"]))
GonTrig <- GonTrig[-c(out1,out2),]
pcaGonTrig <- prcomp(GonTrig)
plot(pcaGonTrig$x[,1], pcaGonTrig$x[,2])

pdf ("PCA-plot.pdf")
# Plot for PCA with all samples highlighting the selected animals
plot(main = "2D PCA-plot with all animals", x = c(min(pcaGonTrig$x[,1]), max(pcaGonTrig$x[,1])), y = c(min(pcaGonTrig$x[,2]), max(pcaGonTrig$x[,2])), type='p', pch=20, xlab='Comp.1(Triglycerides)', ylab='Comp.2(Gonadal Fat)')
  for (num in 1:length(pcaGonTrig$x[,1])){
    if (rownames(GonTrig[num,]) %in% rownames(GonTrigEx) == TRUE){
      points(pcaGonTrig$x[,1][[num]], pcaGonTrig$x[,2][[num]] , type='p', pch=20, col = "orange")
    }else{
      points(pcaGonTrig$x[,1][[num]], pcaGonTrig$x[,2][[num]] , type='p', pch=20)
    }
  }
dev.off()
  

# another function
allPhenotypes <- allPhenotypes[,-c(1,2)]

p <- princomp(na.omit(GonTrig))
loadings = p$loadings[]
p.variance.explained = p$sdev^2 / sum(p$sdev^2)
 
# plot percentage of variance explained for each principal component    
barplot(100*p.variance.explained, las=2, xlab='', ylab='% Variance Explained')

# Plot with bodyweight with all samples highlighting the selected animals





 
 

