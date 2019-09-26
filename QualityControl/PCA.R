# AIL_S1xS2 Analysis on Phenotypes
#
# copyright (c) - Manuel Delpero
# first written April, 2019
#
# Script for the analysis on the selected animals to genotype with gigaMUGA array
# Have a look using PCA if the sample selected cover the variability of the popultaion

setwd("C:/Users/Manuel/Desktop/AIL_S1xS2/RAWDATA")
mriLEAN <- read.csv("MRIlean.txt", sep = "\t", header=TRUE, check.names=FALSE)
mriFAT <- read.csv("MRIfat.txt", sep = "\t", header=TRUE, check.names=FALSE)
allPhenotypes <- read.csv("allPhenotypes.txt",sep = "\t", header=TRUE, check.names=FALSE, na.strings=c("NA", "", "-", "    ", "x"), colClasses="character")
extremes <- read.csv("groups.txt" ,sep = "\t", header=TRUE, check.names=FALSE)

# PCA with all the phenotypes
subsetPCA <- allPhenotypes[which(apply(apply(allPhenotypes,1,is.na),2,sum) == 0), -c(1,2)]
subsetPCAnum <- apply(subsetPCA, 1, function(x){as.numeric(as.character(x))})
pcares <- prcomp(t(subsetPCAnum))
rownames(pcares$x) <- colnames(subsetPCAnum)
selected <- as.character(extremes[,1])

pdf ("PCA-plot_allPhenotypes.pdf")
mcol <- as.numeric(rownames(pcares$x) %in% selected) + 1
plot(main = "PCA-Plot all phenotypes", pcares$x[,1], pcares$x[,2], xlab = "PC1 (0.5285%)", ylab = "PC2 (0.2992%)", col = c("blue", "orange")[mcol], type='p', pch=20)
  legend("topright",
  legend = c("Non selected", "Selected"),
  col = c("blue", "orange"),
  pch = c(20,20,20),
  bty = "n",
  pt.cex = 2,
  cex = 1,
  text.col = "black",
  #horiz = F ,
  #inset = c(0.1, 0.1, 0.1)
)

dev.off()

# PCA with 2 variables (Gon Fat and liver Triglycerides)
GonTrig <- allPhenotypes[, c("Gon", "Triglycerides")]
out1 <- which(is.na(GonTrig[,"Gon"]))
out2 <- which(is.na(GonTrig[,"Triglycerides"]))
GonTrig <- GonTrig[-c(out1,out2),]
pcaGonTrig <- prcomp(GonTrig)
plot(pcaGonTrig$x[,1], pcaGonTrig$x[,2])

pdf ("PCA-plot.pdf")
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
