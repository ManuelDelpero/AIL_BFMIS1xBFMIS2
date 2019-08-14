# Gnese expression analysis using microarrays from Atlas

library(dplyr)
library(preprocessCore)
setwd("C:/Users/Manuel/Desktop/AIL_S1xS2/Microarray Data/DN-2019_8745-Data/Analyse")


arrays <- read.csv("8745_Ergebnisse.SST-RMA-GENE-FULL - Group 2.TXT", header=TRUE, sep="\t")

# QC 
col <- names(arrays)
arrays.signals <- select(arrays, col[grepl("Signal", names(arrays))]) 
boxplot(arrays.signals)
arrays.signalsNorm <- normalize.quantiles(as.matrix(arrays.signals))
boxplot(arrays.signalsNorm)
colnames(arrays.signalsNorm) <- colnames(arrays.signals)

colarrays <- gsub("_20190802_.Clariom_S_Mouse..sst.rma.gene.full.Signal", "", colnames(arrays.signalsNorm))
colarrays <- gsub("_20190730_.Clariom_S_Mouse..sst.rma.gene.full.Signal", "", colarrays)
colarrays <- gsub("X8745_", "", colarrays)
colnames(arrays.signalsNorm) <- colarrays

par(mfrow=c(4,4))
for(x in 1:length(colnames(arrays.signalsNorm))){
  correlations <- as.numeric(cor(arrays.signalsNorm[,x], arrays.signalsNorm[,-x]))
  plot(correlations, ylab="Correlation", xlab="",  xaxt="n")
  
}

heatmap(cor(arrays.signalsNorm))

# Gene expression analysis
col <- colarrays
GF <- arrays.signalsNorm[,c(16:30)]
GFS1 <- GF[, c(1:7)]
GFS2 <- GF[, c(7:15)]
L <- arrays.signalsNorm[,c(1:15)]
LS1 <- L[, c(1:7)]
LS2 <- L[, c(7:15)]

# Differential gene expression analysis in gonadal fat samples
resultsGF <- NULL
for(x in 1:nrow(arrays.signalsNorm)){
  gf1 <- as.numeric(GFS1[x,])
  gf2 <- as.numeric(GFS2[x,])
  values <- c(mean(gf1), mean(gf2), t.test(gf1, gf2)$p.value)
  resultsGF <- rbind(resultsGF, values)
}

rownames(resultsGF) <- arrays[,1]
colnames(resultsGF) <- c("meanGFS1","meanGFS2","Ttest")

length(which(resultsGF[,"Ttest"] < (0.05 / nrow(resultsGF))))
length(which(p.adjust(resultsGF[,"Ttest"], "BH") < 0.05))

# Differential gene expression analysis in liver samples
resultsL <- NULL
for(x in 1:nrow(arrays.signalsNorm)){
  lS1 <- as.numeric(LS1[x,])
  lS2 <- as.numeric(LS2[x,])
  values <- c(mean(lS1), mean(lS2), t.test(lS1, lS2)$p.value)
  resultsL <- rbind(resultsL, values)
}

rownames(resultsL) <- arrays[,1]
colnames(resultsL) <- c("meanLS1","meanLS2","Ttest")

length(which(resultsL[,"Ttest"] < (0.05 / nrow(resultsL))))


