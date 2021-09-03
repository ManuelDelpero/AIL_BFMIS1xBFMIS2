# 
#
# copyright (c) 2018-2021 - Brockmann group - HU Berlin Manuel Delpero
# 
# 

setwd("C:/Users/Manuel/Desktop/AIL_S1xS2/RAWDATA")

data <- read.csv("qPCR_clean.txt", header = TRUE, check.names = FALSE, sep="\t")

data[,"dCT"] <- data[,"meanG"] - data[,"meanH"]
data[,"Transf"] <- 2^(-data[,"dCT"])
genes <- unique(data[,"Target Name"])

results <- matrix(NA, length(as.character(unique(data[,"Target Name"]))), 4,)
colnames(results) <- c("meanS1", "meanS2", "2^-(ddCT)", "p-value")
rownames(results) <- as.character(unique(data[,"Target Name"]))

for (gene in as.character(genes)){
  currentG <- data[grep(gene, data[,"Target Name"]),]
  p <- t.test(currentG[grep("S1", currentG[,"Sample Name"]), "Transf"], currentG[grep("S2", currentG[,"Sample Name"]), "Transf"])[[3]]
  meanS1 <- mean(currentG[grep("S1", currentG[,"Sample Name"]), "dCT"], na.rm = TRUE)
  meanS2 <- mean(currentG[grep("S2", currentG[,"Sample Name"]), "dCT"], na.rm = TRUE)
  fold <- 2^-(meanS1-meanS2)
  results[gene,"meanS1"] <- meanS1
  results[gene,"meanS2"] <- meanS2
  results[gene,"p-value"] <- p
  results[gene,"2^-(ddCT)"] <- fold
}
  
  
  
