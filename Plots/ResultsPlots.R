# Plots for AIL BFMI S1xS2
#
# copyright (c) 2018-2021 - Brockmann group - HU Berlin Manuel Delpero 
# 
# first written december, 2019

setwd("C:/Users/Manuel/Desktop/AIL_S1xS2/RAWDATA")

lodmatrixDOM <- read.table("lodmatrixDOM.txt", header = TRUE, sep = "\t", check.names = FALSE)
lodmatrixADD <- read.table("lodmatrixADD.txt", header = TRUE, sep = "\t", check.names = FALSE)
lodmatrixADDDOM <- read.table("lodmatrixADDDOM_nosum.txt", header = TRUE, sep = "\t", check.names = FALSE)
map <- read.table("map.cleaned.txt", sep="\t")

chromosomes <- c(1:19, "X", "Y")

annotation <- c()
for(chr in chromosomes){
  annotation <- rbind(annotation, map[map[,"chr"] == chr,])
}


lodmatrixDOM <- lodmatrixDOM[rownames(annotation),]
lodmatrixADD <- lodmatrixADD[rownames(annotation),]
lodmatrixADDDOM <- lodmatrixADDDOM[rownames(annotation),]

# Preliminary visual check bodyweight 
bw <- lodmatrixADDDOM[,c(1:32)]
rotate <- function(x) t(apply(x, 2, rev))
image(rotate(bw))

# Gon weight
lodannotmatrix <- cbind(annotation[rownames(lodmatrixADDDOM), ], lodmatrixADDDOM)
chr17 <- lodannotmatrix[which(lodannotmatrix[,"chr"] == 17),]
chr17ord <- chr17[order(chr17[,2], decreasing = FALSE),]
plot(main = "QTL Gon weight [Chr 17]", c(min(as.numeric(chr17ord[, "bp_mm10"])), max(as.numeric(chr17ord[, "bp_mm10"]))), c(0,9), ylab = "-log10 [pvalue]", xlab = "Position [mb]", las = 2, t = "n", xaxt = "n")
points(x = as.numeric(chr17ord[,"bp_mm10"]), y = chr17ord[,"Gon"] , type = "l", col="dodgerblue", lwd = 1)

