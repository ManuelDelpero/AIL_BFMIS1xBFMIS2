# AIL_S1xS2 CTL mapping plots
#
# copyright (c) - Manuel Delpero
# first written july, 2020
# 

setwd("C:/Users/Manuel/Desktop/AIL_S1xS2/RAWDATA")

pheno <- read.table("allPhenotypes_final.txt", sep = "\t", row.names=1)
genotypes <- read.csv("genotypesComplete.txt", header = TRUE, check.names = FALSE, sep="\t", colClasses="character")
lodmatrixC <- read.table("lodscoresCTL.txt", sep = "\t", check.names = FALSE, header = TRUE)
lodmatrixQ <- read.table("lodmatrixADDComplete.txt", sep = "\t", check.names = FALSE, header = TRUE)

markerannot <- read.csv("map.cleaned.txt", header=TRUE, sep="\t", check.names=FALSE)
markerannot <- markerannot[, c(1,2)]
markerannot$Index <- seq.int(nrow(markerannot))
colnames(markerannot) <- c("Chromosome", "Position", "Index")
markerannot <- markerannot[-which(is.na(markerannot[,2])),]

chromosomes <- c(1:19, "X", "Y")

annotation <- c()
for(chr in chromosomes){
  annotation <- rbind(annotation, markerannot[markerannot[,"Chromosome"] == chr,])
}

lodmatrixC <- lodmatrixC[rownames(annotation),]
lodmatrixC <- data.frame(lodmatrixC)
rownames(lodmatrixC) <- rownames(annotation)
lodannotmatrixC <- cbind(annotation, lodmatrixC)
lodmatrixQ <- lodmatrixQ[rownames(annotation),]
lodmatrixQ<- data.frame(lodmatrixQ)
rownames(lodmatrixQ) <- rownames(annotation)
lodannotmatrixQ <- cbind(annotation, lodmatrixQ)


# Order columns by gon weight
index <- order(pheno[,"Gon"], decreasing = FALSE)
sortWeight <- pheno[index,]
sortWeight[,c("Gon", "Leber")] <- sortWeight[,c("Gon", "Leber")] * sortWeight[,"Gewicht"]
# Calculate correlation
data <- pheno[, c("Gon", "Leber", "Gluc172")]
cor(data, use = "pairwise.complete.obs")

# Plot to represent the correlation between the tissues weight
par(cex.lab=1.2, cex.main = 1.3, cex.axis = 1)
mat <- matrix(c(1,2,3), 1, ,byrow = TRUE)
layout(mat, widths = rep.int(3, ncol(mat)))

plot(main="Tissues weight sorted by gonadal fat weight", c(1, nrow(sortWeight)), c(0, max(sortWeight[,"Gon"], na.rm=TRUE)), t = "n", xlab="Individuals", ylab="Weight [g]", ylim = c(0,7), las = 2, xaxt = "n")
  axis(1, at = c(0, 100, 200, 300, 400, 500), c("0", "100", "200", "300", "400", "500"))
  points(sortWeight[,"Gon"], col = "blue" , lwd=0.8 , pch=20 , type="p")
  points(sortWeight[,"Leber"], col = "darkgreen" , lwd=0.8 , pch=20 , type="p")
  points(sortWeight[,"SCF"], col = "orange" , lwd=0.8 , pch=20 , type="p")
  legend("topleft",
  legend = c("Liver", "Gonadal fat", "SCF"),
   col = c("darkgreen", "blue", "orange"),
   pch = c(20,20,20),
   bty = "n",
   pt.cex = 1.4,
   cex = 1,)

GonLiver <- pheno[, c("Gon", "Leber")]
corGonLiver <- cor(GonLiver, use = "pairwise.complete.obs")
rownames(corGonLiver) <- c("Gonadal fat", "Liver")
colnames(corGonLiver) <- c("Gonadal fat", "Liver")
#write.table(corGonLiver, file = "CorrelationGonLiverWeight.txt", quote = FALSE, sep = "\t")

# Correlation plot gonadal fat, liver

chr15TopMarker <- t(genotypes[c("UNC25805470",1),])
#chr15TopMarker <- chr15TopMarker[-which(is.na(chr15TopMarker[,1])),]
sortWeight <- sortWeight[rownames(chr15TopMarker),]
S1 <- names(which(chr15TopMarker[,1] == "A"))
S2 <- names(which(chr15TopMarker[,1] == "B"))
HET <- names(which(chr15TopMarker[,1] == "H"))

plot(main="Correlation plot gonadal fat weight ~ liver weight", c(1,5), c(0,7), t = "n", xlab="Liver weigth [g]", ylab="Gonadal fat weight [g]", las = 2, xaxt = "n")
  mycols = c()
  for (x in 1:nrow(sortWeight)){
    if (rownames(sortWeight[x,]) %in% S1){
	  mycol <- "orange"
	}else if (rownames(sortWeight[x,]) %in% S2){
	  mycol <- "blue"
	}else if (rownames(sortWeight[x,]) %in% HET){
	  mycol <- "black"
	}
  mycols <- c(mycols, mycol)
  }
  points(sortWeight[,"Leber"],sortWeight[,"Gon"], lwd=0.8 , pch=16 , type="p", col = mycols)
  axis(1, at = c(0, 1, 2, 3, 4, 5), c("0", "1", "2", "3", "4", "5"))
  abline(lm(sortWeight[,"Leber"]~sortWeight[,"Gon"]), col="red")
  text(3.5, 6, "r = -0.71, p-value < 2.2e-16")
  legend("topleft",
   legend = c("S1", "HET", "S2"),
   col = c("orange", "black", "blue"),
   pch = c(20,20,20),
   bty = "n",
   pt.cex = 1.4,
   cex = 1,)


# CTL mapping curve across chromosome 15
chr15C <- lodannotmatrixC[which(lodannotmatrixC[,"Chromosome"] == 15),]
chr15C <- chr15C[order(chr15C[,"Position"]),]
chr15Q <- lodannotmatrixQ[which(lodannotmatrixQ[,"Chromosome"] == 15),]
chr15Q <- chr15Q[order(chr15Q[,"Position"]),]

plot(main = "CTL profile liver ~ gonadal fat weight [Chr 15]", c(min(as.numeric(chr15Q[, "Position"])), max(as.numeric(chr15Q[, "Position"]))), c(-8,8), ylab = "-log10 [pvalue]", xlab = "Position [mb]", las = 2, t = "n", xaxt = "n")
  points(x = as.numeric(chr15C[,"Position"]), y = chr15C[,"lodmatrixC"] , type = "l", col="black", lwd = 2)
  points(x = as.numeric(chr15Q[,"Position"]), y = (chr15Q[,"Gon"]) , type = "l", col="blue", lwd = 1)
  points(x = as.numeric(chr15Q[,"Position"]), y = -(chr15Q[,"Leber"]) , type = "l", col="darkgreen", lwd = 1)
  abline(h=4.7, col="green")
  abline(h=4.2, col="orange")
  abline(h=-4.7, col="green")
  abline(h=-4.2, col="orange")
  abline(h=0, col="black")
  axis(1, at = c(0,25000000, 50000000, 75000000, 100000000), c("0", "25", "50", "75", "100"))
  legend("topleft",
   legend = c("QTL mapping liver", "QTL mapping gonadal fat", "CTL mapping"),
   col = c("brown", "blue", "black"),
   pch = c(20,20,20),
   bty = "n",
   pt.cex = 1.4,
   cex = 1,)