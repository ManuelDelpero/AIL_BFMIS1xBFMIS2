# AIL_S1xS2 CTL mapping plots
#
# copyright (c) - Manuel Delpero
# first written july, 2020
# 

setwd("C:/Users/Manuel/Desktop/AIL_S1xS2/RAWDATA")

pheno <- read.table("PhenotypesComplete.txt", sep = "\t", row.names=1)
pheno <- pheno[which(pheno[, "Sex"] == "m"),]
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
sortWeight <- sortWeight[order(sortWeight[, "Gon"]),]
# Calculate correlation
data <- pheno[, c("Gon", "Leber", "Gluc172", "D174", "Triglycerides")]
cor(data, use = "pairwise.complete.obs")

# Plot to represent the correlation between the tissues weight
par(cex.lab=1.5, cex.main = 1.5, cex.axis = 5)
mat <- matrix(c(1,2,3), 1, ,byrow = TRUE)
layout(mat, widths = rep.int(3, ncol(mat)))

par(mar = c(5, 4, 4, 4) + 0.3)  # Leave space for z axis
plot(main="Tissues weight sorted by gonadal fat weight", c(1, nrow(sortWeight)), c(0, max(sortWeight[,"Gon"], na.rm=TRUE)), t = "n", xlab="Individuals", ylab="", ylim = c(0,7), las = 2, xaxt = "n")
  axis(1, at = c(0, 100, 200, 300, 400, 500), c("0", "100", "200", "300", "400", "500"))
  points(sortWeight[,"Gon"], col = "blue" , lwd=0.8 , pch=20 , type="p")
  points(sortWeight[,"Leber"], col = "darkgreen" , lwd=0.8 , pch=20 , type="p")
  points(sortWeight[,"SCF"], col = "orange" , lwd=0.8 , pch=20 , type="p")
  legend("topleft",
  legend = c("Liver", "Gonadal fat", "SCF", "Blood glucose"),
   col = c("darkgreen", "blue", "orange", "red"),
   pch = c(20,20,20),
   bty = "n",
   pt.cex = 1.8,
   cex = 1.5,)
par(new = TRUE)
plot(1:nrow(sortWeight), sortWeight[,"Gluc172"], lwd=0.8 ,pch=20 ,type="p", col = "red", axes = FALSE, bty = "n", xlab = "", ylab = "", ylim = c(0, 800))
  axis(side=4, at = c(pretty(range(sortWeight[,"Gluc172"])), "700", "800"))
  mtext("Blood glucose [mg/dl]", side=4, line=3)
  mtext("Weight[Gr.]", side=2, line=3)
   
GonLiver <- pheno[, c("Gon", "Leber")]
corGonLiver <- cor(GonLiver, use = "pairwise.complete.obs")
rownames(corGonLiver) <- c("Gonadal fat", "Liver")
colnames(corGonLiver) <- c("Gonadal fat", "Liver")
#write.table(corGonLiver, file = "CorrelationGonLiverWeight.txt", quote = FALSE, sep = "\t")

# Correlation plot gonadal fat, liver and genotypes according to the CTL top marker

#chr15TopMarker <- t(genotypes[c("UNC25735466",1),])
#chr15TopMarker <- chr15TopMarker[-which(is.na(chr15TopMarker[,1])),]
#sortWeight <- sortWeight[rownames(chr15TopMarker),]
S1 <- names(which(chr15TopMarker[,1] == "B"))
S2 <- names(which(chr15TopMarker[,1] == "A"))
HET <- names(which(chr15TopMarker[,1] == "H"))
sortWeightS1 <- sortWeight[S1,]
sortWeightS2 <- sortWeight[S2,]
sortWeightHET <- sortWeight[HET,]
#sortWeight <- sortWeight[c(S1, S2, HET),]

plot(main="Correlation plot gonadal fat weight ~ liver weight", c(1,5), c(0,5), t = "n", xlab="Liver weigth [g]", ylab="Gonadal fat weight [g]", las = 2, xaxt = "n")
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
  abline(lm(sortWeightS1[,"Leber"]~sortWeightS1[,"Gon"]), col="orange")
  abline(lm(sortWeightS2[,"Leber"]~sortWeightS2[,"Gon"]), col="blue")
  #text(3.5, 6, "r = -0.71, p-value < 2.2e-16")
  legend("topleft",
   legend = c("BFMI-S1", "HET", "BFMI-S2"),
   col = c("orange", "black", "blue"),
   pch = c(20,20,20),
   bty = "n",
   pt.cex = 1.8,
   cex = 1.5,)
   
# Normal correlation plots
par(cex.lab=1.6, cex.main = 1.5, cex.axis = 1.5)
mat <- matrix(c(1,2,3), 1, ,byrow = TRUE)
layout(mat, widths = rep.int(3, ncol(mat)))

data <- pheno[,c("Gluc172", "Gon", "Leber")]
cor(data, use = "pairwise.complete.obs")

# Liver and gonadal fat
pheno[,c("Gon", "Leber")] <- pheno[,c("Gon", "Leber")] * pheno[,"Gewicht"]
plot(main="Correlation plot GonAT weight ~ liver weight", c(0,6), c(0,6), t = "n", xlab="Liver weight [g]", ylab="GonAT weight [g]", las = 2, xaxt = "n")
  axis(1, at = c(0, 1, 2, 3, 4, 5, 6, 7), c("0", "1", "2", "3", "4", "5", "6", "7"))
  points(pheno[,"Leber"],pheno[,"Gon"], lwd=0.8 , pch=16 , type="p")
  abline(lm(pheno[,"Gon"]~pheno[,"Leber"]), col="red")
  legend("topright", #bg="gray"
  bty = "n",
  cex = 1.5,
  legend = "r = -0.65, p < 2.2e-16",
  )

# Liver and glucose
plot(main="Correlation plot liver weight ~ blood concentration", c(0,6), c(0,550), t = "n", xlab="Liver weight [g]", las = 2, xaxt = "n")
  axis(1, at = c(0, 1, 2, 3, 4, 5, 6, 7), c("0", "1", "2", "3", "4", "5", "6", "7"))
  points(pheno[,"Leber"],pheno[,"Gluc172"], lwd=0.8 , pch=16 , type="p")
  abline(lm(pheno[,"Gluc172"]~pheno[,"Leber"]), col="red")
  legend("topright", #bg="gray"
  bty = "n",
  cex = 1.5,
  legend = "r = 0.74, p < 2.2e-16",
  )

# Gonadal fat and glucose
plot(main="Correlation plot GonAT weight ~ blood glucose concentration", c(0,6), c(0,550), t = "n", xlab="GonAT weight [g]", ylab="", las = 2, xaxt = "n")
  axis(1, at = c(0, 1, 2, 3, 4, 5, 6, 7), c("0", "1", "2", "3", "4", "5", "6", "7"))
  points(pheno[,"Gon"],pheno[,"Gluc172"], lwd=0.8 , pch=16 , type="p")
  abline(lm(pheno[,"Gluc172"]~pheno[,"Gon"]), col="red")
  legend("topright", #bg="gray"
  bty = "n",
  cex = 1.5,
  legend = "r = -0.62, p < 2.2e-16",
  )
  
# Split the two groups and see if there is a difference in correlation
FirstGroup <- sortWeight[1:400,]
SecondGroup <- sortWeight[401:nrow(sortWeight),]

FirstGroup <- FirstGroup[, c("Gon", "Leber", "Gluc172")]
cor(FirstGroup, use = "pairwise.complete.obs")
cor.test(FirstGroup[,"Leber"], FirstGroup[,"Gon"] , use = "pairwise.complete.obs")
cor.test(FirstGroup[,"Leber"], FirstGroup[,"Gluc172"] , use = "pairwise.complete.obs")
cor.test(FirstGroup[,"Gon"], FirstGroup[,"Gluc172"] , use = "pairwise.complete.obs")

SecondGroup <- SecondGroup[, c("Gon", "Leber", "Gluc172")]
cor(SecondGroup, use = "pairwise.complete.obs")
cor.test(SecondGroup[,"Leber"], SecondGroup[,"Gon"] , use = "pairwise.complete.obs")
cor.test(SecondGroup[,"Leber"], SecondGroup[,"Gluc172"] , use = "pairwise.complete.obs")
cor.test(SecondGroup[,"Gon"], SecondGroup[,"Gluc172"] , use = "pairwise.complete.obs")

# Liver and gonadal fat
plot(main="Correlation plot gonadal fat weight ~ liver weight", c(1,5), c(0,6), t = "n", xlab="Liver weight [g]", ylab="Gonadal fat weight [g]", las = 2, xaxt = "n")
  axis(1, at = c(0, 1, 2, 3, 4, 5, 6, 7), c("0", "1", "2", "3", "4", "5", "6", "7"))
  points(FirstGroup[,"Leber"],FirstGroup[,"Gon"], lwd=0.8 , pch=16 , type="p")
  abline(lm(FirstGroup[,"Leber"]~FirstGroup[,"Gon"]), col="red")
  text(3.5, 4, "r = -0.45, p-value < 2.2e-16")
  
# Liver and glucose
plot(main="Correlation plot liver weight ~ blood glucose", c(1,7), c(0,500), t = "n", xlab="Liver weight [g]", ylab="blood glucose [mg/dl]", las = 2, xaxt = "n")
  axis(1, at = c(0, 1, 2, 3, 4, 5, 6, 7), c("0", "1", "2", "3", "4", "5", "6", "7"))
  points(FirstGroup[,"Leber"],FirstGroup[,"Gluc172"], lwd=0.8 , pch=16 , type="p")
  abline(lm(FirstGroup[,"Gluc172"]~FirstGroup[,"Leber"]), col="red")
  text(5, 100, "r = 0.76, p-value < 2.2e-16")
  
# Gonadal fat and glucose
plot(main="Correlation plot gonadal fat weight ~ blood glucose", c(1,7), c(0,500), t = "n", xlab="Gonadal fat weight [g]", ylab="blood glucose [mg/dl]", las = 2, xaxt = "n")
  axis(1, at = c(0, 1, 2, 3, 4, 5, 6, 7), c("0", "1", "2", "3", "4", "5", "6", "7"))
  points(FirstGroup[,"Gon"],FirstGroup[,"Gluc172"], lwd=0.8 , pch=16 , type="p")
  abline(lm(FirstGroup[,"Gluc172"]~FirstGroup[,"Gon"]), col="red")
  text(5, 100, "r = -0.51, p-value < 2.2e-16")
  
# Liver and gonadal fat
plot(main="Correlation plot gonadal fat weight ~ liver weight", c(1,5), c(0,6), t = "n", xlab="Liver weight [g]", ylab="Gonadal fat weight [g]", las = 2, xaxt = "n")
  axis(1, at = c(0, 1, 2, 3, 4, 5, 6, 7), c("0", "1", "2", "3", "4", "5", "6", "7"))
  points(SecondGroup[,"Leber"],SecondGroup[,"Gon"], lwd=0.8 , pch=16 , type="p")
  abline(lm(SecondGroup[,"Leber"]~SecondGroup[,"Gon"]), col="red")
  text(2, 2, "r = -0.09, p-value < 0.3127")
  
# Liver and glucose
plot(main="Correlation plot liver weight ~ blood glucose", c(1,7), c(0,500), t = "n", xlab="Liver weight [g]", ylab="blood glucose [mg/dl]", las = 2, xaxt = "n")
  axis(1, at = c(0, 1, 2, 3, 4, 5, 6, 7), c("0", "1", "2", "3", "4", "5", "6", "7"))
  points(SecondGroup[,"Leber"],SecondGroup[,"Gluc172"], lwd=0.8 , pch=16 , type="p")
  abline(lm(SecondGroup[,"Gluc172"]~SecondGroup[,"Leber"]), col="red")
  text(5, 100, "r = 0.67, p-value < 2.2e-16")
  
# Gonadal fat and glucose
plot(main="Correlation plot gonadal fat weight ~ blood glucose", c(1,7), c(0,500), t = "n", xlab="Gonadal fat weight [g]", ylab="blood glucose [mg/dl]", las = 2, xaxt = "n")
  axis(1, at = c(0, 1, 2, 3, 4, 5, 6, 7), c("0", "1", "2", "3", "4", "5", "6", "7"))
  points(SecondGroup[,"Gon"],SecondGroup[,"Gluc172"], lwd=0.8 , pch=16 , type="p")
  abline(lm(SecondGroup[,"Gluc172"]~SecondGroup[,"Gon"]), col="red")
  text(5, 100, "r = -0.31, p-value < 0.0006452")



# CTL mapping curve across chromosome 15
chr15C <- lodannotmatrixC[which(lodannotmatrixC[,"Chromosome"] == 15),]
chr15C <- chr15C[order(chr15C[,"Position"]),]
chr15Q <- lodannotmatrixQ[which(lodannotmatrixQ[,"Chromosome"] == 15),]
chr15Q <- chr15Q[order(chr15Q[,"Position"]),]

plot(main = "CTL profile gonadal fat weight ~ liver weight  [Chr 15]", c(min(as.numeric(chr15Q[, "Position"])), max(as.numeric(chr15Q[, "Position"]))), c(-8,8), ylab = "-log10 [pvalue]", xlab = "Position [mb]", las = 2, t = "n", xaxt = "n")
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
   pt.cex = 1.8,
   cex = 1.5,)