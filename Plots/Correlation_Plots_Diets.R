# AIL_S1xS2 Correlation and diets plots
#
# copyright (c) - Manuel Delpero
# first written august, 2021
# 

setwd("C:/Users/Manuel/Desktop/AIL_S1xS2/RAWDATA")

pheno <- read.csv("phenotypesCompleteAll.txt", header = TRUE, check.names = FALSE, sep="\t", row.names=1)
GlucSeries <- read.csv("GlucSeries.txt", header = TRUE, check.names = FALSE, sep="\t", row.names = 1)
rownames(GlucSeries) <- gsub("V 888-", "", rownames(GlucSeries))
GlucSeries <- GlucSeries[rownames(pheno),]
pheno <- cbind(pheno, GlucSeries)
pheno <- pheno[which(pheno[, "Sex"] == "m"),]

# Order columns by gon weight
index <- order(pheno[,"Gon"], decreasing = FALSE)
sortWeight <- pheno[index,]
sortWeight[,c("Gon", "Leber")] <- sortWeight[,c("Gon", "Leber")] * sortWeight[,"Gewicht"]
sortWeight <- sortWeight[order(sortWeight[, "Gon"]),]
# Calculate correlation
#data <- phenotypes[, c( "auc","Gon", "Leber", "Gluc166", "D174", "Triglycerides")]
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
   
# correlation plots
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
  
# Liver and Body weight
plot(main="Correlation plot liver weight ~ Body weight", c(0,6), c(0,70), t = "n", xlab="Liver weight [g]", ylab="Body weight [g]", las = 2, xaxt = "n")
  axis(1, at = c(0, 1, 2, 3, 4, 5, 6, 7), c("0", "1", "2", "3", "4", "5", "6", "7"))
  points(pheno[,"Leber"],pheno[,"D172"], lwd=0.8 , pch=16 , type="p")
  abline(lm(pheno[,"D172"]~pheno[,"Leber"]), col="red")
  legend("topright", #bg="gray"
  bty = "n",
  cex = 1.5,
  legend = "r = 0.74, p < 2.2e-16",
  )

# Gon and bodyweight
plot(main="Correlation plot GonAT weight ~ Body weight", c(0,6), c(0,550), t = "n", xlab="GonAT weight [g]", las = 2, xaxt = "n")
  axis(1, at = c(0, 1, 2, 3, 4, 5, 6, 7), c("0", "1", "2", "3", "4", "5", "6", "7"))
  points(pheno[,"Gon"],pheno[,"D172"], lwd=0.8 , pch=16 , type="p")
  abline(lm(pheno[,"D172"]~pheno[,"Gon"]), col="red")
  legend("topright", #bg="gray"
  bty = "n",
  cex = 1.5,
  legend = "r = 0.74, p < 2.2e-16",
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
 
# Gon and glucose
plot(main="Correlation plot GonAT weight ~ blood concentration", c(0,6), c(0,550), t = "n", xlab="GonAT weight [g]", las = 2, xaxt = "n")
  axis(1, at = c(0, 1, 2, 3, 4, 5, 6, 7), c("0", "1", "2", "3", "4", "5", "6", "7"))
  points(pheno[,"Gon"],pheno[,"Gluc172"], lwd=0.8 , pch=16 , type="p")
  abline(lm(pheno[,"Gluc172"]~pheno[,"Gon"]), col="red")
  legend("topright", #bg="gray"
  bty = "n",
  cex = 1.5,
  legend = "r = 0.74, p < 2.2e-16",
  )
  
# Body weight and glucose
plot(main="Correlation plot Body weight ~ blood concentration", c(0,60), c(0,550), t = "n", xlab="Glucose [mg/dl]", las = 2, xaxt = "n")
  axis(1, at = c(0, 1, 2, 3, 4, 5, 6, 7), c("0", "1", "2", "3", "4", "5", "6", "7"))
  points(pheno[,"D172"],pheno[,"Gluc172"], lwd=0.8 , pch=16 , type="p")
  abline(lm(pheno[,"Gluc172"]~pheno[,"D172"]), col="red")
  legend("topright", #bg="gray"
  bty = "n",
  cex = 1.5,
  legend = "r = 0.74, p < 2.2e-16",
  )


# Plots for the diets
GlucMean <- apply(GlucSeries,2,function(x){
  return(mean(x, na.rm = TRUE))
})

plot(main = "", c(0,26), c(0,300), t = "n", xlab = "week", ,ylab="Glucose [mg/dl]", las = 2, xaxt = "n")
  axis(1, at = seq(0, 25, 1), as.character(seq(0, 25, 1)))
  points(c(6,10,14,18,20,21,22,22.5, 23, 23.5, 24, 24.5, 25), GlucMean, lwd=0.8 , pch=16 , type="l")

