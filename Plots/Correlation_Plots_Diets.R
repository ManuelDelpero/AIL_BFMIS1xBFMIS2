# AIL_S1xS2 Correlation and diets plots
#
# copyright (c) - Manuel Delpero
# first written august, 2021
# 

setwd("C:/Users/Manuel/Desktop/AIL_S1xS2/RAWDATA")

pheno <- read.csv("phenotypesCompleteAll.txt", header = TRUE, check.names = FALSE, sep="\t", row.names=1)
Parental_Glucose <- read.csv("Parental_Glucose.txt", header = TRUE, check.names = FALSE, sep="\t", row.names=1)
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
pheno[,c("Gon", "Leber")] <- pheno[,c("Gon", "Leber")] * pheno[,"Gewicht"]
par(cex.lab=1.6, cex.main = 1.5, cex.axis = 1.5)
mat <- matrix(c(1,2,3), 1, ,byrow = TRUE)
layout(mat, widths = rep.int(3, ncol(mat)))

data <- pheno[,c("Gluc172", "Gon", "Leber")]
cor(data, use = "pairwise.complete.obs")

# Liver and gonadal fat
plot(main="GonAT weight ~ liver weight", c(0,6), c(0,6), t = "n", xlab="Liver weight [g]", ylab="GonAT weight [g]", las = 2, xaxt = "n")
  axis(1, at = c(0, 1, 2, 3, 4, 5, 6, 7), c("0", "1", "2", "3", "4", "5", "6", "7"))
  points(pheno[,"Leber"],pheno[,"Gon"], lwd=0.8 , pch=16 , type="p")
  abline(lm(pheno[,"Gon"]~pheno[,"Leber"]), col="red")
  legend("topright", #bg="gray"
  bty = "n",
  cex = 1.5,
  legend = "r = -0.47, p < 2.2e-16",
  )
  
# Liver and Body weight
plot(main="Liver weight ~ body weight", c(0,6), c(0,70), t = "n", xlab="Liver weight [g]", ylab="Body weight [g]", las = 2, xaxt = "n")
  axis(1, at = c(0, 1, 2, 3, 4, 5, 6, 7), c("0", "1", "2", "3", "4", "5", "6", "7"))
  points(pheno[,"Leber"],pheno[,"D172"], lwd=0.8 , pch=16 , type="p")
  abline(lm(pheno[,"D172"]~pheno[,"Leber"]), col="red")
  legend("topright", #bg="gray"
  bty = "n",
  cex = 1.5,
  legend = "r = 0.65, p < 2.2e-16",
  )

# Gon and bodyweight
plot(main="GonAT weight ~ body weight", c(0,6), c(0,550), t = "n", xlab="GonAT weight [g]", ylab="Body weight [g]", las = 2, xaxt = "n")
  axis(1, at = c(0, 1, 2, 3, 4, 5, 6, 7), c("0", "1", "2", "3", "4", "5", "6", "7"))
  points(pheno[,"Gon"],pheno[,"D172"], lwd=0.8 , pch=16 , type="p")
  abline(lm(pheno[,"D172"]~pheno[,"Gon"]), col="red")
  legend("topright", #bg="gray"
  bty = "n",
  cex = 1.5,
  legend = "r = 0.02, p < 0.65",
  )

# Liver and glucose
plot(main="Liver weight ~ blood glucose", c(0,6), c(0,550), t = "n", xlab="Liver weight [g]", ylab="Blood glucose [mg/dl]", las = 2, xaxt = "n")
  axis(1, at = c(0, 1, 2, 3, 4, 5, 6, 7), c("0", "1", "2", "3", "4", "5", "6", "7"))
  points(pheno[,"Leber"],pheno[,"Gluc172"], lwd=0.8 , pch=16 , type="p")
  abline(lm(pheno[,"Gluc172"]~pheno[,"Leber"]), col="red")
  legend("topright", #bg="gray"
  bty = "n",
  cex = 1.5,
  legend = "r = 0.71, p < 2.2e-16",
  )
 
# Gon and glucose
plot(main="GonAT weight ~ blood glucose", c(0,6), c(0,550), t = "n", xlab="GonAT weight [g]", ylab="Blood glucose [mg/dl]", las = 2, xaxt = "n")
  axis(1, at = c(0, 1, 2, 3, 4, 5, 6, 7), c("0", "1", "2", "3", "4", "5", "6", "7"))
  points(pheno[,"Gon"],pheno[,"Gluc172"], lwd=0.8 , pch=16 , type="p")
  abline(lm(pheno[,"Gon"]~pheno[,"Gluc172"]), col="red")
  legend("topright", #bg="gray"
  bty = "n",
  cex = 1.5,
  legend = "r = -0.59, p < 2.2e-16",
  )
  
# Body weight and glucose
plot(main="Body weight ~ blood glucose", c(0,60), c(0,550), t = "n", xlab="Body weight [g]", ylab="Blood glucose [mg/dl]", las = 2, xaxt = "n")
  axis(1, at = seq(0,60,10), as.character(seq(0,60,10)))
  points(pheno[,"D172"],pheno[,"Gluc172"], lwd=0.8 , pch=16 , type="p")
  abline(lm(pheno[,"Gluc172"]~pheno[,"D172"]), col="red")
  legend("topright", #bg="gray"
  bty = "n",
  cex = 1.5,
  legend = "r = 0.31, p < 2.55e-10",
  )


# Plots for the diets
#AIL
GlucMean <- apply(GlucSeries,2,function(x){
  return(mean(x, na.rm = TRUE))
})

GlucSD <- apply(GlucSeries,2,function(x){
  return(sd(x, na.rm = TRUE))
})

par(cex.lab=1.3, cex.main = 1.3, cex.axis = 1.3)
glucose <- pheno[,75:87]
weeks <- as.numeric(gsub("Gluc", "", colnames(glucose)))/7
weeks <- round(weeks, digits = 1)
plot(main = "Blood glucose concentration over 25 weeks in the AIL", c(5,26), c(0,300), t = "n", xlab = "week", ,ylab="Glucose [mg/dl]", las = 2, xaxt = "n")
  axis(1, at = seq(5, 25, 1), as.character(seq(5, 25, 1)))
  points(weeks[1:5], GlucMean[1:5], lwd=3 , pch=16 , type="l", col = "blue")
  points(weeks[5:9], GlucMean[5:9], lwd=3 , pch=16 , type="l", col = "orange")
  points(weeks[9:13], GlucMean[9:13], lwd=3 , pch=16 , type="l", col = "red")
  segments(weeks, GlucMean - GlucSD ,weeks ,GlucMean + GlucSD, lty = 5)
  epsilon <- 0.02
  segments(weeks-epsilon,GlucMean-GlucSD,weeks+epsilon,GlucMean-GlucSD, lty = 5)
  segments(weeks-epsilon,GlucMean+GlucSD,weeks+epsilon,GlucMean+GlucSD, lty = 5)
  #abline(v=20, col="orange", lty = 5)
  #abline(v=23, col="orange", lty = 5)
  #abline(v=25, col="orange", lty = 5)
  #text(x = 12, y = 5, "Standard diet", cex = 1.2)
  #text(x = 21.5, y = 20, "High fat,low carb", cex = 1.2)
  #text(x = 24, y = 60, "High fat, high carb", cex = 1.2)
  #abline(v=20, col="black", lty = 1)
  legend("topleft",
  legend = c("Standard diet", "High fat, low carb diet", "High fat, high carb diet"),
   col = c("blue", "orange", "red"),
   lty = c(1,1,1),
   lwd = c(3,3,3),
   y.intersp=1.4,
   bty = "n",
   pt.cex = 1.4,
   cex = 1.4)

#S1 and S2
plot(main = "Blood glucose concentration over 25 weeks in the S1 line", c(5,26), c(0,300), t = "n", xlab = "week", ,ylab="Glucose [mg/dl]", las = 2, xaxt = "n")
  axis(1, at = seq(5, 25, 1), as.character(seq(5, 25, 1)))
  points(as.numeric(rownames(Parental_Glucose[1:16,]))/7, Parental_Glucose[1:16,1]*10, lwd=3 , pch=16, type="l", col = "blue")
  points(as.numeric(rownames(Parental_Glucose[16:31,]))/7, Parental_Glucose[16:31,1]*10, lwd=3 , pch=16, type="l", col = "red")
  segments(as.numeric(rownames(Parental_Glucose))/7, Parental_Glucose[,1]*10 - Parental_Glucose[,2]*10 ,as.numeric(rownames(Parental_Glucose))/7 ,Parental_Glucose[,1]*10 + Parental_Glucose[,2]*10)
  epsilon <- 0.02
  segments(as.numeric(rownames(Parental_Glucose))/7-epsilon,Parental_Glucose[,1]*10-Parental_Glucose[,2]*10,as.numeric(rownames(Parental_Glucose))/7+epsilon,Parental_Glucose[,1]*10-Parental_Glucose[,2]*10)
  #segments(x-epsilon,y+sd,x+epsilon,y+sd)
  segments(as.numeric(rownames(Parental_Glucose))/7-epsilon,Parental_Glucose[,1]*10+Parental_Glucose[,2]*10,as.numeric(rownames(Parental_Glucose))/7+epsilon,Parental_Glucose[,1]*10+Parental_Glucose[,2]*10)
  #text(x = 12, y = 10, "Standard diet", cex = 1.2)
  #text(x = 23, y = 10, "High fat diet", cex = 1.2)
  #abline(v=20, col="black", lty = 1)
  legend("topleft",
  legend = c("Standard diet", "High fat diet"),
   col = c("blue", "red"),
   lty = c(1,1),
   lwd = c(3,3,3),
   bty = "n",
   pt.cex = 1.4,
   cex = 1.4)
   
plot(main = "Blood glucose concentration over 25 weeks in the S2 line", c(5,26), c(0,300), t = "n", xlab = "week", ,ylab="Glucose [mg/dl]", las = 2, xaxt = "n")
  axis(1, at = seq(5, 25, 1), as.character(seq(5, 25, 1)))
  points(as.numeric(rownames(Parental_Glucose[1:16,]))/7, Parental_Glucose[1:16,3]*10, lwd=3 , pch=16, type="l", col = "blue")
  points(as.numeric(rownames(Parental_Glucose[16:31,]))/7, Parental_Glucose[16:31,3]*10, lwd=3 , pch=16, type="l", col = "red")
  segments(as.numeric(rownames(Parental_Glucose))/7, Parental_Glucose[,3]*10 - Parental_Glucose[,4]*10 ,as.numeric(rownames(Parental_Glucose))/7 ,Parental_Glucose[,3]*10 + Parental_Glucose[,4]*10)
  epsilon <- 0.02
  segments(as.numeric(rownames(Parental_Glucose))/7-epsilon,Parental_Glucose[,3]*10-Parental_Glucose[,4]*10,as.numeric(rownames(Parental_Glucose))/7+epsilon,Parental_Glucose[,3]*10-Parental_Glucose[,4]*10)
  #segments(x-epsilon,y+sd,x+epsilon,y+sd)
  segments(as.numeric(rownames(Parental_Glucose))/7-epsilon,Parental_Glucose[,3]*10+Parental_Glucose[,4]*10,as.numeric(rownames(Parental_Glucose))/7+epsilon,Parental_Glucose[,3]*10+Parental_Glucose[,4]*10)
  #text(x = 12, y = 10, "Standard diet", cex = 1.2)
  #text(x = 23, y = 10, "High fat diet", cex = 1.2)
  #abline(v=20, col="black", lty = 1)
  legend("topleft",
  legend = c("Standard diet", "High fat diet"),
   col = c("blue", "red"),
   lty = c(1,1),
   lwd = c(3,3,3),
   bty = "n",
   pt.cex = 1.4,
   cex = 1.4)