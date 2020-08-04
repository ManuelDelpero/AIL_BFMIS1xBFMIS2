# AIL_S1xS2 Analysis on Phenotypes
#
# copyright (c) - Manuel Delpero
# first written july, 2020
# 

setwd("C:/Users/Manuel/Desktop/AIL_S1xS2/RAWDATA")

pheno <- read.table("allPhenotypes_final.txt", sep = "\t", row.names=1)
lodmatrix <- read.table("lodscoresCTL.txt", sep = "\t", check.names = FALSE, header = TRUE)
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

lodmatrix <- lodmatrix[rownames(annotation),]
lodmatrix <- data.frame(lodmatrix)
rownames(lodmatrix) <- rownames(annotation)
lodannotmatrix <- cbind(annotation, lodmatrix)


# Order columns by gon weight
index <- order(pheno[,"Gon"], decreasing = FALSE)
sortWeight <- pheno[index,]

# Calculate correlation
data <- pheno[, c("Gon", "Leber", "Gluc172")]
cor(data, use = "pairwise.complete.obs")
# Plot to represent the correlation between the tissues weight

plot(main="Tissues weight sorted by gondadal fat weight", c(1, nrow(sortWeight)), c(0, max(sortWeight[,"Gon"], na.rm=TRUE)), t = "n", xlab="Individuals", ylab="Weight (grams)", ylim = c(0,0.14), las = 2, xaxt = "n")
axis(1, at = c(0, 100, 200, 300, 400, 500), c("0", "100", "200", "300", "400", "500"))
lines(sortWeight[,"Gon"], col = "orange" , lwd=1 , pch=20 , type="l")
lines(sortWeight[,"Leber"], col = "blue" , lwd=1 , pch=20 , type="l")
lines(sortWeight[,"SCF"]/sortWeight[, "Gewicht"], col = "purple" , lwd=1 , pch=20 , type="l")
text(360, 0.10, "r = -0.7050693")
   legend("topleft",
   legend = c("Liver", "Gon", "SCF"),
   col = c("blue", "orange", "purple"),
   pch = c(20,20,20),
   bty = "n",
   pt.cex = 1.2,
   cex = 1,)

GonLiver <- pheno[, c("Gon", "Leber")]
corGonLiver <- cor(GonLiver, use = "pairwise.complete.obs")
rownames(corGonLiver) <- c("Gonadal fat", "Liver")
colnames(corGonLiver) <- c("Gonadal fat", "Liver")
#write.table(corGonLiver, file = "CorrelationGonLiverWeight.txt", quote = FALSE, sep = "\t")
   
# CTL mapping curve across chromosome 15

chr15 <- lodannotmatrix[which(lodannotmatrix[,"Chromosome"] == 15),]
chr15 <- chr15[order(chr15[,"Position"]),]
plot(main = "CTL profile liver ~ gonadal fat weight [Chr 15]", c(min(as.numeric(chr15[, "Position"])), max(as.numeric(chr15[, "Position"]))), c(0,5), ylab = "-log10 [pvalue]", xlab = "Position [mb]", las = 2, t = "n", xaxt = "n")
  points(x = as.numeric(chr15[,"Position"]), y = chr15[,"lodmatrix"] , type = "l", col="dodgerblue", lwd = 1)
  abline(h=4.7, col="green")
  abline(h=4.2, col="orange")
  axis(1, at = c(0,25000000, 50000000, 75000000, 100000000), c("0", "25", "50", "75", "100"))