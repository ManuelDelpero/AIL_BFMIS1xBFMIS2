# QTL mapping for the second paper of the BFMIAIL861-S1x861-S2
# Traits to use are: Bodyweight time series

setwd("C:/Users/Manuel/Desktop/AIL_S1xS2/RAWDATA")

genotypes <- read.csv("genotypes.cleaned.txt", header = TRUE, check.names = FALSE, sep="\t", colClasses="character")
phenotypes <- read.csv("PhenotypesComplete.txt", header = TRUE, check.names = FALSE, sep="\t", row.names=1)
markerannot <- read.csv("map.cleaned.txt", header=TRUE, sep="\t", check.names=FALSE)
markerannot <- markerannot[, c(1,2)]
markerannot <- markerannot[order(markerannot[,"bp_mm10"]),]
colnames(markerannot) <- c("Chromosome", "Position")
markerannot <- markerannot[-which(is.na(markerannot[,2])),]
chromosomes <- c(1:19, "X", "Y")

annotation <- c()
for(chr in chromosomes){
  annotation <- rbind(annotation, markerannot[markerannot[,"Chromosome"] == chr,])
}

colnames(annotation) <- c("Chromosome", "Position", "GenTrain Score", "SNP")
colnames(genotypes) <- gsub("AIL", "" , colnames(genotypes)) 
dim(genotypes)
dim(phenotypes)

# Get only the pheno we need
BW <- phenotypes[,grep("D", colnames(phenotypes))]
BW <- BW[,-1]

allPhenotypes <- BW
phenames <- colnames(allPhenotypes)

# Getting rid of the outliers
outliers <- apply(allPhenotypes[, phenames],2, function(x){
  up <- mean(x, na.rm = TRUE) + 3 * sd(x, na.rm = TRUE)
  low <- mean(x, na.rm=TRUE) - 3 * sd(x, na.rm=TRUE)
  x < low | x > up
 })
 
for(x in phenames){
  idx <- which(outliers[,x])
  if(length(idx) > 0) allPhenotypes[idx,x] <- NA
}
phenotypes <- cbind(phenotypes[,c(1, 58)], allPhenotypes)

## Adding genotypes by KASP assay to the genotypes matrix
# Adjust KASP dataset and combine it with the gigaMUGA genotypes
KASPgenotypes <- read.csv("KASPgenotypes.txt", sep = "\t", header = TRUE, check.names = FALSE, colClasses="character")
KASPgenotypes[] <- lapply(KASPgenotypes, function(x) gsub("Homozygous A/A|Homozygous T/T", "A", x))
KASPgenotypes[] <- lapply(KASPgenotypes, function(x) gsub("Homozygous G/G|Homozygous C/C", "B", x))
KASPgenotypes[] <- lapply(KASPgenotypes, function(x) gsub("Heterozygous A/C|Heterozygous A/G|Heterozygous C/T", "H", x))
colnames(KASPgenotypes) <- c("ID", "UNC28010943", "UNC24184030", "UNCHS041907", "JAX00063853", "UNCHS043909", "UNC5812781", "UNC25805470", "UNCHS030444", "UNCHS041714", "UNCHS019508", "UNC27568354")
KASPgenotypes <- t(KASPgenotypes)
colnames(KASPgenotypes) <- KASPgenotypes["ID",]
KASPgenotypes <- KASPgenotypes[-1,]
genotypes[,colnames(KASPgenotypes)] <- "NA"
for (x in 1:nrow(genotypes)){
  if (rownames(genotypes[x,]) %in% rownames(KASPgenotypes)){
    new <- cbind(genotypes[x,1:200], data.frame(as.list(KASPgenotypes[rownames(genotypes[x,]),])))
    colnames(new) <- colnames(genotypes)
    new <- sapply(new, as.character)
    genotypes[x,] <- new
  }
}
genotypes <- genotypes[, order(names(genotypes))]
#write.table(genotypes, file = "genotypesComplete.txt", sep = "\t", quote = FALSE)
dim(genotypes)
dim(phenotypes)

# Convert genotypes to numerical values to map using an additive model and dom model or both
numgeno <- matrix(NA, nrow(genotypes), ncol(genotypes), dimnames=list(rownames(genotypes), colnames(genotypes)))
for(x in 1:nrow(genotypes)){
  h1 <- "A"
  het <- "H"
  h2 <- "B"
  numgeno[x, which(genotypes[x, ] == h1)] <- -1
  numgeno[x, which(genotypes[x, ] == het)] <- 0
  numgeno[x, which(genotypes[x, ] == h2)] <- 1
}


phenotypes <- phenotypes[colnames(genotypes),]
phenotypes <- phenotypes[which(phenotypes[, "Sex"] == "m"),]
numgeno <- numgeno[,rownames(phenotypes)]

# map using an additive model
phenonames <- colnames(phenotypes)[3:34]
pmatrixADD <- matrix(NA, nrow(numgeno), length(phenonames), dimnames= list(rownames(numgeno), phenonames))
for (pname in phenonames){
  cat(pname, " ", "\n")
  pvalues <- apply(numgeno, 1, function(numgeno) {
    numgenoAddd <- as.numeric(unlist(numgeno))
    mdata <- data.frame(cbind(pheno = phenotypes[, pname], sex = phenotypes[, "Sex"], wg = phenotypes[,"WG"], A = numgenoAddd))
    isNA <- which(apply(apply(mdata,1,is.na),2,any))
    if (length(isNA) > 0) mdata <- mdata[-isNA,]
	lm0 <- lm(pheno ~ 1 + sex, data = mdata)
    mmodel <- lm(pheno ~ A + sex, data = mdata)	
    return(anova(mmodel, lm0)[["Pr(>F)"]][2])	
  })
  pmatrixADD[names(pvalues), pname] <- pvalues
}
lodmatrixADDrow <- -log10(pmatrixADD)

lodmatrixADD <- lodmatrixADDrow[rownames(annotation),]
lodmatrixADDannot <- cbind(annotation, lodmatrixADD)

#write.table(lodmatrixADDannot, file = "lodmatrixADD_BW_second.txt", sep = "\t", quote = FALSE)

# some plots
dataset <- lodmatrixADDannot[,grep("D", colnames(lodmatrixADDannot))]
dataset <- cbind(lodmatrixADDannot[,1:3], dataset)
chr15 <- dataset[which(dataset[,"Chromosome"] == 15),]
chr16 <- dataset[which(dataset[,"Chromosome"] == 16),]

par(cex.lab=1.5, cex.main = 1.8, cex.axis = 1.6)
par(mfrow = c(1,2))
#mat <- matrix(c(1,1,2,3), 2, 2, byrow = TRUE)
#layout(mat, widths = rep.int(3, ncol(mat)))
# lodcurve plots for each time point
plot(main = "QTL profile body weight [Chr 15]", c(min(as.numeric(chr15[, "Position"])), max(as.numeric(chr15[, "Position"]))), c(0,9), ylab = "-log10 [pvalue]", xlab = "Position [mb]", las = 2, t = "n", xaxt = "n")
  points(x = as.numeric(chr15[,"Position"]), y = chr15[,"D70"] , type = "l", col="gray1", lwd = 1)
  points(x = as.numeric(chr15[,"Position"]), y = chr15[,"D112"] , type = "l", col="gray45", lwd = 1)
  points(x = as.numeric(chr15[,"Position"]), y = chr15[,"D126"] , type = "l", col="gray60", lwd = 1)
  points(x = as.numeric(chr15[,"Position"]), y = chr15[,"D172"] , type = "l", col="gray90", lwd = 1)
  abline(h=4.7, col="orange")
  abline(h=4.2, col="orange", lty = 2)
  axis(1, at = c(0,25000000, 50000000, 75000000, 100000000, 125000000, 150000000), c("0", "25", "50", "75", "100", "125", "150"))
  legend("topright",
  legend = c("Week 10", "Week 16","Week 18", "Week 25"),
    col = c("gray1", "gray45", "gray60", "gray90"),
    pch = 15,
    pt.cex = 1.7,
    cex = 1, bty = "n" ,
    text.col = "black",
	lwd = c(1,1,1,1)
	)

#mat <- matrix(c(1,1,2,3), 2, 2, byrow = TRUE)
#layout(mat, widths = rep.int(3, ncol(mat)))
# lodcurve plots for each time point
plot(main = "QTL profile body weight [Chr 16]", c(min(as.numeric(chr16[, "Position"])), max(as.numeric(chr16[, "Position"]))), c(0,12), ylab = "-log10 [pvalue]", xlab = "Position [mb]", las = 2, t = "n", xaxt = "n")
  points(x = as.numeric(chr16[,"Position"]), y = chr16[,"D70"] , type = "l", col="gray1", lwd = 1)
  points(x = as.numeric(chr16[,"Position"]), y = chr16[,"D112"] , type = "l", col="gray45", lwd = 1)
  points(x = as.numeric(chr16[,"Position"]), y = chr16[,"D126"] , type = "l", col="gray60", lwd = 1)
  points(x = as.numeric(chr16[,"Position"]), y = chr16[,"D172"] , type = "l", col="gray90", lwd = 1)
  abline(h=4.7, col="orange")
  abline(h=4.2, col="orange", lty = 2)
  axis(1, at = c(0,25000000, 50000000, 75000000, 100000000, 125000000, 150000000), c("0", "25", "50", "75", "100", "125", "150"))
  legend("topright",
  legend = c("Week 10", "Week 16","Week 18", "Week 25"),
    col = c("gray1", "gray45", "gray60", "gray90"),
    pch = 15,
    pt.cex = 1.7,
    cex = 1, bty = "n" ,
    text.col = "black",
	lwd = c(1,1,1,1)
	)
	
# Effect plots for top marker chr 16
topmarker <- t(genotypes["UNCHS041907",])
genopheno <- cbind(topmarker, phenotypes[,phenonames])
genopheno[,1] <- gsub("A" ,"S1", genopheno[,1])
genopheno[,1] <- gsub("H" ,"HET", genopheno[,1])
genopheno[,1] <- gsub("B", "S2", genopheno[,1])
colnames(genopheno) <- c("Genotype", 21, 28, 35, 42, 49, 56, 63, 70, 77, 84, 91, 98, 105, 112, 119, 125, 126, 133, 139, 140, 142, 144, 147, 150, 154, 157, 160, 163, 166, 169, 172, 174)	 
timepoints <- c(28, 42, 56, 70, 84, 98, 112, 126, 140, 160, 174)  
 
plot(main="Body weight at Chr 16 top marker", c(25,180), c(10,60), ylab="[g]", xlab="Age [weeks]", yaxs = "i", las = 2, t = "n", xaxt="n")
  axis(1, at = c(28, 42, 55, 70, 84, 98, 112, 126, 140, 160, 174), c("4", "6", "8", "10", "12", "14", "16", "18", "20", "22", "25"), lwd = 1, cex.axis=1.6)
  meansBFMI <- c()
  meansHET <- c()
  meansB6 <- c()
  for (x in timepoints){
   bptBFMI <- boxplot(at = x+1, as.numeric(as.character(genopheno[which(genopheno == "S1"), as.character(x)]))  ~ genopheno[which(genopheno == "S1"),"Genotype"], width=6, col = "gray20", axes = FALSE, add=TRUE, notch= TRUE, outcex=0.5, boxwex = 1.5) 
   meanBFMI <- bptBFMI$stats[3,] 
   meansBFMI <- c(meansBFMI, meanBFMI)
   bptHET <- boxplot(at = x, as.numeric(as.character(genopheno[which(genopheno == "HET"), as.character(x)]))  ~ genopheno[which(genopheno == "HET"),"Genotype"], width=6, col = "gray50", axes = FALSE, add=TRUE, notch= TRUE, outcex=0.5, boxwex = 1.5) 
   meanHET <- bptHET$stats[3,]
   meansHET <- c(meansHET, meanHET)
   bptB6 <- boxplot(at = x-1, as.numeric(as.character(genopheno[which(genopheno == "S2"), as.character(x)]))  ~ genopheno[which(genopheno == "S2"),"Genotype"], width=6, col = "gray88", axes = FALSE, add=TRUE, notch= TRUE, outcex=0.5, boxwex = 1.5) 
   meanB6 <- bptB6$stats[3,]
   meansB6 <- c(meansB6, meanB6)
  }
  lines(c(28, 42, 55, 70, 84, 98, 112, 126, 140, 160, 174), meansBFMI, col="gray0", lwd=1)
  lines(c(28, 42, 55, 70, 84, 98, 112, 126, 140, 160, 174), meansHET, col="gray50", lwd=1)
  lines(c(28, 42, 55, 70, 84, 98, 112, 126, 140, 160, 174), meansB6, col="gray88", lwd=1) 
  legend("topleft",  
  legend = c("BFMI861-S1", "HET", "BFMI861-S2"), 
   col = c("gray20", "gray50", "gray88"),
   pch = 15, pt.cex = 1.7, cex = 1, bty = "n"
   )

# Effect plots for top marker chr 15
topmarker <- t(genotypes["UNCHS040893",])
genopheno <- cbind(topmarker, phenotypes[,phenonames])
genopheno[,1] <- gsub("B" ,"S1", genopheno[,1])
genopheno[,1] <- gsub("H" ,"HET", genopheno[,1])
genopheno[,1] <- gsub("A", "S2", genopheno[,1])
colnames(genopheno) <- c("Genotype", 21, 28, 35, 42, 49, 56, 63, 70, 77, 84, 91, 98, 105, 112, 119, 125, 126, 133, 139, 140, 142, 144, 147, 150, 154, 157, 160, 163, 166, 169, 172, 174)	 
timepoints <- c(28, 42, 56, 70, 84, 98, 112, 126, 140, 160, 174)  
 
plot(main="Body weight at Chr 15 top marker", c(25,180), c(10,60), ylab="[g]", xlab="Age [weeks]", yaxs = "i", las = 2, t = "n", xaxt="n")
  axis(1, at = c(28, 42, 55, 70, 84, 98, 112, 126, 140, 160, 174), c("4", "6", "8", "10", "12", "14", "16", "18", "20", "22", "25"), lwd = 1, cex.axis=1.6)
  meansBFMI <- c()
  meansHET <- c()
  meansB6 <- c()
  for (x in timepoints){
   bptBFMI <- boxplot(at = x+1, as.numeric(as.character(genopheno[which(genopheno == "S1"), as.character(x)]))  ~ genopheno[which(genopheno == "S1"),"Genotype"], width=6, col = "gray20", axes = FALSE, add=TRUE, notch= TRUE, outcex=0.5, boxwex = 1.5) 
   meanBFMI <- bptBFMI$stats[3,] 
   meansBFMI <- c(meansBFMI, meanBFMI)
   bptHET <- boxplot(at = x, as.numeric(as.character(genopheno[which(genopheno == "HET"), as.character(x)]))  ~ genopheno[which(genopheno == "HET"),"Genotype"], width=6, col = "gray50", axes = FALSE, add=TRUE, notch= TRUE, outcex=0.5, boxwex = 1.5) 
   meanHET <- bptHET$stats[3,]
   meansHET <- c(meansHET, meanHET)
   bptB6 <- boxplot(at = x-1, as.numeric(as.character(genopheno[which(genopheno == "S2"), as.character(x)]))  ~ genopheno[which(genopheno == "S2"),"Genotype"], width=6, col = "gray88", axes = FALSE, add=TRUE, notch= TRUE, outcex=0.5, boxwex = 1.5) 
   meanB6 <- bptB6$stats[3,]
   meansB6 <- c(meansB6, meanB6)
  }
  lines(c(28, 42, 55, 70, 84, 98, 112, 126, 140, 160, 174), meansBFMI, col="gray0", lwd=1)
  lines(c(28, 42, 55, 70, 84, 98, 112, 126, 140, 160, 174), meansHET, col="gray50", lwd=1)
  lines(c(28, 42, 55, 70, 84, 98, 112, 126, 140, 160, 174), meansB6, col="gray88", lwd=1) 
  legend("topleft",  
  legend = c("BFMI861-S1", "HET", "BFMI861-S2"), 
   col = c("gray20", "gray50", "gray88"),
   pch = 15, pt.cex = 1.7, cex = 1, bty = "n"
   )
