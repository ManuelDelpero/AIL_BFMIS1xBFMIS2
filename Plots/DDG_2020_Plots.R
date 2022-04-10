# Plots for DDG 2022

setwd("C:/Users/Manuel/Desktop/AIL_S1xS2/RAWDATA")

genotypes <- read.csv("genotypes.cleaned.txt", header = TRUE, check.names = FALSE, sep="\t", colClasses="character")
GlucSeries <- read.csv("GlucSeries.txt", header = TRUE, check.names = FALSE, sep="\t", row.names = 1)
rownames(GlucSeries) <- gsub("V 888-", "", rownames(GlucSeries))
Plasma <- read.csv("Plasma.txt", header = TRUE, check.names = FALSE, sep="\t", row.names=1)
Plasma <- Plasma[,c(1,2,3,4,5,6)]
Pankreas_Insulin = read.csv("Pankreas_Insulin.txt", sep = "\t", check.names = FALSE, header = TRUE, row.names = 1)
phenotypes <- read.csv("phenotypesCompleteAll.txt", header = TRUE, check.names = FALSE, sep="\t", row.names=1)
annotation <- read.csv("map.cleaned.txt", header=TRUE, sep="\t", check.names=FALSE)
annotation <- annotation[, c(1,2,3,6)]
colnames(annotation) <- c("Chromosome", "Position", "GenTrain Score", "SNP")
colnames(genotypes) <- gsub("AIL", "" , colnames(genotypes)) 
dim(genotypes)
dim(phenotypes)

GlucSeries <- GlucSeries[rownames(phenotypes),]
Plasma <- Plasma[rownames(phenotypes),]
phenotypes <- cbind(phenotypes, GlucSeries)
Pankreas_Insulin <- Pankreas_Insulin[rownames(phenotypes),]
phenotypes <- cbind(phenotypes, Pankreas_Insulin, Plasma)

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

# Take only the males
phenotypes <- phenotypes[which(phenotypes[, "Sex"] == "m"),] 
Plasma <- Plasma[rownames(phenotypes),]
phenotypes <- cbind(phenotypes, Plasma)
numgeno <- numgeno[,rownames(phenotypes)]
#phenotypes <- phenotypes[colnames(genotypes),]
wg <- phenotypes[, "WG"]
grandmother <- phenotypes[, "Grandma"]
Leber <- phenotypes[, "Leber"]
Gon <- phenotypes[, "Gon"]

phenonames <- c("Pankreas_Insulin", "Insulin")

# map 
phenotypes[,"Insulin"] <- log2(phenotypes[,"Insulin"])
phenotypes[,"Pankreas_Insulin"] <- log2(phenotypes[,"Pankreas_Insulin"])

# Additive model
pmatrixADD <- matrix(NA, nrow(numgeno), length(phenonames), dimnames= list(rownames(numgeno), phenonames))
for (pname in phenonames){
  cat(pname, " ", "\n")
  pvalues <- apply(numgeno, 1, function(numgeno) {
    numgenoAddd <- as.numeric(unlist(numgeno))
    mdata <- data.frame(cbind(pheno = phenotypes[, pname], sex = phenotypes[, "Sex"], wg = phenotypes[,"WG"], A = numgenoAddd))
    isNA <- which(apply(apply(mdata,1,is.na),2,any))
    if (length(isNA) > 0) mdata <- mdata[-isNA,]
	lm0 <- lm(pheno ~ 1, data = mdata)
    mmodel <- lm(pheno ~ A, data = mdata)	
    return(anova(mmodel, lm0)[["Pr(>F)"]][2])	
  })
  pmatrixADD[names(pvalues), pname] <- pvalues
}
lodmatrixADDrow <- -log10(pmatrixADD)

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

lodmatrixDOMrow <- lodmatrixDOMrow[rownames(annotation),]
lodmatrixADDrow <- lodmatrixADDrow[rownames(annotation),]
lodmatrixADDDOMrow <- lodmatrixADDDOMrow[rownames(annotation),]


# chr 17
lodannotmatrix <- cbind(annotation[rownames(lodmatrixADDrow), ], lodmatrixADDrow)
chr3 <- lodannotmatrix[which(lodannotmatrix[,"Chromosome"] == 17),]
datasetRow <- chr3[, c("Chromosome", "Position", "Insulin", "Pankreas_Insulin")]
datasetRow <- datasetRow[which(dataset[,"Position"] < 26389801),]
datasetRow[which.max(datasetRow[,4]),4] = 4

plot(main = "Gatlgq", c(min(as.numeric(datasetRow[, "Position"])), max(as.numeric(datasetRow[, "Position"]))), c(-0.5,12), ylab = "LOD score", xlab = "Position [mb]", las = 2, t = "n", xaxt = "n")
  points(x = as.numeric(datasetRow[,"Position"]), y = datasetRow[,"Insulin"], t ='l', col="cornflowerblue", lwd = 2)
  lines(x = as.numeric(datasetRow[,"Position"]), y = datasetRow[,"Pankreas_Insulin"], t = "l", col="red", lwd = 2)
  points(x = as.numeric(datasetRow[,"Position"]), y = rep(-0.5,nrow(datasetRow)), pch = "|", col = "black", lwd = 1)
  abline(h=4.7, col="orange" )
  abline(h=4.2, col="orange", lty = 2)
  axis(1, at = c(3250036, 10000000, 20000000, 30000000), c("3.2", "10", "20", "30"))
  legend("topleft",
  legend = c("Pankreas Insulin", "Plasma Insulin"),
    bty = "n",
    col = c("cornflowerblue", "red"),
    lty=c(1,1,1,1),
    pt.cex = 1.2,
    pt.bg = "lightsteelblue1",
    cex = 1.2,
	lwd = c(3,3,3),
    text.col = "black")