# QTL mapping using just the females
#
# copyright (c) 2018-2021 - Brockmann group - HU Berlin Manuel Delpero 
# 
# first written april, 2020

setwd("C:/Users/Manuel/Desktop/AIL_S1xS2/RAWDATA")
#setwd("/home/manuel/AIL_S1xS2/DATA")
genotypes <- read.csv("genotypes.cleaned.txt", header = TRUE, check.names = FALSE, sep="\t", colClasses="character")
phenotypes <- read.csv("allPhenotypes.txt", header = TRUE, check.names = FALSE, sep="\t", row.names=1)
annotation <- read.csv("map.cleaned.txt", header=TRUE, sep="\t", check.names=FALSE)
annotation <- annotation[, c(1,2,3,6)]
colnames(annotation) <- c("Chromosome", "Position", "GenTrain Score", "SNP")
colnames(genotypes) <- gsub("AIL", "" , colnames(genotypes)) 

## QTL mapping using also the genotypes by KASP assay
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
write.table(genotypes, file = "genotypesComplete.txt", sep = "\t", quote = FALSE)

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
phenonames <- colnames(phenotypes[,c(3:57, 61, 62, 65:74)])

# Getting rid of the outliers
outliers <- apply(phenotypes[, phenonames],2, function(x){
  up <- mean(x, na.rm = TRUE) + 3 * sd(x, na.rm = TRUE)
  low <- mean(x, na.rm=TRUE) - 3 * sd(x, na.rm=TRUE)
  x < low | x > up
 })
 
for(x in phenonames){
  idx <- which(outliers[,x])
  if(length(idx) > 0) phenotypes[idx,x] <- NA
}

# Adjust the tissue weight for the total bodyweight
tissues <- colnames(phenotypes[, c(46:55)])
for (x in tissues){
  phenotypes[, x] <- phenotypes[, x] / phenotypes[, "Gewicht"]
}
phenotypes <- phenotypes[colnames(genotypes),]
# Get the females
phenotypes <- phenotypes[which(phenotypes[, "Sex"] == "f"),]
numgeno <- numgeno[rownames(KASPgenotypes) , rownames(phenotypes)]

# Dom dev + Add model without using the sum of LODS and no covariates (real dom)
pmatrixADDDOM <- matrix(NA, nrow(numgeno), length(phenonames), dimnames= list(rownames(numgeno), phenonames))
for (pname in phenonames){
  cat(pname, " ", "\n")
  pvalues <- apply(numgeno, 1, function(numgeno) {
    numgenoDomm <- as.numeric(as.numeric(unlist(numgeno)) != 0)
    numgenoAddd <- as.numeric(unlist(numgeno))
    mdata <- data.frame(cbind(pheno = phenotypes[, pname], sex = phenotypes[, "Sex"], A = numgenoAddd, D = numgenoDomm))
    isNA <- which(apply(apply(mdata,1,is.na),2,any))
    if (length(isNA) > 0) mdata <- mdata[-isNA, ]
      lm0 <- lm(pheno ~ 1 + sex, data = mdata)
      mmodel <- lm(pheno ~ D + A + sex, data = mdata)
      return(anova(mmodel, lm0)[["Pr(>F)"]][2])	
  })
  pmatrixADDDOM[names(pvalues), pname] <- pvalues
}
lodmatrixADDDOM <- -log10(pmatrixADDDOM)
write.table(lodmatrixADDDOM, file = "lodmatrixADDDOM_nosumF.txt", quote = FALSE, sep = "\t")
lodannotmatrix <- cbind(annotation[rownames(lodmatrixADDDOM), ], lodmatrixADDDOM)

# Dominance dev model
pmatrixDOM <- matrix(NA, nrow(numgeno), length(phenonames), dimnames= list(rownames(numgeno), phenonames))
for (pname in phenonames){
  cat(pname, " ", "\n")
  pvalues <- apply(numgeno, 1, function(numgeno) {
    numgenoDomm <- as.numeric(as.numeric(unlist(numgeno)) != 0)
    mdata <- data.frame(cbind(pheno = phenotypes[, pname], sex = phenotypes[, "Sex"], D = numgenoDomm))
    isNA <- which(apply(apply(mdata,1,is.na),2,any))
    if (length(isNA) > 0) mdata <- mdata[-isNA, ]
      lm0 <- lm(pheno ~ 1 + sex, data = mdata)
      mmodel <- lm(pheno ~ D + sex, data = mdata)
      return(anova(mmodel, lm0)[["Pr(>F)"]][2])	
  })
  pmatrixDOM[names(pvalues), pname] <- pvalues
}
lodmatrixDOM <- -log10(pmatrixDOM)
write.table(lodmatrixDOM, file = "lodmatrixDOM_nosumF.txt", quote = FALSE, sep = "\t")

# Additive model
pmatrixADD <- matrix(NA, nrow(numgeno), length(phenonames), dimnames= list(rownames(numgeno), phenonames))
for (pname in phenonames){
  cat(pname, " ", "\n")
  pvalues <- apply(numgeno, 1, function(numgeno) {
    numgenoAddd <- as.numeric(unlist(numgeno))
    mdata <- data.frame(cbind(pheno = phenotypes[, pname], sex = phenotypes[, "Sex"], A = numgenoAddd))
    isNA <- which(apply(apply(mdata,1,is.na),2,any))
    if (length(isNA) > 0) mdata <- mdata[-isNA, ]
      lm0 <- lm(pheno ~ 1 + sex, data = mdata)
      mmodel <- lm(pheno ~ A + sex, data = mdata)
      return(anova(mmodel, lm0)[["Pr(>F)"]][2])	
  })
  pmatrixADD[names(pvalues), pname] <- pvalues
}
lodmatrixADD <- -log10(pmatrixADD)     # QTL on chr7 for liver weight and a small QTL on chr 16 for bodyweight
write.table(lodmatrixADD, file = "lodmatrixADDF.txt", quote = FALSE, sep = "\t")