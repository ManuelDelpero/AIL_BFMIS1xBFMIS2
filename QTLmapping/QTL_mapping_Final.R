# Modelling using AIL BFMI S1xS2
#
# copyright (c) 2018-2021 - Brockmann group - HU Berlin Manuel Delpero & Danny Arends
# 
# first written december, 2019

setwd("C:/Users/Manuel/Desktop/AIL_S1xS2/RAWDATA")
#setwd("/home/manuel/AIL_S1xS2/DATA")
genotypes <- read.csv("genotypes.cleaned.txt", header = TRUE, check.names = FALSE, sep="\t", colClasses="character")
phenotypes <- read.csv("allPhenotypes.txt", header = TRUE, check.names = FALSE, sep="\t", row.names=1)
annotation <- read.csv("map.cleaned.txt", header=TRUE, sep="\t", check.names=FALSE)
annotation <- annotation[, c(1,2,3,6)]
colnames(annotation) <- c("Chromosome", "Position", "GenTrain Score", "SNP")
colnames(genotypes) <- gsub("AIL", "" , colnames(genotypes)) 

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

# Dom + Add model without using the sum of LODS and no covariates
pmatrixADDDOM <- matrix(NA, nrow(genotypes), length(phenonames), dimnames= list(rownames(genotypes), phenonames))
for (pname in phenonames[58:68]){
  cat(pname, " ", "\n")
  pvalues <- apply(numgeno, 1, function(numgeno) {
    numgenoDomm <- as.numeric(as.numeric(unlist(numgeno)) != 0)
	numgenoAddd <- as.numeric(unlist(numgeno))
    mdata <- data.frame(cbind(pheno = phenotypes[, pname], A = numgenoAddd, D = numgenoDomm))
    isNA <- which(apply(apply(mdata,1,is.na),2,any))
    if (length(isNA) > 0) mdata <- mdata[-isNA, ]
	lm0 <- lm(pheno ~ 1, data = mdata)
    mmodel <- lm(pheno ~ D + A, data = mdata)
    return(anova(mmodel, lm0)[["Pr(>F)"]][2])	
  })
  pmatrixADDDOM[names(pvalues), pname] <- pvalues
}
lodmatrixADDDOM <- -log10(pmatrixADDDOM)
write.table(lodmatrixADDDOM, file = "lodmatrixADDDOM_nosum.txt", quote = FALSE, sep = "\t")
lodannotmatrix <- cbind(annotation[rownames(lodmatrixADDDOM), ], lodmatrixADDDOM)

# Dominance variation model
pmatrixDOM <- matrix(NA, nrow(genotypes), length(phenonames), dimnames= list(rownames(genotypes), phenonames))
for (pname in phenonames){
  cat(pname, " ", "\n")
  pvalues <- apply(numgeno, 1, function(numgeno) {
    numgenoDomm <- as.numeric(as.numeric(unlist(numgeno)) != 0)
    mdata <- data.frame(cbind(pheno = phenotypes[, pname], D = numgenoDomm))
    isNA <- which(apply(apply(mdata,1,is.na),2,any))
    if (length(isNA) > 0) mdata <- mdata[-isNA, ]
	lm0 <- lm(pheno ~ 1, data = mdata)
    mmodel <- lm(pheno ~ D, data = mdata)
    return(anova(mmodel, lm0)[["Pr(>F)"]][2])	
  })
  pmatrixDOM[names(pvalues), pname] <- pvalues
}
lodmatrixDOM <- -log10(pmatrixDOM)
write.table(lodmatrixDOM, file = "lodmatrixDOM_nosum.txt", quote = FALSE, sep = "\t")

# Additive model
pmatrixADD <- matrix(NA, nrow(genotypes), length(phenonames), dimnames= list(rownames(genotypes), phenonames))
for (pname in phenonames){
  cat(pname, " ", "\n")
  pvalues <- apply(numgeno, 1, function(numgeno) {
	numgenoAddd <- as.numeric(unlist(numgeno))
    mdata <- data.frame(cbind(pheno = phenotypes[, pname], A = numgenoAddd))
    isNA <- which(apply(apply(mdata,1,is.na),2,any))
    if (length(isNA) > 0) mdata <- mdata[-isNA, ]
	lm0 <- lm(pheno ~ 1, data = mdata)
    mmodel <- lm(pheno ~ A, data = mdata)
    return(anova(mmodel, lm0)[["Pr(>F)"]][2])	
  })
  pmatrixADD[names(pvalues), pname] <- pvalues
}
lodmatrixADD <- -log10(pmatrixADD)
write.table(lodmatrixADD, file = "lodmatrixADD.txt", quote = FALSE, sep = "\t")

