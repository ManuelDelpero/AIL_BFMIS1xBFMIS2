# Modelling using AIL BFMI S1xS2
#
# copyright (c) 2018-2021 - Brockmann group - HU Berlin Manuel Delpero & Danny Arends
# 
# first written december, 2019

setwd("C:/Users/Manuel/Desktop/AIL_S1xS2/RAWDATA")

genotypes <- read.csv("genomatrix.clean.txt", header = TRUE, check.names = FALSE, sep="\t", colClasses="character")
phenotypes <- read.csv("allPhenotypes_final.txt", header = TRUE, check.names = FALSE, sep="\t", row.names=1)
markerannot <- read.csv("map.cleaned.txt", header=TRUE, sep="\t", check.names=FALSE)
markerannot <- markerannot[, c(1,2,3,6)]
colnames(markerannot) <- c("Chromosome", "Position", "GenTrain Score", "SNP")
colnames(genotypes) <- gsub("V 888-", "" , colnames(genotypes)) 
phenotypes <- phenotypes[colnames(genotypes),]

# Covariates we could/need to include in the model, we test them on their pvalue
wg <- phenotypes[, "WG"]
grandmother <- phenotypes[, "Grandma"]

# Convert genotypes to numerical values to map using an additive model and dom model or both
numgeno <- matrix(NA, nrow(genotypes), ncol(genotypes), dimnames=list(rownames(genotypes), colnames(genotypes)))
for(x in 1:nrow(genotypes)){
  alleles <- na.omit(unique(unlist(strsplit(as.character(genotypes[x,]), ""))))
  h1 <- paste0(alleles[1], alleles[1])
  het <- c(paste0(alleles[1], alleles[2]), paste0(alleles[2], alleles[1]))
  h2 <- paste0(alleles[2], alleles[2])
  numgeno[x, which(genotypes[x, ] == h1)] <- -1
  numgeno[x, which(genotypes[x, ] %in% het)] <- 0
  numgeno[x, which(genotypes[x, ] == h2)] <- 1
}
phenonames <- colnames(phenotypes[,c(3:57)])
 
# AD model
# QTL mapping testing if we need to include the covariates
pmatrix <- matrix(NA, nrow(genotypes), length(phenonames), dimnames= list(rownames(genotypes), phenonames))
for (pname in phenonames){
  pheno <- phenotypes[, pname]
  p.gmother <- anova(lm(pheno ~ phenotypes[, "Grandma"]))["Pr(>F)"][1,]
  p.wg <- anova(lm(pheno ~ phenotypes[, "WG"]))["Pr(>F)"][1,]
  myfactors <- c()
  if(p.gmother < 0.05) myfactors <- c(myfactors, "grandmother")
  if(p.wg < 0.05) myfactors <- c(myfactors, "wg")
  myformula <- paste0("pheno ~ ", paste(c(myfactors, "numgeno"), collapse = " + "))
  cat(pname, " ", "\n")
  
  pvalues <- apply(numgeno, 1, function(numgeno, pheno, wg, grandmother, myformula) {
    numgeno <- as.numeric(unlist(numgeno))
    mmodel <- lm(formula(myformula))
    return(anova(mmodel)["Pr(>F)"]["numgeno",])
  }, pheno = phenotypes[,pname], wg = as.factor(phenotypes[,"WG"]), grandmother = as.factor(phenotypes[,"Grandma"]), myformula = myformula)
  pmatrix[names(pvalues), pname] <- pvalues
}
lodmatrixAD <- -log10(pmatrix)

# Dom model
pmatrix <- matrix(NA, nrow(genotypes), length(phenonames), dimnames= list(rownames(genotypes), phenonames))
for (pname in phenonames){
  pheno <- phenotypes[, pname]
  p.gmother <- anova(lm(pheno ~ phenotypes[, "Grandma"]))["Pr(>F)"][1,]
  p.wg <- anova(lm(pheno ~ phenotypes[, "WG"]))["Pr(>F)"][1,]
  myfactors <- c()
  if(p.gmother < 0.05) myfactors <- c(myfactors, "grandmother")
  if(p.wg < 0.05) myfactors <- c(myfactors, "wg")
  myformula <- paste0("pheno ~ ", paste(c(myfactors, "numgeno"), collapse = " + "))
  cat(pname, " ", "\n")
  
  pvalues <- apply(numgeno, 1, function(numgeno, pheno, wg, grandmother, myformula) {
    numgeno <- as.numeric(as.numeric(unlist(numgeno)) != 0)
    mmodel <- lm(formula(myformula))
    return(anova(mmodel)["Pr(>F)"]["numgeno",])
  }, pheno = phenotypes[,pname], wg = as.factor(phenotypes[,"WG"]), grandmother = as.factor(phenotypes[,"Grandma"]), myformula = myformula)
  pmatrix[names(pvalues), pname] <- pvalues
}
lodmatrixDOM <- -log10(pmatrix)

# Dom + Add model
pmatrixADD <- matrix(NA, nrow(genotypes), length(phenonames), dimnames= list(rownames(genotypes), phenonames))
pmatrixDOM <- matrix(NA, nrow(genotypes), length(phenonames), dimnames= list(rownames(genotypes), phenonames))
for (pname in phenonames){
  pheno <- phenotypes[, pname]
  p.gmother <- anova(lm(pheno ~ phenotypes[, "Grandma"]))["Pr(>F)"][1,]
  p.wg <- anova(lm(pheno ~ phenotypes[, "WG"]))["Pr(>F)"][1,]
  myfactors <- c()
  if(p.gmother < 0.05) myfactors <- c(myfactors, "grandmother")
  if(p.wg < 0.05) myfactors <- c(myfactors, "wg")
  myformula <- paste0("pheno ~ ", paste(c(myfactors, "numgeno"), collapse = " + "))
  cat(pname, " ", "\n")
  
  pvalues <- apply(numgeno, 1, function(numgeno, pheno, wg, grandmother, myformula) {
    numgenoDom <- as.numeric(as.numeric(unlist(numgeno)) != 0)
	numgenoAdd <- as.numeric(unlist(numgeno))
	myformula <- paste0("pheno ~ ", paste(c(myfactors, "numgenoDom", "numgenoAdd"), collapse = " + "))
    mmodel <- lm(formula(myformula))
    return(anova(mmodel)["Pr(>F)"][c("numgenoDom", "numgenoAdd"),])
  }, pheno = phenotypes[,pname], wg = as.factor(phenotypes[,"WG"]), grandmother = as.factor(phenotypes[,"Grandma"]), myformula = myformula)
  pmatrixDOM[names(pvalues[1,]), pname] <- pvalues[1,]
  pmatrixADD[names(pvalues[2,]), pname] <- pvalues[2,]
}
lodmatrixDOM <- -log10(pmatrixDOM)
lodmatrixADD <- -log10(pmatrixADD)
<<<<<<< HEAD
lodmatrixADDDOM <- lodmatrixDOM + lodmatrixADD
=======
lodmatrixADDDOM <- lodmatrixDOM + lodmatrixADD
>>>>>>> 43d1ba8b6a1aaa16af5664b369118b524373dff2
