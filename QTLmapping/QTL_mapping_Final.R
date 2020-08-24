# Modelling using AIL BFMI S1xS2
#
# copyright (c) 2018-2021 - Brockmann group - HU Berlin Manuel Delpero & Danny Arends
# 
# first written december, 2019

setwd("C:/Users/Manuel/Desktop/AIL_S1xS2/RAWDATA")
#setwd("/home/manuel/AIL_S1xS2/DATA")
genotypes <- read.csv("genotypes.cleaned.txt", header = TRUE, check.names = FALSE, sep="\t", colClasses="character")
#phenotypes <- read.csv("allPhenotypes.txt", header = TRUE, check.names = FALSE, sep="\t", row.names=1)
phenotypes <- read.csv("GlucSeries.txt", header = TRUE, check.names = FALSE, sep="\t", row.names=1)
rownames(phenotypes) <- gsub("V 888-", "", rownames(phenotypes))
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
#write.table(genotypes, file = "genotypesComplete.txt", sep = "\t", quote = FALSE)


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

#phenonames <- colnames(phenotypes[,c(3:57, 61, 62, 65:74)])
phenonames <- colnames(phenotypes)
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

# Calculate lost of weight during diet switch
switchLost <- phenotypes[,"D125"] - phenotypes[,"D126"]
phenotypes <- cbind(phenotypes, switchLost)
phenotypes <- phenotypes[colnames(genotypes),]
phenonames <- colnames(phenotypes[,c(3:57, 61, 62, 65:75)])
#write.table(phenotypes, file = "PhenotypesComplete.txt", sep = "\t", quote = FALSE, row.names = TRUE)

sex <- read.csv("allPhenotypes.txt", header = TRUE, check.names = FALSE, sep="\t", row.names=1)
sex <- sex[rownames(phenotypes), "Sex"]
phenotypes <- cbind(sex, phenotypes)
colnames(phenotypes) <- c("Sex","Gluc42" , "Gluc70" , "Gluc98" , "Gluc125" ,"Gluc139" ,"Gluc147" ,"Gluc154" ,"Gluc157" ,"Gluc160" ,"Gluc163" ,"Gluc166" ,"Gluc169" ,"Gluc172")
phenotypes <- phenotypes[colnames(genotypes),]


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
write.table(lodmatrixADDDOM, file = "lodmatrixADDDOMCompleteGlucSeries.txt", quote = FALSE, sep = "\t")
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
write.table(lodmatrixDOM, file = "lodmatrixDOMCompleteGlucSeries.txt", quote = FALSE, sep = "\t")

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
lodmatrixADD <- -log10(pmatrixADD)
write.table(lodmatrixADD, file = "lodmatrixADDCompleteGlucSeries.txt", quote = FALSE, sep = "\t")

# Figure out the variance explained by each QTL considering the direction of the effect
lodmatrixDOM <- read.csv("lodmatrixDOMComplete.txt", header = TRUE, sep = "\t", check.names = FALSE)
lodmatrixADD <- read.csv("lodmatrixADDComplete.txt", header = TRUE, sep = "\t", check.names = FALSE)
lodmatrixADDDOM <- read.table("lodmatrixADDDOMComplete.txt", header = TRUE, sep = "\t", check.names = FALSE)

# get just the top markers for the phenotypes
phenos <- c("Gon", "Leber")
TopMarkersGon <- c("UNC5791802", "UNC20599050", "UNC25806117", "UNC27977378")
TopMarkersLiver <- "UNCHS043909"

# Var explained by the Dominance or additive effect
varExplained <- matrix(NA, 5, 2) 
rownames(varExplained) <- c(TopMarkersGon, TopMarkersLiver)
colnames(varExplained) <- phenos
for (pname in phenos){
  for (x in rownames(varExplained)){
    if ((lodmatrixADD[x, pname] > lodmatrixDOM[x, pname]) && (lodmatrixADD[x, pname] > lodmatrixADDDOM[x, pname])){ 
	  pheno <- phenotypes[, pname]
	  numgenoAddd <- as.numeric(unlist(numgeno[x,]))
      mdata <- data.frame(cbind(pheno, sex = phenotypes[, "Sex"], A = numgenoAddd))
      isNA <- which(apply(apply(mdata,1,is.na),2,any))
      if (length(isNA) > 0) mdata <- mdata[-isNA, ]
      mmodel <- lm(pheno ~ sex + A, data = mdata)
      sumSQ <- anova(mmodel)["A", "Sum Sq"]
	  var <- sumSQ / sum((pheno - mean(pheno, na.rm=TRUE))^2, na.rm=TRUE)
      var <- round(var * 100, digits=1)	  
    }else if ((lodmatrixDOM[x, pname] > lodmatrixADD[x, pname]) && (lodmatrixDOM[x, pname] > lodmatrixADDDOM[x, pname])){
      pheno <- phenotypes[, pname]
	  numgenoDomm <- as.numeric(as.numeric(unlist(numgeno[x,])) != 0)
      mdata <- data.frame(cbind(pheno, sex = phenotypes[, "Sex"], D = numgenoDomm))
      isNA <- which(apply(apply(mdata,1,is.na),2,any))
      if (length(isNA) > 0) mdata <- mdata[-isNA, ]
      mmodel <- lm(pheno ~ sex + D, data = mdata)
      sumSQ <- anova(mmodel)["D", "Sum Sq"]
	  var <- sumSQ / sum((pheno - mean(pheno, na.rm=TRUE))^2, na.rm=TRUE)
	  var <- round(var * 100, digits=1)	  
	}else var <- NA
    varExplained[x,pname] <- var
  }
}
	 

# Function to calculate the variance explained by the additive and dominance deviation effect at the same time giving the marker name and the phenotype as arguments
DomAddVariance <- function(marker, pname){
  numgenoAddd <- as.numeric(unlist(numgeno[marker,]))
  numgenoDomm <- as.numeric(as.numeric(unlist(numgeno[marker,])) != 0)
  mdata <- data.frame(cbind(pheno = phenotypes[, pname], sex = phenotypes[, "Sex"], A = numgenoAddd, D = numgenoDomm))
  isNA <- which(apply(apply(mdata,1,is.na),2,any))
  if (length(isNA) > 0) mdata <- mdata[-isNA, ]
  mmodel <- lm(pheno ~ A + D, data = mdata)
  sumSQ <- anova(mmodel)[c("A","D"), "Sum Sq"]
  pheno <- phenotypes[, pname]
  var <- sumSQ / sum((pheno - mean(pheno, na.rm=TRUE))^2, na.rm=TRUE)
  var <- round(var * 100, digits=1)
  return(paste(var[1], "Additive", var[2], "Dominance deviation"))
}

# Get variance giving manually the genotypes coded as we want
getVarianceExplained <- function(genotypes, phenotypes, pheno.col = "d77", marker = "gUNC5036315"){
  genotype <- as.factor(t(genotypes[marker,]))
  littersize <- as.factor(phenotypes[, "WG"])
  subfamily <- as.factor(phenotypes[, "Mutter"])
  sex <- as.factor(phenotypes[, "Sex"])
  phenotype <- as.numeric(phenotypes[, pheno.col])
  model <- lm(phenotype ~ genotype)
  tryCatch(res  <- anova(model), error = function(e){ res <<- NA })
  cat(names(model$coefficients),"\n")
  cat(model$coefficients,"\n")
  varExplained  <- res[, "Sum Sq"] / sum((phenotype - mean(phenotype, na.rm=TRUE))^2, na.rm=TRUE)
  names(varExplained) <- c("marker", "others")
  return(round(varExplained * 100, digits=1))
}

getVarianceExplained(genotypes, phenotypes, pheno.col = "Gon", marker = "UNCHS041907")