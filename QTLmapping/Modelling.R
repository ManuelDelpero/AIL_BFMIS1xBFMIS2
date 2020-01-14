# Modelling using AIL BFMI S1xS2
#
# copyright (c) 2018-2021 - Brockmann group - HU Berlin Manuel Delpero & Danny Arends
# 
# first written december, 2019

setwd("C:/Users/Manuel/Desktop/AIL_S1xS2/RAWDATA")

genotypes <- read.csv("genotypes.cleaned.txt", header = TRUE, check.names = FALSE, sep="\t", colClasses="character")
phenotypes <- read.csv("allPhenotypes.txt", header = TRUE, check.names = FALSE, sep="\t", row.names=1)
annotation <- read.csv("map.cleaned.txt", header=TRUE, sep="\t", check.names=FALSE)
annotation <- annotation[, c(1,2,3,6)]
colnames(annotation) <- c("Chromosome", "Position", "GenTrain Score", "SNP")
colnames(genotypes) <- gsub("AIL", "" , colnames(genotypes)) 
phenotypes <- phenotypes[colnames(genotypes),]

# Covariates we could/need to include in the model, we test them on their pvalue
wg <- phenotypes[, "WG"]
grandmother <- phenotypes[, "Grandma"]

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

# Dom + Add model with sum of LODS
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
lodmatrixADDDOM <- lodmatrixDOM + lodmatrixADD
write.table(lodmatrixDOM, file = "lodmatrixDOM.txt", quote = FALSE, sep = "\t")
write.table(lodmatrixADD, file = "lodmatrixADD.txt", quote = FALSE, sep = "\t")
write.table(lodmatrixADDDOM, file = "lodmatrixADDDOM.txt", quote = FALSE, sep = "\t")

# Dom + Add model without using the sum of LODS
pmatrixADDDOM <- matrix(NA, nrow(genotypes), length(phenonames), dimnames= list(rownames(genotypes), phenonames))
for (pname in phenonames){
  pheno <- phenotypes[, pname]
  p.gmother <- anova(lm(pheno ~ phenotypes[, "Grandma"]))["Pr(>F)"][1,]
  p.wg <- anova(lm(pheno ~ phenotypes[, "WG"]))["Pr(>F)"][1,]
  myfactors <- c()
  if(p.gmother < 0.05) myfactors <- c(myfactors, "grandmother")
  if(p.wg < 0.05) myfactors <- c(myfactors, "wg")
  if((p.wg > 0.05) && (p.gmother > 0.05)) myfactors <- "1"
  cat(pname, " ", "\n")
  
  pvalues <- apply(numgeno, 1, function(numgeno) {
    numgenoDomm <- as.numeric(as.numeric(unlist(numgeno)) != 0)
	numgenoAddd <- as.numeric(unlist(numgeno))
    mdata <- data.frame(cbind(pheno = phenotypes[, pname], wg = phenotypes[, "WG"], grandmother = phenotypes[, "Grandma"], A = numgenoAddd, D = numgenoDomm))
    isNA <- which(apply(apply(mdata,1,is.na),2,any))
    if (length(isNA) > 0) mdata <- mdata[-isNA, ]
	myformula0 <- paste0("pheno ~ ", paste(myfactors, collapse = " + "))
	myformula <- paste0("pheno ~ ", paste(c(myfactors, "D", "A"), collapse = " + "))
	lm0 <- lm(formula = (myformula0), data = mdata)
    mmodel <- lm(formula(myformula), data = mdata)
    return(anova(mmodel, lm0)[["Pr(>F)"]][2])	
  })
  pmatrixADDDOM[names(pvalues), pname] <- pvalues
}
lodmatrixADDDOM <- -log10(pmatrixADDDOM)
write.table(lodmatrixADDDOM, file = "lodmatrixADDDOM_nosum.txt", quote = FALSE, sep = "\t")
lodannotmatrix <- cbind(annotation[rownames(lodmatrixADDDOM), ], lodmatrixADDDOM)

# Some plots
lines(lodmatrixDOM[,"D140"], main = "D140", col = as.numeric(as.factor(annotation[,"Chromosome"])), las = 2)
lines(lodmatrixADD[,"D140"], main = "D140", col = as.numeric(as.factor(annotation[,"Chromosome"])), las = 2)
plot(lodmatrixADDDOM[,"D140"], main = "D140", col = as.numeric(as.factor(annotation[,"chr"])), las = 2)
chr15 <- lodannotmatrix[which(lodannotmatrix[,"chr"] == 15),]
chr15ord <- chr15[order(chr15[,2], decreasing = FALSE),]
plot(main = "QTL profile bodyweight [Chr 15]", c(min(as.numeric(chr15[, "bp_mm10"])), max(as.numeric(chr15[, "bp_mm10"]))), c(0,9), ylab = "-log10 [pvalue]", xlab = "Position [mb]", las = 2, t = "n", xaxt = "n")
points(x = as.numeric(chr15ord[,"bp_mm10"]), y = chr15ord[,"D140"] , type = "l", col="dodgerblue", lwd = 1)