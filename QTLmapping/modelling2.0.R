# Modelling using AIL BFMI S1xS2
#
# copyright (c) 2018-2021 - Brockmann group - HU Berlin Manuel Delpero & Danny Arends
# 
# first written december, 2019

setwd("C:/Users/Manuel/Desktop/AIL_S1xS2/RAWDATA")

genotypes <- read.csv("genotypes.cleaned.txt", header = TRUE, check.names = FALSE, sep="\t", colClasses="character")
phenotypes <- read.csv("PhenotypesComplete.txt", header = TRUE, check.names = FALSE, sep="\t", row.names=1)
annotation <- read.csv("map.cleaned.txt", header=TRUE, sep="\t", check.names=FALSE)
annotation <- annotation[, c(1,2,3,6)]
colnames(annotation) <- c("Chromosome", "Position", "GenTrain Score", "SNP")
colnames(genotypes) <- gsub("AIL", "" , colnames(genotypes)) 
dim(genotypes)
dim(phenotypes)

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

phenonames <- colnames(phenotypes[,c(3:57, 61:72)])


#phenonames <- colnames(phenotypes[,c(3:57, 61, 62, 65:75)])
#write.table(phenotypes, file = "PhenotypesComplete.txt", sep = "\t", quote = FALSE, row.names = TRUE)
# Covariates we could/need to include in the model
phenotypes <- phenotypes[colnames(genotypes),]
wg <- phenotypes[, "WG"]
grandmother <- phenotypes[, "Grandma"]
Leber <- phenotypes[, "Leber"]
Gon <- phenotypes[, "Gon"]
phenonames <- c("Gon", "Leber", "Gluc172", "D174", "Triglycerides", "ITTauc")

## Choose the best model for each phenotype
# map first using the raw models for each phenotype
phenonames <- c("Gon", "Leber", "Gluc172", "D174", "Triglycerides", "ITTauc")

# Dom dev + Add model (real dom)
pmatrixADDDOM <- matrix(NA, nrow(numgeno), length(phenonames), dimnames= list(rownames(numgeno), phenonames))
for (pname in phenonames){
  cat(pname, " ", "\n")
  pvalues <- apply(numgeno, 1, function(numgeno) {
    numgenoDomm <- as.numeric(as.numeric(unlist(numgeno)) != 0)
    numgenoAddd <- as.numeric(unlist(numgeno))
    mdata <- data.frame(cbind(pheno = phenotypes[, pname], sex = phenotypes[, "Sex"], liver = phenotypes[, "Leber"], gon = phenotypes[, "Gon"], wg = phenotypes[,"WG"], A = numgenoAddd, D = numgenoDomm))
    isNA <- which(apply(apply(mdata,1,is.na),2,any))
    if (length(isNA) > 0) mdata <- mdata[-isNA, ]
	lm0 <- lm(pheno ~ 1 + sex, data = mdata)
    mmodel <- lm(pheno ~ D + A + sex, data = mdata)
    return(anova(mmodel, lm0)[["Pr(>F)"]][2])	
  })
  pmatrixADDDOM[names(pvalues), pname] <- pvalues
}
lodmatrixADDDOMrow <- -log10(pmatrixADDDOM)
#lodannotmatrix <- cbind(annotation[rownames(lodmatrixADDDOM), ], lodmatrixADDDOM)

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
lodmatrixDOMrow <- -log10(pmatrixDOM)

# Additive model
pmatrixADD <- matrix(NA, nrow(numgeno), length(phenonames), dimnames= list(rownames(numgeno), phenonames))
for (pname in phenonames){
  cat(pname, " ", "\n")
  pvalues <- apply(numgeno, 1, function(numgeno) {
    numgenoAddd <- as.numeric(unlist(numgeno))
    mdata <- data.frame(cbind(pheno = phenotypes[, pname], sex = phenotypes[, "Sex"], liver = phenotypes[, "Leber"], gon = phenotypes[, "Gon"], wg = phenotypes[,"WG"], A = numgenoAddd))
    isNA <- which(apply(apply(mdata,1,is.na),2,any))
    if (length(isNA) > 0) mdata <- mdata[-isNA,]
	lm0 <- lm(pheno ~ 1 + sex, data = mdata)
    mmodel <- lm(pheno ~ A + sex, data = mdata)	
    return(anova(mmodel, lm0)[["Pr(>F)"]][2])	
  })
  pmatrixADD[names(pvalues), pname] <- pvalues
}
lodmatrixADDrow <- -log10(pmatrixADD)

## start defining the right models for each phenotype
## model for Blood glucose 
resLiver <- anova(lm(phenotypes[,"Gluc172"] ~ phenotypes[,"Leber"]))
resGon <- anova(lm(phenotypes[,"Gluc172"] ~ phenotypes[,"Gon"]))
phenotype <- as.numeric(phenotypes[, "Gluc172"])

# Define glucose variance explained by liver
varExplainedLiver  <- resLiver[, "Sum Sq"] / sum((phenotype - mean(phenotype, na.rm=TRUE))^2, na.rm=TRUE)
round(varExplainedLiver * 100, digits=1) 

# Define glucose variance explained by Gonadal fat
varExplainedGon  <- resGon[, "Sum Sq"] / sum((phenotype - mean(phenotype, na.rm=TRUE))^2, na.rm=TRUE)
round(varExplainedGon * 100, digits=1) 

# Colinearity between gonadal fat an liver therefore is difficult to define the exact variance explained of glucose 
# by each phenotype but only the range (up to 60.7 % for liver and up to 39.2 % for gonadal fat) 
# liver and gonadal fat explains most of the variation of glucose
# Using the raw models we identified overlapping QTLs for gonadal fat, liver and glucose and also liver triglycerides 
# we assume that the QTLs for glucose that overlaps with the other QTLs for liver and gonadal fat have an indirect effect on glucose
# we know that liver and gonadal fat can influence the blood glucose
# Final decision: Map glucose including liver and gonadal fat as covariates to identify QTLs directly responsible for the glucose levels and not indirectly

## Models for gonadal fat weight and liver
# See if the liver influences the gonadal fat and the other way round with causal modelling 
anova(lm(phenotypes[, "Gon"] ~ phenotypes[,"Leber"]))
anova(lm(phenotypes[, "Leber"] ~ phenotypes[,"Gon"]))
hist(residuals(lm(phenotypes[, "Gon"] ~ 1)))
hist(residuals(lm(phenotypes[, "Gon"] ~ phenotypes[,"Leber"])))
hist(residuals(lm(phenotypes[, "Leber"] ~ 1)))

numgenoDomm <- as.numeric(as.numeric(unlist(numgeno["UNCHS043909",])) != 0)
numgenoAddd <- as.numeric(unlist(numgeno["UNCHS043909",]))
mdata <- data.frame(cbind(sex = phenotypes[, "Sex"], liver = phenotypes[, "Leber"], gon = phenotypes[, "Gon"], wg = phenotypes[,"WG"], A = numgenoAddd, D = numgenoDomm))
isNA <- which(apply(apply(mdata,1,is.na),2,any))
if (length(isNA) > 0) mdata <- mdata[-isNA, ]

lm1 <- lm(gon ~ sex + liver + A + D, data = mdata)
lm0 <- lm(gon ~ sex + liver, data = mdata)
anova(lm1,lm0)
# When we add the liver as a covariate the LOD score of the QTL on chr17 for gonadal fat drops from 8.59 to 3.65 

lm1 <- lm(liver ~ sex + gon + A + D, data = mdata)
lm0 <- lm(liver ~ sex + gon, data = mdata)
anova(lm1,lm0)
# The same happen for liver, adding gonadal fat as a covariate the Lod score drops from 7.5 to 2.9

# They both influence each other
# Final decision: map gonadal fat and liver using the raw models since we are not able do define exactly which one is influecning the other one

# Liver triglycerides
numgenoDomm <- as.numeric(as.numeric(unlist(numgeno["JAX00632487",])) != 0)
numgenoAddd <- as.numeric(unlist(numgeno["JAX00632487",]))
mdata <- data.frame(cbind(sex = phenotypes[, "Sex"], liver = phenotypes[, "Leber"], trig = phenotypes[, "Triglycerides"], wg = phenotypes[,"WG"], A = numgenoAddd, D = numgenoDomm))
isNA <- which(apply(apply(mdata,1,is.na),2,any))
if (length(isNA) > 0) mdata <- mdata[-isNA, ]

lm1 <- lm(trig ~ sex + liver + A , data = mdata)
lm0 <- lm(trig ~ sex + liver, data = mdata)
anova(lm1,lm0)
# Final decision: correct liver triglycerides using liver weight

pmatrixADDDOM <- matrix(NA, nrow(numgeno), length(phenonames), dimnames= list(rownames(numgeno), phenonames))
for (pname in phenonames){
  cat(pname, " ", "\n")
  pvalues <- apply(numgeno, 1, function(numgeno) {
    numgenoDomm <- as.numeric(as.numeric(unlist(numgeno)) != 0)
    numgenoAddd <- as.numeric(unlist(numgeno))
    mdata <- data.frame(cbind(pheno = phenotypes[, pname], sex = phenotypes[, "Sex"], liver = phenotypes[, "Leber"], gon = phenotypes[, "Gon"], wg = phenotypes[,"WG"], A = numgenoAddd, D = numgenoDomm))
    isNA <- which(apply(apply(mdata,1,is.na),2,any))
    if (length(isNA) > 0) mdata <- mdata[-isNA, ]
	if (pname == "Gluc172"){ 
      lm0 <- lm(pheno ~ 1 + sex + liver + gon, data = mdata)
      mmodel <- lm(pheno ~ D + A + sex + gon + liver, data = mdata)	
	} else if (pname == "Triglycerides"){
	  lm0 <- lm(pheno ~ 1 + sex + liver, data = mdata)
      mmodel <- lm(pheno ~ D + A + sex + liver, data = mdata)	
    } else {
	  lm0 <- lm(pheno ~ 1 + sex, data = mdata)
      mmodel <- lm(pheno ~ D + A + sex, data = mdata)
    }
    return(anova(mmodel, lm0)[["Pr(>F)"]][2])	
  })
  pmatrixADDDOM[names(pvalues), pname] <- pvalues
}
lodmatrixADDDOM <- -log10(pmatrixADDDOM)

# Dominance dev model
pmatrixDOM <- matrix(NA, nrow(numgeno), length(phenonames), dimnames= list(rownames(numgeno), phenonames))
for (pname in phenonames){
  cat(pname, " ", "\n")
  pvalues <- apply(numgeno, 1, function(numgeno) {
    numgenoDomm <- as.numeric(as.numeric(unlist(numgeno)) != 0)
    mdata <- data.frame(cbind(pheno = phenotypes[, pname], liver = phenotypes[,"Leber"], gon = phenotypes[,"Gon"], sex = phenotypes[, "Sex"], D = numgenoDomm))
    isNA <- which(apply(apply(mdata,1,is.na),2,any))
    if (length(isNA) > 0) mdata <- mdata[-isNA, ]
    if (pname == "Gluc172"){ 
      lm0 <- lm(pheno ~ 1 + sex + liver + gon, data = mdata)
      mmodel <- lm(pheno ~ D + sex + liver + gon, data = mdata)	
	} else if (pname == "Triglycerides"){
	  lm0 <- lm(pheno ~ 1 + sex + liver, data = mdata)
      mmodel <- lm(pheno ~ D + sex + liver, data = mdata)	
    } else {
	  lm0 <- lm(pheno ~ 1 + sex, data = mdata)
      mmodel <- lm(pheno ~ D + sex, data = mdata)
    }
    return(anova(mmodel, lm0)[["Pr(>F)"]][2])	
  })
  pmatrixDOM[names(pvalues), pname] <- pvalues
}
lodmatrixDOM <- -log10(pmatrixDOM)

# Additive model
pmatrixADD <- matrix(NA, nrow(numgeno), length(phenonames), dimnames= list(rownames(numgeno), phenonames))
for (pname in phenonames){
  cat(pname, " ", "\n")
  pvalues <- apply(numgeno, 1, function(numgeno) {
    numgenoAddd <- as.numeric(unlist(numgeno))
    mdata <- data.frame(cbind(pheno = phenotypes[, pname], sex = phenotypes[, "Sex"], liver = phenotypes[, "Leber"], gon = phenotypes[, "Gon"], wg = phenotypes[,"WG"], A = numgenoAddd))
    isNA <- which(apply(apply(mdata,1,is.na),2,any))
    if (length(isNA) > 0) mdata <- mdata[-isNA,]
	if (pname == "Gluc172"){ 
      lm0 <- lm(pheno ~ 1 + sex + liver + gon, data = mdata)
      mmodel <- lm(pheno ~ A + sex + liver + gon, data = mdata)
	} else if (pname == "Triglycerides"){
	  lm0 <- lm(pheno ~ 1 + sex + liver, data = mdata)
      mmodel <- lm(pheno ~ A + sex + liver, data = mdata)		  
    } else {
	  lm0 <- lm(pheno ~ 1 + sex, data = mdata)
      mmodel <- lm(pheno ~ A + sex, data = mdata)
    }	
    return(anova(mmodel, lm0)[["Pr(>F)"]][2])	
  })
  pmatrixADD[names(pvalues), pname] <- pvalues
}
lodmatrixADD <- -log10(pmatrixADD)

## make plots with row models and adjusted models

# make some transparancy colors, note: always pass alpha on the 0-255 scale

makeTransparent<-function(someColor, alpha=100)
{
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
    blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}
makeTransparent("red")

# Plots for row models 
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

par(cex.lab=1.2, cex.main = 1.3, cex.axis = 1)
mat <- matrix(c(1,1,2,2), 2, 2, byrow = TRUE)
layout(mat, widths = rep.int(3, ncol(mat)))

chrs <- c(1:19,"X")
gap <- 80000000
map.sorted <- NULL
chr.lengths <- c()
chr.starts <- c(0)
chrmids <- c()
i <- 1
for(chr in chrs){
  onChr <- which(markerannot[,"Chromosome"] == chr)
  map.sorted <- rbind(map.sorted, markerannot[onChr,])
  chr.lengths <- c(chr.lengths, max(markerannot[onChr, "Position"]))
  chr.starts <- c(chr.starts, chr.starts[i] + max(markerannot[onChr, "Position"]) + gap)
  i <- i + 1
}

chr.start <- chr.starts[-5]
chr.ends <- chr.start + chr.lengths
names(chr.starts) <- chrs
names(chr.lengths) <- chrs

for (x in chrs){
  chrmid <- as.numeric(chr.lengths[x]/2) + as.numeric(chr.starts[x])
  chrmids <- c(chrmids, chrmid)
}

	
plot(x = c(-gap, tail(chr.starts,1)), y = c(0,9), t = 'n', xlab="Chromosome", ylab="-log10[P]",xaxt='n', xaxs="i", yaxs="i", las=2, main=paste0("Manhattan plots - Raw models"))
phenotype <- "Gon"
for(chr in chrs){
  onChr <- rownames(map.sorted[map.sorted[,"Chromosome"] == chr,])
  currentADDDOM <- lodmatrixADDDOMrow[onChr, phenotype]
  currentDOM <- lodmatrixDOMrow[onChr, phenotype]
  currentADD <- lodmatrixADDrow[onChr, phenotype]
  if (chr == "X")
    points(x=chr.starts[chr] + map.sorted[onChr,"Position"], y = lodmatrixADDrow[onChr, phenotype], t ='p', pch = 16, cex = 1.5, col= "cornflowerblue")
  for (p in 1:length(currentADDDOM)){
    pos <- chr.starts[chr] + map.sorted[onChr,"Position"]
    if ((currentADDDOM[p] >  currentDOM[p]) && (currentADDDOM[p] > currentADD[p])){
      if (chr %in% seq(1,20,2))
        points(x=pos[p], y = currentADDDOM[p], t ='p', pch = 18, cex = 1.2, col= "cornflowerblue")
      else (points(x=pos[p], y = currentADDDOM[p], t ='p', pch = 18, cex = 1.2, col= "cornflowerblue"))
    }
    if ((currentDOM[p] >  currentADDDOM[p]) && (currentDOM[p] > currentADD[p])){
      if (chr %in% seq(1,20,2))
        points(x=pos[p], y = currentDOM[p], t ='p', pch = 18, cex = 1,2, col= "cornflowerblue")
      else (points(x=pos[p], y = currentDOM[p], t ='p', pch = 18, cex = 1.2, col= "cornflowerblue"))
    }
      if ((currentADD[p] >  currentADDDOM[p]) && (currentADD[p] > currentDOM[p])){
        if (chr %in% seq(1,20,2))
          points(x=pos[p], y = currentADD[p], t ='p', pch = 18, cex = 1.2, col= "cornflowerblue")
	else (points(x=pos[p], y = currentADD[p], t ='p', pch = 18, cex = 1.2,col= "cornflowerblue"))
    }
  }
}
phenotype <- "Leber"
for(chr in chrs){
  onChr <- rownames(map.sorted[map.sorted[,"Chromosome"] == chr,])
  currentADDDOM <- lodmatrixADDDOMrow[onChr, phenotype]
  currentDOM <- lodmatrixDOMrow[onChr, phenotype]
  currentADD <- lodmatrixADDrow[onChr, phenotype]
  if (chr == "X")
    points(x=chr.starts[chr] + map.sorted[onChr,"Position"], y = lodmatrixADDrow[onChr, phenotype], t ='p', pch = 16, cex = 1.5, col= "lightgreen")
  for (p in 1:length(currentADDDOM)){
    pos <- chr.starts[chr] + map.sorted[onChr,"Position"]
    if ((currentADDDOM[p] >  currentDOM[p]) && (currentADDDOM[p] > currentADD[p])){
      if (chr %in% seq(1,20,2))
        points(x=pos[p], y = currentADDDOM[p], t ='p', pch = 18, cex = 1.2, col= "lightgreen")
      else (points(x=pos[p], y = currentADDDOM[p], t ='p', pch = 18, cex = 1.2, col= "lightgreen"))
    }
    if ((currentDOM[p] >  currentADDDOM[p]) && (currentDOM[p] > currentADD[p])){
      if (chr %in% seq(1,20,2))
        points(x=pos[p], y = currentDOM[p], t ='p', pch = 18, cex = 1.2, col= "lightgreen")
      else (points(x=pos[p], y = currentDOM[p], t ='p', pch = 18, cex = 1.2, col= "lightgreen"))
    }
      if ((currentADD[p] >  currentADDDOM[p]) && (currentADD[p] > currentDOM[p])){
        if (chr %in% seq(1,20,2))
          points(x=pos[p], y = currentADD[p], t ='p', pch = 18, cex = 1.2, col= "lightgreen")
	else (points(x=pos[p], y = currentADD[p], t ='p', pch = 18, cex = 1.2,col= "lightgreen"))
    }
  }
}
phenotype <- "Gluc172"
for(chr in chrs){
  onChr <- rownames(map.sorted[map.sorted[,"Chromosome"] == chr,])
  currentADDDOM <- lodmatrixADDDOMrow[onChr, phenotype]
  currentDOM <- lodmatrixDOMrow[onChr, phenotype]
  currentADD <- lodmatrixADDrow[onChr, phenotype]
  if (chr == "X")
    lines(x=chr.starts[chr] + map.sorted[onChr,"Position"], y = lodmatrixADDrow[onChr, phenotype], t ='p', pch = 16, cex = 1.5, col = rgb(0, 0, 255, max = 255, alpha = 125, names = "blue50"))
  for (p in 1:length(currentADDDOM)){
    pos <- chr.starts[chr] + map.sorted[onChr,"Position"]
    if ((currentADDDOM[p] >  currentDOM[p]) && (currentADDDOM[p] > currentADD[p])){
      if (chr %in% seq(1,20,2))
        lines(x=pos[p], y = currentADDDOM[p], t ='p', pch = 18, cex = 1.2, col= makeTransparent("red"))
      else (lines(x=pos[p], y = currentADDDOM[p], t ='p', pch = 18, cex = 1.2, col= makeTransparent("red")))
    }
    if ((currentDOM[p] >  currentADDDOM[p]) && (currentDOM[p] > currentADD[p])){
      if (chr %in% seq(1,20,2))
        lines(x=pos[p], y = currentDOM[p], t ='p', pch = 18, cex = 1.2, col= makeTransparent("red"))
      else (lines(x=pos[p], y = currentDOM[p], t ='p', pch = 18, cex = 1.2, col= makeTransparent("red")))
    }
      if ((currentADD[p] >  currentADDDOM[p]) && (currentADD[p] > currentDOM[p])){
        if (chr %in% seq(1,20,2))
          lines(x=pos[p], y = currentADD[p], t ='p', pch = 18, cex = 1.2, col= makeTransparent("red"))
	else (lines(x=pos[p], y = currentADD[p], t ='p', pch = 18, cex = 1.2,col= makeTransparent("red")))
    }
  }
}
#rect(xc[1:4], c(0,0,0,0), xc[5:8], yc[1:4], border = "red")
axis(1, chrs, at = chrmids)
abline(h = 4.3, col="orange",lty=3)
#abline(v = c(chr.start,chr.ends), col = "red") 
abline(h= 4.7, col="green",lty=3)
axis(1, chrs, at = chrmids)
legend("topright", #bg="gray"
  bty = "n",
  legend = c("Gonadal fat", "Liver", "Blood glucose"),
  pch = c(15, 15),
  col = c("cornflowerblue", "lightgreen", "red"))
  
  
# Plots for full models
chrs <- c(1:19,"X")
gap <- 80000000
map.sorted <- NULL
chr.lengths <- c()
chr.starts <- c(0)
chrmids <- c()
i <- 1
for(chr in chrs){
  onChr <- which(markerannot[,"Chromosome"] == chr)
  map.sorted <- rbind(map.sorted, markerannot[onChr,])
  chr.lengths <- c(chr.lengths, max(markerannot[onChr, "Position"]))
  chr.starts <- c(chr.starts, chr.starts[i] + max(markerannot[onChr, "Position"]) + gap)
  i <- i + 1
}

chr.start <- chr.starts[-5]
chr.ends <- chr.start + chr.lengths
names(chr.starts) <- chrs
names(chr.lengths) <- chrs

for (x in chrs){
  chrmid <- as.numeric(chr.lengths[x]/2) + as.numeric(chr.starts[x])
  chrmids <- c(chrmids, chrmid)
}
	
lodmatrixDOM <- lodmatrixDOM[rownames(annotation),]
lodmatrixADD <- lodmatrixADD[rownames(annotation),]
lodmatrixADDDOM <- lodmatrixADDDOM[rownames(annotation),]

plot(x = c(-gap, tail(chr.starts,1)), y = c(0,9), t = 'n', xlab="Chromosome", ylab="-log10[P]",xaxt='n', xaxs="i", yaxs="i", las=2, main=paste0("Manhattan plots - Full model blood glucose"))
phenotype <- "Gon"
for(chr in chrs){
  onChr <- rownames(map.sorted[map.sorted[,"Chromosome"] == chr,])
  currentADDDOM <- lodmatrixADDDOM[onChr, phenotype]
  currentDOM <- lodmatrixDOM[onChr, phenotype]
  currentADD <- lodmatrixADD[onChr, phenotype]
  if (chr == "X")
    points(x=chr.starts[chr] + map.sorted[onChr,"Position"], y = lodmatrixADDrow[onChr, phenotype], t ='p', pch = 16, cex = 1.5, col= "cornflowerblue")
  for (p in 1:length(currentADDDOM)){
    pos <- chr.starts[chr] + map.sorted[onChr,"Position"]
    if ((currentADDDOM[p] >  currentDOM[p]) && (currentADDDOM[p] > currentADD[p])){
      if (chr %in% seq(1,20,2))
        points(x=pos[p], y = currentADDDOM[p], t ='p', pch = 18, cex = 1.2, col= "cornflowerblue")
      else (points(x=pos[p], y = currentADDDOM[p], t ='p', pch = 18, cex = 1.2, col= "cornflowerblue"))
    }
    if ((currentDOM[p] >  currentADDDOM[p]) && (currentDOM[p] > currentADD[p])){
      if (chr %in% seq(1,20,2))
        points(x=pos[p], y = currentDOM[p], t ='p', pch = 18, cex = 1,2, col= "cornflowerblue")
      else (points(x=pos[p], y = currentDOM[p], t ='p', pch = 18, cex = 1.2, col= "cornflowerblue"))
    }
      if ((currentADD[p] >  currentADDDOM[p]) && (currentADD[p] > currentDOM[p])){
        if (chr %in% seq(1,20,2))
          points(x=pos[p], y = currentADD[p], t ='p', pch = 18, cex = 1.2, col= "cornflowerblue")
	else (points(x=pos[p], y = currentADD[p], t ='p', pch = 18, cex = 1.2,col= "cornflowerblue"))
    }
  }
}
phenotype <- "Leber"
for(chr in chrs){
  onChr <- rownames(map.sorted[map.sorted[,"Chromosome"] == chr,])
  currentADDDOM <- lodmatrixADDDOM[onChr, phenotype]
  currentDOM <- lodmatrixDOM[onChr, phenotype]
  currentADD <- lodmatrixADD[onChr, phenotype]
  if (chr == "X")
    points(x=chr.starts[chr] + map.sorted[onChr,"Position"], y = lodmatrixADDrow[onChr, phenotype], t ='p', pch = 16, cex = 1.5, col= "lightgreen")
  for (p in 1:length(currentADDDOM)){
    pos <- chr.starts[chr] + map.sorted[onChr,"Position"]
    if ((currentADDDOM[p] >  currentDOM[p]) && (currentADDDOM[p] > currentADD[p])){
      if (chr %in% seq(1,20,2))
        points(x=pos[p], y = currentADDDOM[p], t ='p', pch = 18, cex = 1.2, col= "lightgreen")
      else (points(x=pos[p], y = currentADDDOM[p], t ='p', pch = 18, cex = 1.2, col= "lightgreen"))
    }
    if ((currentDOM[p] >  currentADDDOM[p]) && (currentDOM[p] > currentADD[p])){
      if (chr %in% seq(1,20,2))
        points(x=pos[p], y = currentDOM[p], t ='p', pch = 18, cex = 1.2, col= "lightgreen")
      else (points(x=pos[p], y = currentDOM[p], t ='p', pch = 18, cex = 1.2, col= "lightgreen"))
    }
      if ((currentADD[p] >  currentADDDOM[p]) && (currentADD[p] > currentDOM[p])){
        if (chr %in% seq(1,20,2))
          points(x=pos[p], y = currentADD[p], t ='p', pch = 18, cex = 1.2, col= "lightgreen")
	else (points(x=pos[p], y = currentADD[p], t ='p', pch = 18, cex = 1.2,col= "lightgreen"))
    }
  }
}
phenotype <- "Gluc172"
for(chr in chrs){
  onChr <- rownames(map.sorted[map.sorted[,"Chromosome"] == chr,])
  currentADDDOM <- lodmatrixADDDOM[onChr, phenotype]
  currentDOM <- lodmatrixDOM[onChr, phenotype]
  currentADD <- lodmatrixADD[onChr, phenotype]
  if (chr == "X")
    points(x=chr.starts[chr] + map.sorted[onChr,"Position"], y = lodmatrixADDrow[onChr, phenotype], t ='p', pch = 16, cex = 1.5, col= "red")
  for (p in 1:length(currentADDDOM)){
    pos <- chr.starts[chr] + map.sorted[onChr,"Position"]
    if ((currentADDDOM[p] >  currentDOM[p]) && (currentADDDOM[p] > currentADD[p])){
      if (chr %in% seq(1,20,2))
        points(x=pos[p], y = currentADDDOM[p], t ='p', pch = 18, cex = 1.2, col = makeTransparent("red"))
      else (points(x=pos[p], y = currentADDDOM[p], t ='p', pch = 18, cex = 1.2, col = makeTransparent("red")))
    }
    if ((currentDOM[p] >  currentADDDOM[p]) && (currentDOM[p] > currentADD[p])){
      if (chr %in% seq(1,20,2))
        points(x=pos[p], y = currentDOM[p], t ='p', pch = 18, cex = 1.2, col = makeTransparent("red"))
      else (points(x=pos[p], y = currentDOM[p], t ='p', pch = 18, cex = 1.2, col = makeTransparent("red")))
    }
      if ((currentADD[p] >  currentADDDOM[p]) && (currentADD[p] > currentDOM[p])){
        if (chr %in% seq(1,20,2))
          points(x=pos[p], y = currentADD[p], t ='p', pch = 18, cex = 1.2, col = makeTransparent("red"))
	else (points(x=pos[p], y = currentADD[p], t ='p', pch = 18, cex = 1.2, col = makeTransparent("red")))
    }
  }
}
#rect(xc[1:4], c(0,0,0,0), xc[5:8], yc[1:4], border = "red")
axis(1, chrs, at = chrmids)
abline(h = 4.3, col="orange",lty=3)
#abline(v = c(chr.start,chr.ends), col = "red") 
abline(h= 4.7, col="green",lty=3)
axis(1, chrs, at = chrmids)
legend("topright", #bg="gray"
  bty = "n",
  legend = c("Gonadal fat", "Liver", "Blood glucose"),
  pch = c(15, 15),
  col = c("cornflowerblue", "lightgreen", "red"))

## Pie charts with variance explained chr 3
# Glucose 
numgenoDomm <- as.numeric(as.numeric(unlist(numgeno["JAX00063853",])) != 0)
numgenoAddd <- as.numeric(unlist(numgeno["JAX00063853",]))
mdata <- data.frame(cbind(glucose = phenotypes[, "Gluc172"], sex = phenotypes[, "Sex"], liver = phenotypes[, "Leber"], gon = phenotypes[, "Gon"], wg = phenotypes[,"WG"], A = numgenoAddd, D = numgenoDomm))
isNA <- which(apply(apply(mdata,1,is.na),2,any))
if (length(isNA) > 0) mdata <- mdata[-isNA, ]
lm1 <- lm(gon ~ sex + liver + A + D, data = mdata)
lm0 <- lm(gon ~ sex + liver, data = mdata)
anova(lm1,lm0)
res <- anova(lm(glucose ~ sex + liver + gon + A + D, data = mdata))
phenotype <- as.numeric(mdata[,"glucose"])
varExplained <- res[,"Sum Sq"] / sum((phenotype - mean(phenotype, na.rm = TRUE))^2, na.rm = TRUE)
varExplained <- round(varExplained * 100, digits =1)
varExplained <- c(20.8, 1.9, 20.7, 21.2, 33.7, 1.6)
slices <- c(varExplained)
lbls <- c("Sex (20.8 %)","Gonadal fat (1.9 %)", "Shared gonadal fat and liver (20.7 %)", "Liver (21.2 %)","Unknown (33.7 %)", "QTL chr 15 (1.6 %)")
pie(slices, labels = lbls, main="Blood glucose variance explained", init.angle = 120, col = c("lightblue", "darkorange", "darkorange2", "darkorange4", "white", "darkolivegreen3"))

# Liver
numgenoDomm <- as.numeric(as.numeric(unlist(numgeno["UNCHS043909",])) != 0)
numgenoAddd <- as.numeric(unlist(numgeno["UNCHS043909",]))
mdata <- data.frame(cbind(glucose = phenotypes[, "Gluc172"], sex = phenotypes[, "Sex"], liver = phenotypes[, "Leber"], gon = phenotypes[, "Gon"], wg = phenotypes[,"WG"], A = numgenoAddd, D = numgenoDomm))
isNA <- which(apply(apply(mdata,1,is.na),2,any))
if (length(isNA) > 0) mdata <- mdata[-isNA, ]
lm1 <- lm(liver ~ sex + A + D, data = mdata)
lm0 <- lm(liver ~ sex, data = mdata)
anova(lm1,lm0)
res <- anova(lm(liver ~ sex + A + D, data = mdata))
phenotype <- as.numeric(mdata[,"liver"])
varExplained <- res[,"Sum Sq"] / sum((phenotype - mean(phenotype, na.rm = TRUE))^2, na.rm = TRUE)
varExplained <- round(varExplained * 100, digits =1)
varExplained <- c(31.6, 4.2, 64.3)
slices <- c(varExplained)
lbls <- c("Sex (31.6 %)", "QTL chr 17 (4.2 %)", "Unknown (64.3 %)")
pie(slices, labels = lbls, main="Liver variance explained", col = c("lightblue", "orange", "white"))

# Gonadal fat QTL chr 17
numgenoDomm <- as.numeric(as.numeric(unlist(numgeno["UNCHS043909",])) != 0)
numgenoAddd <- as.numeric(unlist(numgeno["UNCHS043909",]))
mdata <- data.frame(cbind(glucose = phenotypes[, "Gluc172"], sex = phenotypes[, "Sex"], liver = phenotypes[, "Leber"], gon = phenotypes[, "Gon"], wg = phenotypes[,"WG"], A = numgenoAddd, D = numgenoDomm))
isNA <- which(apply(apply(mdata,1,is.na),2,any))
if (length(isNA) > 0) mdata <- mdata[-isNA, ]
lm1 <- lm(gon ~ sex + A + D, data = mdata)
lm0 <- lm(gon ~ sex, data = mdata)
anova(lm1,lm0)
res <- anova(lm(gon ~ sex + A + D, data = mdata))
phenotype <- as.numeric(mdata[,"gon"])
varExplained <- res[,"Sum Sq"] / sum((phenotype - mean(phenotype, na.rm = TRUE))^2, na.rm = TRUE)
varExplained <- round(varExplained * 100, digits =1)
varExplained <- c(73.1, 2.1, 24.9)
slices <- c(varExplained)
lbls <- c("Sex (73.1 %)", "QTL chr 17 (2.1 %)", "Unknown (24.9 %)")
pie(slices, labels = lbls, main="Gonadal fat variance explained", init.angle = 180, col = c("lightblue", "orange", "white"))

# Gonadal fat QTL chr 3 
numgenoDomm <- as.numeric(as.numeric(unlist(numgeno["UNC5791802",])) != 0)
numgenoAddd <- as.numeric(unlist(numgeno["UNC5791802",]))
mdata <- data.frame(cbind(glucose = phenotypes[, "Gluc172"], sex = phenotypes[, "Sex"], liver = phenotypes[, "Leber"], gon = phenotypes[, "Gon"], wg = phenotypes[,"WG"], A = numgenoAddd, D = numgenoDomm))
isNA <- which(apply(apply(mdata,1,is.na),2,any))
if (length(isNA) > 0) mdata <- mdata[-isNA, ]
lm1 <- lm(gon ~ sex + A + D, data = mdata)
lm0 <- lm(gon ~ sex, data = mdata)
anova(lm1,lm0)
res <- anova(lm(gon ~ sex + A , data = mdata))
phenotype <- as.numeric(mdata[,"gon"])
varExplained <- res[,"Sum Sq"] / sum((phenotype - mean(phenotype, na.rm = TRUE))^2, na.rm = TRUE)
varExplained <- round(varExplained * 100, digits =1)
slices <- c(varExplained)
lbls <- c("QTL chr 3 (14.4 %)", "Unknown (85.6 %)")
pie(slices, labels = lbls, main="Gonadal fat variance explained", col = c("orange", "white"))
