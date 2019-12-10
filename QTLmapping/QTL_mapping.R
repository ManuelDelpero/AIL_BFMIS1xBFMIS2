# QTL mapping
#
# copyright (c) 2018-2021 - Brockmann group - HU Berlin Manuel Delpero & Danny Arends
# 
# first written december, 2019

setwd("C:/Users/Manuel/Desktop/AIL_S1xS2/RAWDATA")

genotypes <- read.csv("genotypes.clean.txt", header = TRUE, sep = "\t", check.names = FALSE, colClasses = "character")
phenotypes <- read.csv("allPhenotypes.txt", header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
markerannot <- read.csv("snp_map.karl.txt", header = TRUE, sep = ",", check.names = FALSE, row.names = 1)

chromosomes <- c(1:19, "X", "Y")

annotation <- c()
for(chr in chromosomes){
  annotation <- rbind(annotation, markerannot[markerannot[,"chr"] == chr,])
}

# Make sure that the ordering between phenotypes and genotypes matches 
phenotypes <- phenotypes[-192,]
genotypes <- genotypes[,-192]
colnames(genotypes) <- gsub("AIL", "", colnames(genotypes))
phenotypes <- phenotypes[colnames(genotypes),]
genotypes <- genotypes[, rownames(phenotypes)]
write.table(genotypes, "OrderedGenotypes.txt", sep = "\t", quote=FALSE)

# Covariates we could/need to include in the model, we test them on their pvalue
wg <- as.factor(phenotypes[, "WG"])
mother <- phenotypes[, "Mutter"]
phenonames <- colnames(phenotypes)[-c(1,2, 58,59)]

pmatrix <- matrix(NA, nrow(genotypes), length(phenonames), dimnames= list(rownames(genotypes), phenonames))
for (pname in phenonames){
  pheno <- phenotypes[, pname]
  p.mother <- anova(lm(pheno ~ mother))["Pr(>F)"]["mother",]
  p.wg <- anova(lm(pheno ~ wg))["Pr(>F)"]["wg",]
  myfactors <- c()
  if(p.mother < 0.05) myfactors <- c(myfactors, "mother")
  if(p.wg < 0.05) myfactors <- c(myfactors, "wg")
  myformula <- paste0("pheno ~ ", paste(c(myfactors, "geno"), collapse = " + "))
  cat(pname, " ", myformula, "\n")
  
  pvalues <- apply(genotypes, 1, function(geno, pheno, wg, mother, myformula) {
    mmodel <- lm(formula(myformula))
    return(anova(mmodel)["Pr(>F)"]["geno",])
  }, pheno = phenotypes[,pname], wg = wg, mother = mother, myformula = myformula)
  pmatrix[names(pvalues), pname] <- pvalues
}
lodmatrix <- -log10(pmatrix)

write.table(lodmatrix, "lodmatrix.txt", sep = "\t", quote=FALSE)


# Manhattan plots
for (pname in  phenonames){   #phenonames) {
  plot(lodmatrix[,pname], main = pname, col = as.numeric(as.factor(annotation[,"chr"])))
  abline(h = -log10(0.05 / nrow(lodmatrix)), col="orange")
  abline(h = -log10(0.01 / nrow(lodmatrix)), col="green")
}

signmatrix <- lodmatrix[which(apply(lodmatrix, 1, function(x){ any(x > -log10(0.05 / nrow(lodmatrix))) })),]
signannotmatrix <- cbind(annotation[rownames(signmatrix), ], signmatrix)
write.table(signannotmatrix, "signannotmatrixxx.txt", sep = "\t", quote=FALSE)


