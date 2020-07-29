library(ctl)

setwd("C:/Users/Manuel/Desktop/AIL_S1xS2/RAWDATA")

pheno <- read.table("allPhenotypes_final.txt", sep = "\t", row.names=1)
geno <- read.table("genotypes.cleaned.txt", sep = "\t", colClasses="character")
map <- read.table("snp_map.karl.txt", sep = ",", colClasses="character", header= TRUE, row.names=1)
map <- map[which(rownames(map) %in% rownames(geno)),]

geno <- geno[rownames(map),]

rownames(pheno) <- paste0("AIL", rownames(pheno))

pheno <- pheno[colnames(geno),]
geno <- t(geno)
genot <- apply(geno, 2, function(x){as.numeric(factor(x, levels=c("A", "H", "B")))})
rownames(genot) <- rownames(geno)


pheno <- pheno[,-which(colnames(pheno) %in% c("Sex", "ID", "WG", "Mutter", "Grandma"))]
pheno <- pheno[, c("Gon","Leber")]
phenot <- apply(pheno,2,as.numeric)

res <- CTLscan(genotypes=genot, phenotypes = phenot, adjust = FALSE)  # chr15 there are 2 CTLs

# plot to represent the ctl curve across chr15
sig <- map[names(which(res$Leber$ctl[, "Gon"] > 4)),]
chr15 <- rownames(map[which(map[,"chr"] == 15),])
lods <- res$Leber$ctl[, "Gon"]
lods <- data.frame(lods)
lods15 <- lods[chr15,]
lines(lods15)


