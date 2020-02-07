# Plots to represent the gene expression results for AIL BFMI S1xS2
#
# copyright (c) 2018-2021 - Brockmann group - HU Berlin Manuel Delpero 
# 
# first written december, 2019

setwd("C:/Users/Manuel/Desktop/AIL_S1xS2/RAWDATA")

exprliver <- read.table("liverExpressions.txt", sep = "\t", header = TRUE, check.names = FALSE)
exprGonadalfat <- read.table("GonExpressions.txt", sep = "\t", header = TRUE, check.names = FALSE)
exprSkeletalmuscle <- read.table("SkExpressions.txt", sep = "\t", header = TRUE, check.names = FALSE)
exprPankreas <- read.table("PankreasExpressions.txt", sep = "\t", header = TRUE, check.names = FALSE)

Diffexprliver <- read.table("liver_significant_ann.txt", sep = "\t", header = TRUE, check.names = FALSE)
DiffexprGonadalfat <- read.table("gonadalfat_significant_ann.txt", sep = "\t", header = TRUE, check.names = FALSE)
DiffexprSkeletalmuscle <- read.table("skeletalmuscle_significant_ann.txt", sep = "\t", header = TRUE, check.names = FALSE)
DiffexprPankreas <- read.table("pankreas_significant_ann.txt", sep = "\t", header = TRUE, check.names = FALSE)

## Volcano Plots
# Gonadal fat
col <- as.numeric(abs(log2(exprGonadalfat[,"logFC"])) > 0.2630344) + as.numeric(exprGonadalfat[,"p.value"] < 0.10)+ as.numeric(exprGonadalfat[,"p.value"] < 0.10) 
plot(exprGonadalfat[, "logFC"], -log10(exprGonadalfat[, "p.value"]))
abline(v=0.2630344)
abline(v=-0.2630344)
abline(h=-log10(0.05/nrow(exprGonadalfat)))


plot(exprliver[, "logFC"], -log10(exprliver[, "p.value"]))
plot(skeletalmuscle[, "logFC"], -log10(skeletalmuscle[, "p.value"]))
plot(pankreas[, "logFC"], -log10(pankreas[, "p.value"]))

