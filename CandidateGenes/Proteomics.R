# QC Anlysis on protemoics data from BGI

setwd("C:/Users/Manuel/Desktop/AIL_S1xS2/Proteomics/report/data")

data <- read.csv("Quantification/dia-protein-Summary.txt", sep = "\t", check.names = FALSE, header = TRUE, row.names = 1)
annotation <- read.csv("Identification/annotation_allprotein.txt", sep = "\t", check.names = FALSE, header = TRUE, row.names = 1)
geneExpr <- read.csv("C:/Users/Manuel/Desktop/AIL_S1xS2/RAWDATA/ExprGon.txt", sep = "\t", check.names = FALSE, header = TRUE)

# combine intensities with annotation and clean data
annot <- annotation[,c("Nr_Identity", "Swissprot Description")]
annot <- annot[grep("GN", annot[,2]),]
data <- data[rownames(annot),]
data <- data[,-31]
Genes <- strsplit(as.character(annot[,2])," ",)
Genes <- gsub("GN=", "", unlist(lapply(Genes, function(x) x[grepl("GN", x)])))
idx <- which(duplicated(Genes))
Genes <- Genes[-idx]
data <- data[-idx,]
annot <- annot[-idx,]

dataAnnot <- cbind(annot, data)

rownames(dataAnnot) <- Genes

GenesExpr <- as.character(geneExpr[,"mgi_symbol"])
geneExpr <- geneExpr[-which(duplicated(GenesExpr)),]
GenesExpr <- as.character(geneExpr[,"mgi_symbol"])
geneExpr <- geneExpr[-which(is.na(GenesExpr)),]

rownames(geneExpr) <- as.character(geneExpr[,"mgi_symbol"])

# some qc
corr <- cor(data, use = "pairwise.complete.obs")
heatmap(corr)

# Calculate correlation between protein and gene expression
geneExpr <- geneExpr[rownames(dataAnnot),]

geneExprS1 <- data.frame(meanS1 = geneExpr[,"mean(s1)"])
S1FatDataAnnot <- dataAnnot[,grep("Fat_S1", colnames(dataAnnot))]
S1FatProteins <- data.frame(Means=rowMeans(S1FatDataAnnot))
cor.test(S1FatProteins[,1], geneExprS1[,1])

geneExprS2 <- data.frame(meanS2 = geneExpr[,"mean(s2)"])
S2FatDataAnnot <- dataAnnot[,grep("Fat_S2", colnames(dataAnnot))]
S2FatProteins <- data.frame(Means=rowMeans(S2FatDataAnnot))
cor.test(S2FatProteins[,1], geneExprS2[,1])

# PCA with proteomics data
subsetPCAnum <- data[which(apply(apply(data,1,is.na),2,sum) == 0),]
pcares <- prcomp(t(subsetPCAnum))
summary(pcares)
rownames(pcares$x) <- colnames(subsetPCAnum)

selected <- rownames(pcares$x)[grep("S1", rownames(pcares$x))]
mcol <- as.numeric(rownames(pcares$x) %in% selected) + 1
plot(main = "PCA-Plot all phenotypes", pcares$x[,1], pcares$x[,2], xlab = "PC1 (0.56%)", ylab = "PC2 (0.23%)", col = c("blue", "red")[mcol], type='p', pch=20)
  legend("topright",
  legend = c("S1", "S2"),
  col = c("red", "blue"),
  pch = c(20,20,20),
  bty = "n",
  pt.cex = 2,
  cex = 1,
  text.col = "black",
  #horiz = F ,
  #inset = c(0.1, 0.1, 0.1)
)