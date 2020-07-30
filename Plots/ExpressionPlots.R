# Plots to represent the gene expression results for AIL BFMI S1xS2
#
# copyright (c) 2018-2021 - Brockmann group - HU Berlin Manuel Delpero 
# 
# first written december, 2019

setwd("C:/Users/Manuel/Desktop/AIL_S1xS2/RAWDATA")

exprliver <- read.csv("liverExpressions.txt", sep = "\t", header = TRUE, check.names = FALSE)
exprGonadalfat <- read.csv("GonExpressions.txt", sep = "\t", header = TRUE, check.names = FALSE)
exprSkeletalmuscle <- read.csv("SkExpressions.txt", sep = "\t", header = TRUE, check.names = FALSE)
exprPankreas <- read.csv("PankreasExpressions.txt", sep = "\t", header = TRUE, check.names = FALSE)

Diffexprliver <- read.csv("liver_significant_ann.txt", sep = "\t", header = TRUE, check.names = FALSE)
DiffexprGonadalfat <- read.csv("gonadalfat_significant_ann.txt", sep = "\t", header = TRUE, check.names = FALSE)
DiffexprSkeletalmuscle <- read.csv("skeletalmuscle_significant_ann.txt", sep = "\t", header = TRUE, check.names = FALSE)
DiffexprPankreas <- read.csv("pankreas_significant_ann.txt", sep = "\t", header = TRUE, check.names = FALSE)

## Volcano Plots
regions <- read.table("QTLregions2212020.txt", sep = "\t", header = TRUE)

# Keep the regions that overlap between the traits that show correlation (f.i. gonadal adipose tissue weight and liver) 
regionsGon <- regions[c(30, 44, 45, 46, 47, 56),] # Gon fat
regionsLiver <- regions[47,] # Liver

library(biomaRt)

bio.mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")

getregion <- function(bio.mart, chr, startpos, endpos) {
  region <- paste0(chr, ":",startpos, ":", endpos)
  cat("function: ", " has region: ", region, "\n")
  res.biomart <- getBM(attributes = c("ensembl_gene_id",                                                    # Things that we want to get from biomart
                                      "chromosome_name", "start_position", "end_position", "strand", 
                                      "external_gene_name", "mgi_id", "mgi_symbol", "mgi_description"), 
                       filters = c("chromosomal_region", "biotype"),                                        # Things that we will use to query biomart
                       values = list(region, "protein_coding"),                                             # The thing that we are querying
                       mart = bio.mart)

  cat("function: ", " has ", nrow(res.biomart), "\n")
  return (res.biomart) 
}

# Get genes in regions for Gon fat
genesGon <- vector("list", nrow(regionsGon))
for(x in 1:nrow(regionsGon)){																			
  if(!is.na(regionsGon[x, "Chr"])){
    genesGon[[x]] <- getregion(bio.mart, regionsGon[x, "Chr"], regionsGon[x, "StartPos"], regionsGon[x, "StopPos"])
	cat(x, " has ", nrow(genesGon[[x]]), "genes\n")
    fname <- paste0("genes_in_", regionsGon[x, "Chr"],"-", regionsGon[x, "StartPos"], ":", regionsGon[x, "StopPos"], as.character(regionsGon[x, "Phenotype"]) , ".txt") 
    #write.table(genes[[x]], file = fname, sep="\t", quote = FALSE, row.names = FALSE)
  }else{
 	cat(x, " has NA region\n")
  }
}

# Get genes in regions for liver
genesLiver <- vector("list", nrow(regionsLiver))
for(x in 1:nrow(regionsLiver)){																			
  if(!is.na(regionsLiver[x, "Chr"])){
    genesLiver[[x]] <- getregion(bio.mart, regionsLiver[x, "Chr"], regionsLiver[x, "StartPos"], regionsLiver[x, "StopPos"])
	cat(x, " has ", nrow(genesLiver[[x]]), "genes\n")
    fname <- paste0("genes_in_", regionsLiver[x, "Chr"],"-", regionsLiver[x, "StartPos"], ":", regionsLiver[x, "StopPos"], as.character(regionsLiver[x, "Phenotype"]) , ".txt") 
    #write.table(genes[[x]], file = fname, sep="\t", quote = FALSE, row.names = FALSE)
  }else{
 	cat(x, " has NA region\n")
  }
}

# figure out all the unique genes, result: matrix with 4 columns: name, chromosome, start position, end position
uniquegenesGon <- NULL
for(x in genesGon){ 
  if(!is.null(x)){
    subset <- x[ ,c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "mgi_symbol")]
	uniquegenesGon <- rbind(uniquegenesGon, subset) 
  }
}
uniquegenesGon <- uniquegenesGon[!duplicated(uniquegenesGon),] 
table(uniquegenesGon[ ,"chromosome_name"])

uniquegenesLiver <- NULL
for(x in genesLiver){ 
  if(!is.null(x)){
    subset <- x[ ,c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "mgi_symbol")]
	uniquegenesLiver <- rbind(uniquegenesLiver, subset) 
  }
}
uniquegenesLiver <- uniquegenesLiver[!duplicated(uniquegenesLiver),] 
table(uniquegenesLiver[ ,"chromosome_name"])

# Volvano plot function	(expression values and annotation of the genes in the QTL region as arguments)	 
VolcanoPlot <- function(x,y){       
  sig <- x[which((x[,"logFC"] < -0.2 | x[,"logFC"] > 0.2) & (x[,"p.value"] < 0.05/5000)),]
  int1 <- x[which((x[,"logFC"] < -0.2 | x[,"logFC"] > 0.2) & (x[,"p.value"] > 0.05/5000)),]
  int2 <- x[which((x[,"logFC"] > -0.2 | x[,"logFC"] < 0.2) & (x[,"p.value"] < 0.05/5000) & (!(x[,1] %in% sig[,1]))),]
  not <- x[which((x[,"logFC"] > -0.2 | x[,"logFC"] < 0.2) & (x[,"p.value"] > 0.05/5000) & (!(x[,1] %in% int1[,1]))),]
  myPchsig <- ifelse(sig[,1] %in% y[,1], 16, 1)
  myCexsig <- ifelse(sig[,1] %in% y[,1], 2, 1)
  myPchint1 <- ifelse(int1[,1] %in% y[,1], 16, 1)
  myPchint2 <- ifelse(int2[,1] %in% y[,1], 16, 1)
  title <- readline(prompt="Enter name of the tissue: ")
  plot(main = paste0("Volcano Plot - ", title), x = c(-1,1), y = c(0,15), ylab="-log10[pvalue]", xlab="log2 Fold change", ylim = c(0, 17), xlim = c(-0.8,0.8), t = "n")
    points(sig[, "logFC"], -log10(sig[, "p.value"]), col = "green", pch= myPchsig, cex = myCexsig)
    points(not[, "logFC"], -log10(not[, "p.value"]), col = "red", pch= 1, cex = 1)
    points(int1[, "logFC"], -log10(int1[, "p.value"]), col = "orange", pch= myPchint1, cex = 1)
    points(int2[, "logFC"], -log10(int2[, "p.value"]), col = "orange", pch= myPchint2, cex = 1)
    textToPlot1 <- sig[which(sig[,"ensembl_gene_id"] %in% y[,"ensembl_gene_id"]),]
	textToPlot2 <- int1[which(int1[,"ensembl_gene_id"] %in% y[,"ensembl_gene_id"]),]
	textToPlot3 <- int2[which(int2[,"ensembl_gene_id"] %in% y[,"ensembl_gene_id"]),]
    names <- c(as.character(textToPlot1[,3]), as.character(textToPlot2[,3]), as.character(textToPlot3[,3]))
    posx <- c(textToPlot1[, "logFC"], textToPlot2[, "logFC"], textToPlot3[, "logFC"])
    posy <- c(-log10(textToPlot1[, "p.value"]), -log10(textToPlot2[, "p.value"]), -log10(textToPlot3[, "p.value"]))
    text(posx, posy, names)  
    legend("topright",  bg="gray",
      legend=c("Significant", "Interesting", "Not significant","In QTLs","Out of QTLs"), 
      col=c("green","orange", "red","black", "black"), pch=c(16, 16, 16, 16, 1),
      bty="n", border=F, ncol=2)
}

par(cex.lab=1.2, cex.main = 1.8, cex.axis = 1.4)
par(mfrow = c(1,2))

VolcanoPlot(exprGonadalfat, uniquegenesGon)
VolcanoPlot(exprliver, uniquegenesGon)
