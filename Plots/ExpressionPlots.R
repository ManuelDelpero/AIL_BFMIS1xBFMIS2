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
regions <- regions[c(39, 40, 41, 49, 50, 51, 52 , 55, 56),]

library(biomaRt)

bio.mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")

getregion <- function(bio.mart, chr, startpos, endpos) {
  region <- paste0(chr, ":",startpos, ":", endpos)
  cat("function: ", " has region: ", region, "\n")
  res.biomart <- getBM(attributes = c("ensembl_gene_id",                                                    # Things that we want to get from biomart
                                      "chromosome_name", "start_position", "end_position", "strand", 
                                      "external_gene_name", "mgi_id", "mgi_symbol", "mgi_description"), 
                       filters = c("chromosomal_region", "biotype"),                                        # Things that we will use to query biomart
                       values = list(
region, "protein_coding"),                                             # The thing that we are querying
                       mart = bio.mart)

  cat("function: ", " has ", nrow(res.biomart), "\n")
  return (res.biomart) 
}

# Get genes in regions
genes <- vector("list", nrow(regions))
for(x in 1:nrow(regions)){																			
  if(!is.na(regions[x, "Chr"])){
    genes[[x]] <- getregion(bio.mart, regions[x, "Chr"], regions[x, "StartPos"], regions[x, "StopPos"])
	cat(x, " has ", nrow(genes[[x]]), "genes\n")
    fname <- paste0("genes_in_", regions[x, "Chr"],"-", regions[x, "StartPos"], ":", regions[x, "StopPos"], as.character(regions[x, "Phenotype"]) , ".txt") 
    #write.table(genes[[x]], file = fname, sep="\t", quote = FALSE, row.names = FALSE)
  }else{
 	cat(x, " has NA region\n")
  }
}

# figure out all the unique genes, result: matrix with 4 columns: name, chromosome, start position, end position
uniquegenes <- NULL
for(x in genes){ 
  if(!is.null(x)){
    subset <- x[ ,c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand")]
	uniquegenes <- rbind(uniquegenes, subset) 
  }
}
uniquegenes <- uniquegenes[!duplicated(uniquegenes),] 
table(uniquegenes[ ,"chromosome_name"])

bamfiles <- c("/halde/BFMI_Alignment_Mar19/merged_sorted_860-S12.bam",  # 860-S12  (high coverage)
             "/halde/BFMI_Alignment_Mar19/merged_sorted_861-S1.bam",    # 861-S1 (medium coverage)
             "/halde/BFMI_Alignment_Mar19/merged_sorted_861-S2.bam")    # 861-S2 (medium coverage)


# Gonadal fat
out <- which(exprGonadalfat[, "logFC"] < -1)
exprGonadalfat = exprGonadalfat[-out,]
myPch <- ifelse(exprGonadalfat[,1] %in% uniquegenes[,1], 16, 1)
colz <- as.numeric(exprGonadalfat[,"logFC"] > 0.2) + as.numeric(exprGonadalfat[,"logFC"] < -0.2) + as.numeric(exprGonadalfat[,"p.value"] < 0.05/5000) + as.numeric(exprGonadalfat[,"p.value"] < 0.01/5000) +1
plot(exprGonadalfat[, "logFC"], -log10(exprGonadalfat[, "p.value"]), col=as.numeric(colz == 4) + 2, pch= myPch, ylab="log10(pvalue)", xlab="log2 Fold change")


plot(exprliver[, "logFC"], -log10(exprliver[, "p.value"]))
plot(skeletalmuscle[, "logFC"], -log10(skeletalmuscle[, "p.value"]))
plot(pankreas[, "logFC"], -log10(pankreas[, "p.value"]))

