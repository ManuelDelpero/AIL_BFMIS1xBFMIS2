# Analysis in diff. express genes after noise in BFMI S1xS2
#
# copyright (c) 2018-2020 - Brockmann group - HU Berlin, Manuel Delpero
# last written october, 2019

setwd("C:/Users/Manuel/Desktop/AIL_S1xS2/RAWDATA")

expressions <- read.table(file = "expression_Noise.txt", sep = "\t", header = TRUE, check.names = FALSE)

# Select the different tissues
liver <- which(grepl("L", colnames(expressions)))
gonadalfat <- which(grepl("G", colnames(expressions)))
skeletalmuscle <- which(grepl("S", colnames(expressions)))
pankreas <- which(grepl("P", colnames(expressions)))

# Diff. gene expression analysis
getSignificant <- function(expressions, Tissue = "G", adjust = "BH", p.val = 0.05){
  S1P <- which(grepl("S1", colnames(expressions)) & grepl(Tissue, colnames(expressions)))
  S2P <- which(grepl("S2", colnames(expressions)) & grepl(Tissue, colnames(expressions)))

  res <- t(apply(expressions[, c(S1P, S2P)],1, function(x) {
    s1v <- x[1:length(S1P)]
    s2v <- x[(length(S1P)+1):(length(S1P) + length(S2P))]
    pval <- tryCatch({t.test(s1v, s2v)$p.value}, error = function(e) { return(NA) })
    return(c(mean(s1v), mean(s2v), sd(s1v), sd(s2v), mean(s1v) / mean(s2v), log2(mean(s1v) / mean(s2v)), pval))
  }))
  colnames(res) <- c("mean(s1)", "mean(s2)", "sd(s1)", "sd(s2)", "FC", "logFC", "p.value")
  significant <- res[which(p.adjust(res[,"p.value"], adjust) < p.val),]
  rownames(significant) <- gsub("_at", "", rownames(significant))
  return(significant)
}
# Annotation
annotate <- function(significant){
  library(biomaRt)
  bio.mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")

  res.biomart <- getBM(attributes = c("ensembl_gene_id", "mgi_id", "mgi_symbol", "mgi_description", "chromosome_name", "start_position", "end_position", "strand"), 
                       filters = c("ensembl_gene_id"), 
                       values = rownames(significant), mart = bio.mart)
  rownames(res.biomart) <- res.biomart[,1]
  annotated <- cbind(res.biomart[rownames(significant), ], significant)
  return(annotated)
}

DiffExprGon <- annotate(getSignificant(expressions, "G"))
DiffExprLiver <- annotate(getSignificant(expressions, "L"))
DiffExprPank <- annotate(getSignificant(expressions, "P"))
DiffExprMuscle <- annotate(getSignificant(expressions, "S"))

#write.table(DiffExprLiver, "liverExpressions.txt", sep = "\t", quote = FALSE, row.names = FALSE)
#write.table(DiffExprGon, "GonExpressions.txt", sep = "\t", quote = FALSE, row.names = FALSE)
#write.table(DiffExprMuscle, "SkExpressions.txt", sep = "\t", quote = FALSE, row.names = FALSE)
#write.table(DiffExprPank, "PankreasExpressions.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Get genes in the QTL regions and check whichones are diff expressed in all the tissues
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

regions <- read.table("QTLregions2212020.txt", sep = "\t", header = TRUE)

# Keep the regions that overlap between the traits that show correlation (f.i. gonadal adipose tissue weight and liver) 
regions <- regions[c(25, 44, 47),]

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





