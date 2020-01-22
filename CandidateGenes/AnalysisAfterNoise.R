# Analysis in diff. express genes BFMI S1xS2
#
# copyright (c) 2015-2020 - Brockmann group - HU Berlin, Manuel Delpero
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
DiffExprLiver <- getSignificant(expressions, "L")
DiffExprPank <- getSignificant(expressions, "P")
DiffExprMuscle <- getSignificant(expressions, "S")