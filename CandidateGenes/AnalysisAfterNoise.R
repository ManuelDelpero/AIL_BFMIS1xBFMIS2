# Analysis in diff. express genes BFMI S1xS2
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
getSignificant <- function(expressions, Tissue = "G", adjust = "BH", p.val = 1){
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

write.table(DiffExprLiver, "liverExpressions.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(DiffExprGon, "GonExpressions.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(DiffExprMuscle, "SkExpressions.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(DiffExprPank, "PankreasExpressions.txt", sep = "\t", quote = FALSE, row.names = FALSE)

setwd("C:/Users/Manuel/Desktop/AIL_S1xS2/RAWDATA/SNPsGenesGonLiver")

# Figured out which genes in the QTL regions for the gonadal fat weight, liver weight and glucose are diff. expressed
myvep <- read.csv("VEP.txt", sep = "\t", header=TRUE, colClasses="character")
colnames(myvep) <- c("Uploaded_variation", "Location", "Allele", "Consequence", "IMPACT", "SYMBOL", "Gene", "Feature_type", "Feature", "BIOTYPE", "EXON", "INTRON", "HGVSc", "HGVSp", "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons", "Existing_variation", "DISTANCE", "STRAND", "FLAGS", "SYMBOL_SOURCE", "HGNC_ID", "TSL", "APPRIS", "REFSEQ_MATCH", "SIFT", "CLIN_SIG", "SOMATIC", "PHENO", "MOTIF_NAME", "MOTIF_POS", "HIGH_INF_POS", "MOTIF_SCORE_CHANGE")

genes <- myvep[, "Gene"]
Diffexprgenesg <- DiffExprGon[which(rownames(DiffExprGon) %in% genes) ,]
Diffexprgenesg <- Diffexprgenes[, c("chromosome_name", "mgi_symbol")]
Diffexprgenesg <- Diffexprgenes[order(as.numeric(Diffexprgenes[, "chromosome_name"])),]
write.table(Diffexprgenes, file = "GenesDiffExpr.txt", sep = "\t", quote = FALSE)

Diffexprgenesl <- DiffExprLiver[which(rownames(DiffExprLiver) %in% genes) ,]
Diffexprgenesl <- Diffexprgenes[, c("chromosome_name", "mgi_symbol")]
Diffexprgenesl <- Diffexprgenes[order(as.numeric(Diffexprgenes[, "chromosome_name"])),]

