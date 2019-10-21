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

DiffExprGon <- getSignificant(GonFatExpr, G)

if (x == 1){
genes <- read.table("gonadalfat_significant_ann.txt", sep = "\t", header = TRUE)
uniquegenes <-  genes[, c("external_gene_name", "chromosome_name", "start_position", "end_position", "strand")] 

# Bam files to use for the SNPs calling function
bamfiles <- c("/halde/BFMI_Alignment_Mar19/merged_sorted_860-S12.bam",  # 860-S12  (high coverage)
             "/halde/BFMI_Alignment_Mar19/merged_sorted_861-S1.bam",    # 861-S1 (medium coverage)
             "/halde/BFMI_Alignment_Mar19/merged_sorted_861-S2.bam")    # 861-S2 (medium coverage)

# Snps in diff. exp. genes
  for(x in 1:nrow(uniquegenes)){ 
    startpos <- uniquegenes[x, 3]
    endpos <- uniquegenes[x, 4]
    if(uniquegenes[x, 5] == 1){
      startpos <- startpos-500
    }else{
      endpos <- endpos +500
    }
    callSNPs(bamfiles, uniquegenes[x, 2], startpos, endpos, uniquegenes[x, 1]) 
  }
}