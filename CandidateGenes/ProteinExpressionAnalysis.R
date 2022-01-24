# Diff expression analysis using proteomics data from parental lines BFMI-S1 and BFMI-S2

setwd("C:/Users/Manuel/Desktop/AIL_S1xS2/Proteomics/report/data/")

expressions = read.csv("Quantification/dia-protein-Summary.txt", sep = "\t", check.names = FALSE, header = TRUE, row.names = 1)
annotation <- read.csv("Identification/annotation_allprotein.txt", sep = "\t", check.names = FALSE, header = TRUE, row.names = 1)

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
  return(significant)
}

res <- getSignificant(expressions, Tissue = "Fat")
resLiver <- getSignificant(expressions, Tissue = "Liver")

# combine res with annotation
annot <- annotation[,c("Nr_Identity", "Swissprot Description")]
annot <- annot[grep("GN", annot[,2]),]
annotGon <- annot[rownames(res),]
annotLiver <- annot[rownames(resLiver),]


resAnnotGon <- cbind(annotGon, res)
resAnnotLiver <- cbind(annotLiver, resLiver)

write.csv(resAnnotGon, file = "Identification/DiffExprProteins_GonFat.csv", quote = FALSE)
write.csv(resAnnotGon, file = "Identification/DiffExprProteins_Liver.csv", quote = FALSE)