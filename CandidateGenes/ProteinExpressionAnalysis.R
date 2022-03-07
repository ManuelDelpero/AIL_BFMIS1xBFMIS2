# Diff expression analysis using proteomics data from parental lines BFMI-S1 and BFMI-S2

setwd("C:/Users/Manuel/Desktop/AIL_S1xS2/Proteomics/report/data/")

expressions = read.csv("Quantification/dia-protein-Summary.txt", sep = "\t", check.names = FALSE, header = TRUE, row.names = 1)
annotation <- read.csv("Identification/annotation_allprotein.txt", sep = "\t", check.names = FALSE, header = TRUE, row.names = 1)

# Diff. gene expression analysis
getSignificant <- function(expressions, Tissue = "G", adjust = "BH", p.val = 0.01){
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

resAll <- getSignificant(expressions, Tissue = "Fat", p.val = 1)
resLiverAll <- getSignificant(expressions, Tissue = "Liver", p.val = 1)

# combine res with annotation
annot <- annotation[,c("Nr_Identity", "Swissprot Description")]
annot <- annot[grep("GN", annot[,2]),]
annotGon <- annot[rownames(res),]
annotLiver <- annot[rownames(resLiver),]


resAnnotGon <- cbind(annotGon, res)
resAnnotLiver <- cbind(annotLiver, resLiver)

resAnnotGon <- resAnnotGon[,-2]
resAnnotLiver <- resAnnotLiver[,-2]


#write.table(resAnnotGon, file = "Identification/DiffExprProteins_GonFat.csv", quote = FALSE, sep = ",")
#write.table(resAnnotLiver, file = "Identification/DiffExprProteins_Liver.csv", quote = FALSE, sep = ",")


a <- lapply(strsplit(unlist(rownames(resAnnotLiver)),"|"), '[', 4:9)
unlist(lapply(a, paste, collapse = ""))
#write.table(unlist(lapply(a, paste, collapse = "")), file = "UniprotID_Liver.txt", quote = FALSE, row.names = FALSE)

a <- lapply(strsplit(unlist(rownames(res)),"|"), '[', 4:9)
unlist(lapply(a, paste, collapse = ""))
#write.table(unlist(lapply(a, paste, collapse = "")), file = "UniprotID_Gon.txt", quote = FALSE, row.names = FALSE)

# Volvano plot function	(all expression values (x) and expression values of only diff expressed (z))	 
VolcanoPlot <- function(x,z){       
  sig <- x[which((x[,"logFC"] < -0.2 | x[,"logFC"] > 0.2) & (x[,"p.value"] < max(z[,"p.value"]))),]
  int1 <- x[which((x[,"logFC"] < -0.2 | x[,"logFC"] > 0.2) & (x[,"p.value"] > max(z[,"p.value"]))),]
  int2 <- x[which((x[,"logFC"] > -0.2 | x[,"logFC"] < 0.2) & (x[,"p.value"] < max(z[,"p.value"])) & (!(x[,1] %in% sig[,1]))),]
  not <- x[which((x[,"logFC"] > -0.2 | x[,"logFC"] < 0.2) & (x[,"p.value"] > max(z[,"p.value"])) & (!(x[,1] %in% int1[,1]))),]
  title <- readline(prompt="Enter name of the tissue: ")
  plot(main = paste0("Volcano Plot - ", title), x = c(-1,1), y = c(0,15), ylab="-log10[pvalue]", xlab="log2 Fold change", ylim = c(0, 17), xlim = c(-0.8,0.8), t = "n")
    points(sig[, "logFC"], -log10(sig[, "p.value"]), col = "green")
    points(not[, "logFC"], -log10(not[, "p.value"]), col = "red", cex = 1)
    points(int1[, "logFC"], -log10(int1[, "p.value"]), col = "orange", cex = 1)
    points(int2[, "logFC"], -log10(int2[, "p.value"]), col = "orange", cex = 1)
    legend("topright",  bg="gray",
      legend=c("Significant", "Remarkable", "Not significant"), 
      col=c("green","orange", "red"), pch=c(1, 1, 1),
      bty="n", border=F, ncol=2)
}

VolcanoPlot(resAll, res)