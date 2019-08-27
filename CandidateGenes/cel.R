setwd("C:/Users/Manuel/Desktop/AIL_S1xS2/Microarray Data/DN-2019_8745-Data")
arraymapping <- read.table("mapping.txt", sep = '\t', header=TRUE, colClasses = "character", row.names=1)

setwd("C:/Users/Manuel/Desktop/AIL_S1xS2/Microarray Data/DN-2019_8745-Data/Rohdaten")

library(affy)

dat <- ReadAffy(cdfname ='clariomsmousemmensgcdf') # Use the clariomsmouse CDF from http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/23.0.0/enst.asp
eset <- mas5(dat)

library("affyPLM")

expressions <- log2(assayData(eset)$exprs)
expressions <- normalize.quantiles(expressions)
rownames(expressions) <- rownames(assayData(eset)$exprs)
colnames(expressions) <- colnames(assayData(eset)$exprs)

ids <- unlist(lapply(strsplit(colnames(expressions), "_"),"[",2))

liver <- which(grepl("L", ids))
gonadalfat <- which(grepl("G", ids))
skeletalmuscle <- which(grepl("S", ids))
pankreas <- which(grepl("P", ids))

colnames(expressions) <- arraymapping[ids, 1]

corM <- cor(expressions)
heatmap(corM, scale = "none")

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

liver <- getSignificant(expressions, "L", p.val = 1.1)
gonadalfat <- getSignificant(expressions, "G", p.val = 1.1)
skeletalmuscle <- getSignificant(expressions, "S", p.val = 1.1)
pankreas <- getSignificant(expressions, "P", p.val = 1.1)

write.table(liver, "liver_all.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gonadalfat, "gonadalfat_all.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(skeletalmuscle, "skeletalmuscle_all.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(pankreas, "pankreas_all.txt", sep = "\t", quote = FALSE, row.names = FALSE)

liver <- getSignificant(expressions, "L", p.val = 0.01)
gonadalfat <- getSignificant(expressions, "G", p.val = 0.01)
skeletalmuscle <- getSignificant(expressions, "S", p.val = 0.01)
pankreas <- getSignificant(expressions, "P", p.val = 0.01)

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

liver <- annotate(liver)
gonadalfat <- annotate(gonadalfat)
skeletalmuscle <- annotate(skeletalmuscle)
pankreas <- annotate(pankreas)

write.table(liver, "liver_significant_ann.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gonadalfat, "gonadalfat_significant_ann.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(skeletalmuscle, "skeletalmuscle_significant_ann.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(pankreas, "pankreas_significant_ann.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Volcano plots
plot(gonadalfat[, "logFC"], -log10(gonadalfat[, "p.value"]))
plot(liver[, "logFC"], -log10(liver[, "p.value"]))
plot(skeletalmuscle[, "logFC"], -log10(skeletalmuscle[, "p.value"]))
plot(pankreas[, "logFC"], -log10(pankreas[, "p.value"]))

subset_gonadalfat <- gonadalfat[which(abs(gonadalfat[, "logFC"]) > 0.4),]
subset_liver <- liver[which(abs(liver[, "logFC"]) > 0.25),]
subset_skeletalmuscle <- skeletalmuscle[which(abs(skeletalmuscle[, "logFC"]) > 0.25),]
subset_pankreas <- gonadalfat[which(abs(pankreas[, "logFC"]) > 0.25),]
