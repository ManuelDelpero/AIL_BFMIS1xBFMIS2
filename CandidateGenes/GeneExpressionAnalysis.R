setwd("C:/Users/Manuel/Desktop/AIL_S1xS2/Microarray Data/DN-2019_8745-Data")
arraymapping <- read.table("mapping.txt", sep = '\t', header=TRUE, colClasses = "character", row.names=1)

setwd("C:/Users/Manuel/Desktop/AIL_S1xS2/Microarray Data/DN-2019_8745-Data/Rohdaten")

library(affy)
library("gplots")

dat <- ReadAffy(cdfname ='ClariomSMouse_Mm_ENST') # Use the clariomsmouse CDF from http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/23.0.0/enst.asp
eset <- mas5(dat)

library("affyPLM")

# QC 
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

# Heatmap

col_labels <- c("orange","orange", "orange","orange","blue", "blue", 
                "blue","blue", "blue","orange","orange","orange","blue",
				"blue","blue", "blue", "orange", "blue", "blue", "orange",
				"orange", "orange","orange","orange","orange", "blue", "blue",
				"blue", "blue", "blue", "blue","blue", "blue", "orange","orange",
				"blue","blue","blue","orange","orange", "orange","orange","orange",
				"orange", "orange","orange","orange","orange", "orange","blue","blue",
				"blue","blue","blue","blue","blue","blue")

# But due to the way heatmap.2 works - we need to fix it to be in the 
# order of the data!  
colors = unique(c(seq(0.94, 1,length = 30)))

my_palette <- rev(heat.colors(29))

distance= dist(corM, method ="euclidean")    
hcluster = hclust(distance, method ="ward.D")
dend1 <- as.dendrogram(hcluster)
cols_branches <- c("orange", "blue")
dend1 <- color_branches(dend1, k = 4, col = cols_branches)
col_labels <- col_labels[order(order.dendrogram(dend1))]


heatmap.2(corM,
  main = "Correlation", # heat map title
  notecol="black",      # change font color of cell labels to black
  density.info="none",  # turns off density plot inside color legend
  trace="none", 
  col=my_palette,
  breaks = colors,  # turns off trace lines inside the heat map
  #margins =c(12,9),     # widens margins around plot
  dendrogram="col",     # only draw a row dendrogram
  #Rowv="NA",
  ColSideColors = col_labels,
  key.xlab="correlation coefficient"
  )
legend("bottomleft",
  bty = "n",
  legend = c("S1", "S2"),
  pch = c(15, 15),
  col = c("orange", "blue"),
  pt.cex=1.5,
  cex=1.5)


liverexpr <- expressions[, which(grepl("L", colnames(expressions)))]
LiverS1expr <- liverexpr[, which(grepl("S1", colnames(liverexpr)))]
LiverS1expr <- data.frame(apply(LiverS1expr, 1 , mean))
rownames(LiverS1expr) <- gsub("_at", "", rownames(LiverS1expr))
Gonexpr <- expressions[, which(grepl("G", colnames(expressions)))]
GonS1expr <- Gonexpr[, which(grepl("S1", colnames(Gonexpr)))]
GonS1expr <- data.frame(apply(GonS1expr, 1 , mean))
rownames(GonS1expr) <- gsub("_at", "", rownames(GonS1expr))
LiverS1expr_ann <- annotate(LiverS1expr)
GonS1expr_ann <- annotate(GonS1expr)

liverexpr <- expressions[, which(grepl("L", colnames(expressions)))]
LiverS2expr <- liverexpr[, which(grepl("S2", colnames(liverexpr)))]
LiverS2expr <- data.frame(apply(LiverS2expr, 1 , mean))
rownames(LiverS2expr) <- gsub("_at", "", rownames(LiverS2expr))
Gonexpr <- expressions[, which(grepl("G", colnames(expressions)))]
GonS2expr <- Gonexpr[, which(grepl("S2", colnames(Gonexpr)))]
GonS2expr <- data.frame(apply(GonS2expr, 1 , mean))
rownames(GonS2expr) <- gsub("_at", "", rownames(GonS2expr))
LiverS2expr_ann <- annotate(LiverS2expr)
GonS2expr_ann <- annotate(GonS2expr)

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

liver <- annotate(getSignificant(expressions, "L", p.val = 1.1))
gonadalfat <- annotate(getSignificant(expressions, "G", p.val = 1.1))
skeletalmuscle <- annotate(getSignificant(expressions, "S", p.val = 1.1))
pankreas <- annotate(getSignificant(expressions, "P", p.val = 1.1))

write.table(liver, "liver_all_ann.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gonadalfat, "gonadalfat_all_ann.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(skeletalmuscle, "skeletalmuscle_all.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(pankreas, "pankreas_all.txt", sep = "\t", quote = FALSE, row.names = FALSE)

liver <- getSignificant(expressions, "L", p.val = 0.01)
gonadalfat <- getSignificant(expressions, "G", p.val = 0.01)
skeletalmuscle <- getSignificant(expressions, "S", p.val = 0.01)
pankreas <- getSignificant(expressions, "P", p.val = 0.01)

liver <- annotate(liver)
gonadalfat <- annotate(gonadalfat)
skeletalmuscle <- annotate(skeletalmuscle)
pankreas <- annotate(pankreas)

# Figure out which genes are up or down regulated to use in KEGG
upRegS1 <- c()
downRegS1 <- c()
upRegS2 <- c()
downRegS2 <- c()
upS1 <- c()
upS2 <- c()
downS1 <- c()
downS2 <- c()

for (x in 1:nrow(gonadalfat)){
  if (gonadalfat[x, "mean(s1)"] > gonadalfat[x, "mean(s2)"]){
    upS1 <- gonadalfat[x, c("ensembl_gene_id")]
    downS2 <- gonadalfat[x, c("ensembl_gene_id")]
  }else{
    upS2 <- gonadalfat[x, c("ensembl_gene_id")]
	downS1 <- gonadalfat[x, c("ensembl_gene_id")]
    }
  upRegS1 <- c(upRegS1, upS1)
  downRegS1 <- c(downRegS1, downS1)
  upRegS2 <- c(upRegS2, upS2)
  downRegS2 <- c(downRegS2, downS2)
}	
  
write.table(liver, "liver_significant_ann.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gonadalfat, "gonadalfat_significant_ann.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(skeletalmuscle, "skeletalmuscle_significant_ann.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(pankreas, "pankreas_significant_ann.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Volcano plots
plot(gonadalfat[, "logFC"], -log10(gonadalfat[, "p.value"]))
plot(liver[, "logFC"], -log10(liver[, "p.value"]))
plot(skeletalmuscle[, "logFC"], -log10(skeletalmuscle[, "p.value"]))
plot(pankreas[, "logFC"], -log10(pankreas[, "p.value"]))

subset_gonadalfat <- gonadalfat[which(abs(gonadalfat[, "logFC"]) > 0.25),]
subset_liver <- liver[which(abs(liver[, "logFC"]) > 0.25),]
subset_skeletalmuscle <- skeletalmuscle[which(abs(skeletalmuscle[, "logFC"]) > 0.25),]
subset_pankreas <- gonadalfat[which(abs(pankreas[, "logFC"]) > 0.25),]



   
 

