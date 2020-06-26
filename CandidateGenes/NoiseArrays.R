#
# Figuring out Noise versus Signal, get rid of the noise in microarray data
#

setwd("C:/Users/Manuel/Desktop/AIL_S1xS2/Microarray Data/DN-2019_8745-Data")
arraymapping <- read.table("mapping.txt", sep = '\t', header=TRUE, colClasses = "character", row.names=1)

setwd("C:/Users/Manuel/Desktop/AIL_S1xS2/Microarray Data/DN-2019_8745-Data/Rohdaten")

library(affy)

dat <- ReadAffy(cdfname ='clariomsmousemmenstcdf') # Use the clariomsmouse CDF from http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/23.0.0/enst.asp
eset <- mas5(dat)

library("affyPLM")
rawexpr <- assayData(eset)$exprs

sampleexpr <- rawexpr[,1]

plot(c(0, 1000), c(0,1000), t = 'n')
hist(sampleexpr, breaks = 10000, add = TRUE)

getNoiseLvl <- function(sampleexpr){
  require(mixtools)
  # Do a histogram
  histdata <- hist(sampleexpr, breaks = 15000#, plot = FALSE)
  out <- normalmixEM(sampleexpr, c(.5,.5))
  # Get the maximum value of the histogram -> This is the mean of the normal for the noise level
  estmean <- histdata$breaks[which.max(histdata$counts)]

  # Return the max value so we can use it later on to delete intensities that are below this value
  return(c(estmean * 2, min(out$mu) + min(out$sigma)))
}

cutoffs <- apply(rawexpr, 2, getNoiseLvl)

# Use the cutoffs and remove everything below it, set them to 1 <- important use 1 not NA
for(i in 1:ncol(cutoffs)){
  rawexpr[rawexpr[,i] < cutoffs[2, i],i] <- 1
}

# Do the log2 transformation, and afterwards replace the 0s with NA
tenderexpr <- log2(rawexpr)
tenderexpr[tenderexpr == 0] <- NA

plot(apply(tenderexpr, 2, function(x){sum(!is.na(x))}))

# Do normalization as normal, and the rest of the analysis
expressions <- normalize.quantiles(tenderexpr)
rownames(expressions) <- rownames(assayData(eset)$exprs)
colnames(expressions) <- colnames(assayData(eset)$exprs)

ids <- unlist(lapply(strsplit(colnames(expressions), "_"),"[",2))

liver <- which(grepl("L", ids))
gonadalfat <- which(grepl("G", ids))
skeletalmuscle <- which(grepl("S", ids))
pankreas <- which(grepl("P", ids))

colnames(expressions) <- arraymapping[ids, 1]
rownames(expressions) <- gsub("_at", "", rownames(expressions))

corM <- cor(expressions, use="pair")
heatmap(corM, scale = "none")

cesCluster <- c("ENSMUSG00000071047", "ENSMUSG00000078964", "ENSMUSG00000057400", "ENSMUSG00000056973", "ENSMUSG00000061959", "ENSMUSG00000031725", "ENSMUSG00000074156", "ENSMUSG00000058019")

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

expressions <- annotate(expressions)
write.table(expressions, file = "expression_Noise.txt", sep = "\t", row.names = FALSE)
