# KASP for AIL BFMI S1xS2
#
# copyright (c) 2018-2021 - Brockmann group - HU Berlin Manuel Delpero 
# 
# first written december, 2019

setwd("C:/Users/Manuel/Desktop/AIL_S1xS2/RAWDATA")
library(biomaRt)

map <- read.table("map.cleaned.txt", sep="\t")
lods <- read.table("lods.txt", sep = "\t",  header = TRUE)
chromosomes <- c(1:19, "X", "Y")
map <- map[,c(1,2)]

# Get the main QTLs for each phenotype
res <- c()
for (x in phenotypes){
  Lodscores <- lods[x,]
  ord <- sort(Lodscores, decreasing = TRUE)
  nqtl <- 0
  while ((any(ord > 4)) && (nqtl < 3)){
    info <- cbind(map[names(ord[1]),], rownames(ord[1]), ord[1][[1]])   # get rid of the markers as rownames!!
    markers <- rownames(map[which(map[, "chr"] == info[, "chr"]),])    
    ord <- ord[-which(names(ord) %in% markers)]
    nqtl <- nqtl + 1
    res <- rbind(res, info)
  }
}
colnames(res) <- c("Chr", "Pos", "Trait", "Lod")

res <- res[-which(is.na(res[, "Chr"] )),]

write.table(res, file = "KASP_regions.txt", sep = "\t", quote = FALSE, row.names = FALSE)
QTLs <- read.table("KASP_regions.txt", sep = "\t",  header = TRUE)

# Keep just the unique markers
QTLs <- QTLs[-which(duplicated(QTLs[,5])),]
write.table(res, file = "UniqueKASP_regions.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Get the regions for each QTL
regions <- c()
for (x in 1:nrow(QTLs)){
  startpos <- QTLs[, "Pos"] -120
  endpos <- QTLs[, "Pos"] + 120
  chr <- QTLs[, "Chr"]
  regions <- cbind(startpos, endpos, chr)
}

colnames(regions) <- c("StartPos", "StopPos", "Chr")
rownames(regions) <- QTLs[,"Trait"]

bio.mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")

getregion <- function(bio.mart, chr, startpos, endpos) {
  region <- paste0(chr, ":",startpos, ":", endpos)
  cat("function: ", " has region: ", region, "\n")
  res.biomart <- getSequence(attributes = c("start_position", "end_position",  
                                      ), 
                       filters = c(),                                       
                       values = list(region, ,                                           
                       mart = bio.mart)

  cat("function: ", " has ", nrow(res.biomart), "\n")
  return (res.biomart) 
}

# Get seq in regions
QTLregions <- vector("list", nrow(regions))
for(x in 1:nrow(regions)){																				
  if(!is.na(regions[x, "Chr"])){
    genes[[x]] <- getregion(bio.mart, regions[x, "Chr"], regions[x, "StartPos"], regions[x, "StopPos"])
    fname <- paste0("genes_in_", regions[x, "Chr"],"-", regions[x, "StartPos"], ":", regions[x, "StopPos"], ".txt")
    write.table(genes[[x]], file = fname, sep="\t", quote = FALSE)
  }else{
 	cat(x, " has NA region\n")
  }
}

  