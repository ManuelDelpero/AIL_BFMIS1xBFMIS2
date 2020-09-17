setwd("C:/Users/Manuel/Desktop/AIL_S1xS2/RAWDATA/SNPsGenesGonLiver/SNPsGenesMetS")

mdata <- read.table("outputVEP.vcf", header = TRUE)

gts <- apply(mdata[,10:ncol(mdata)],2,function(x){unlist(lapply(strsplit(x, ":"),"[",1)) })
gts[gts=="./."] <- NA

vepp <- unlist(lapply(strsplit(as.character(mdata[,8]), ";"), "[",15))
veppsplit <- strsplit(gsub("CSQ=", "", vepp), ",")
rsIDs <- unlist(lapply(strsplit(unlist(lapply(veppsplit, "[", 1)), "|",fixed=TRUE), "[",18))
domains <- unlist(lapply(strsplit(unlist(lapply(veppsplit, "[", 1)), "|",fixed=TRUE), "[",24))
types <- unlist(lapply(lapply(veppsplit, strsplit, "|", fixed=TRUE), function(x){ 
  sift <- unlist(lapply(strsplit(unlist(lapply(veppsplit, "[", 1)), "|",fixed=TRUE), "[",24))
  type <- unlist(lapply(x,"[",2))
  gene <- unlist(lapply(x,"[",4))
  change <- paste0(unlist(lapply(x,"[",9)), " ", unlist(lapply(x,"[",15)), " ", unlist(lapply(x,"[",16)))
  transcript <- paste0(unlist(lapply(x,"[",7)), " ", unlist(lapply(x,"[",8)))
  for(x in 1:length(type)){
    #if(type[x] %in% c("missense_variant")) { return(paste0(type[x], "\t", " ", change[x], " ", gene[x], "")) }
	if(type[x] %in% c("missense_variant")) { return(paste0(type[x], "\t", " ",gene[x], "", sift[x])) }
    if(type[x] %in% c("regulatory_region_variant")) { return(paste0("regulatory\t", transcript[x], "")) }
    if(type[x] %in% c("splice_acceptor_variant")) { return(paste0("splice_acceptor\t", gene[x], "")) }
    if(type[x] %in% c("splice_donor_variant")) { return(paste0("splice_donor\t", gene[x], "")) }
    if(type[x] %in% c("3_prime_UTR_variant")) { return(paste0("3_prime\t", gene[x], "")) }
    if(type[x] %in% c("5_prime_UTR_variant")) { return(paste0("5_prime\t", gene[x], "")) }
    #if(type[x] %in% c("synonymous_variant")) { return(paste0("synonymous\t", gene[x], "")) }
    if(type[x] %in% c("stop_gained")) { return(paste0("stop_gained\t", gene[x], "")) }
    if(type[x] %in% c("stop_lost")) { return(paste0("stop_lost\t", gene[x], "")) }
    if(type[x] %in% c("start_lost")) { return(paste0("start_lost\t", gene[x], "")) }
  }
  return("\t")
}))

mexcel <- cbind("RSID" = rsIDs, mdata[, c(1,2,4,5,6)], "TYPE" = types, "DOMAIN" = domains, gts)
#write.table(mexcel, file="annotationSNPsCandidateGenes_all.txt", sep = "\t", quote = FALSE, row.names=FALSE, na = "")

# Only SNPs that are found in the S1 and not in the S2
mexcel <- mexcel[which((mexcel[, "BFMI861.S1"] == "1/1" & mexcel[, "BFMI861.S2"] == "0/0" )),]

write.table(mexcel, file="annotationSNPsCandidateGenes_all.txt", sep = "\t", quote = FALSE, row.names=FALSE, na = "")

# Filter the candidates, get only the important SNPs
mexcel <- mexcel[which(!mexcel[, "TYPE"] == "\t"),]
write.table(mexcel, file="annotationSNPsCandidateGenes_Filtered.txt", sep = "\t", quote = FALSE, row.names=FALSE, na = "")

annotGenes <- read.table("annotationSNPsCandidateGenes_all.txt", sep = "\t", header = TRUE, check.names= FALSE)
candidates <- as.character(unique(annotGenes[,8]))
candidates <- candidates[-grep("ENSMU", candidates)]


annotGenes[which(annotGenes[, "TYPE"] %in% c("missense_variant", "splice_donor", "splice_acceptor", "stop_gained", "stop_lost", "start_lost")) , ]

# Check which candidate genes are diff expressed
candidates <- read.table("CandidatesGenes_Complete.txt", sep = "\t", header = FALSE, check.names= FALSE)
setwd("C:/Users/Manuel/Desktop/AIL_S1xS2/RAWDATA")
GonExpressions <- read.csv("DiffExprGon.txt", sep = "\t", header = TRUE, check.names= FALSE)
GonDiffExpressedCandidates <- GonExpressions[which(as.character(GonExpressions[, 3]) %in% as.character(candidates[,1])),]
write.table(GonDiffExpressedCandidates, file="GonDiffExpressedCandidates.txt", sep = "\t", quote = FALSE, row.names=FALSE, na = "")

# Select the expressions for every candidate in every tissue
GonExpressions <- read.csv("ExprGon.txt", sep = "\t", header = TRUE, check.names= FALSE)
LiverExpressions <- read.csv("ExprLiver.txt", sep = "\t", header = TRUE, check.names= FALSE)
Musclexpressions <- read.csv("ExprMuscle.txt", sep = "\t", header = TRUE, check.names= FALSE)
PankreasExpressions <- read.csv("ExprPank.txt", sep = "\t", header = TRUE, check.names= FALSE)

GonExpressedCandidates <- GonExpressions[which(as.character(GonExpressions[, 3]) %in% as.character(candidates[,1])), ]
row.names(GonExpressedCandidates) <- as.character(GonExpressedCandidates[, "mgi_symbol"] )
GonExpressedCandidates <- GonExpressedCandidates[as.character(candidates[,1]),]
write.table(GonExpressedCandidates[, "p.value"], file = "GonPvalue.txt", quote = FALSE, row.names = FALSE)


LiverExpressedCandidates <- LiverExpressions[which(as.character(LiverExpressions[, 3]) %in% as.character(candidates[,1])),]
row.names(LiverExpressedCandidates) <- as.character(LiverExpressedCandidates[, "mgi_symbol"] )
LiverExpressedCandidates <- LiverExpressedCandidates[as.character(candidates[,1]),]
write.table(LiverExpressedCandidates[, "p.value"], file = "LiverPvalue.txt", quote = FALSE, row.names = FALSE)

MuscleExpressedCandidates <- Musclexpressions[which(as.character(Musclexpressions[, 3]) %in% as.character(candidates[,1])),]
row.names(MuscleExpressedCandidates) <- as.character(MuscleExpressedCandidates[, "mgi_symbol"] )
MuscleExpressedCandidates <- MuscleExpressedCandidates[as.character(candidates[,1]),]
write.table(MuscleExpressedCandidates[, "p.value"], file = "MusclePvalue.txt", quote = FALSE, row.names = FALSE)

PankreasExpressedCandidates <- PankreasExpressions[which(as.character(PankreasExpressions[, 3]) %in% as.character(candidates[,1])),]
row.names(PankreasExpressedCandidates) <- as.character(PankreasExpressedCandidates[, "mgi_symbol"] )
PankreasExpressedCandidates <- PankreasExpressedCandidates[as.character(candidates[,1]),]
write.table(PankreasExpressedCandidates[, "p.value"], file = "PankPvalue.txt", quote = FALSE, row.names = FALSE)