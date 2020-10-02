# Decision tree for prioritization of the candidate genes in the associated regions
#
# copyright (c) 2018-2021 - Brockmann group - HU Berlin Manuel Delpero
# 

setwd("C:/Users/Manuel/Desktop/AIL_S1xS2/RAWDATA/SNPsGenesGonLiver/SNPsGenesMetS")

FullList <- read.table("annotationSNPsCandidateGenes_allTrig.txt", sep = "\t", check.names = FALSE, header=TRUE)
setwd("C:/Users/Manuel/Desktop/AIL_S1xS2/RAWDATA")
Diffexprliver <- read.csv("DiffExprLiver.txt", sep = "\t", header = TRUE, check.names = FALSE)
DiffexprGonadalfat <- read.csv("DiffExprGon.txt", sep = "\t", header = TRUE, check.names = FALSE)
GenesInfo <- read.csv("genesInfoTrig.txt", sep = "\t", header = TRUE, check.names = FALSE)
#DiffexprSkeletalmuscle <- read.csv("DiffExprMuscle.txt", sep = "\t", header = TRUE, check.names = FALSE)
#DiffexprPankreas <- read.csv("DiffExprPankreas.txt", sep = "\t", header = TRUE, check.names = FALSE)
FullList[, "GENE"] <- as.character(FullList[, "GENE"])
FullList[, "GENE"] <- gsub(" ", "", FullList[, "GENE"])
Candidates <- as.character(unique(FullList[, "GENE"]))
Candidates <- Candidates[-grep("ENSM", Candidates)]
Candidates <- Candidates[-1]
#FullList <- FullList[which(FullList[, "GENE"] %in% Candidates),]
FullList[, "DOMAIN"] <- as.character(FullList[, "DOMAIN"])
FullList[which(FullList[, "DOMAIN"] == ""), "DOMAIN"] <- NA
GenesInfo[, "mgi_symbol"] <-  as.character(GenesInfo[, "mgi_symbol"])

#write.table(FullList, file="fulllist.txt", sep = "\t", quote = FALSE, row.names=FALSE, na = "")

# Download genes in interesting pathways
library("StarBioTrek")
species="mmusculus"
pathwaydb="kegg"
path<-GetData(species,pathwaydb)

# Select interesting pathways
fattyAcidsPathways1 <- path[grep("Fatty", path)]
fattyAcidsPathways2 <- path[grep("Fat", path)]
CholersterolPathways <- path[grep("Cholest", path)]
LipidPathways <- path[grep("lipid", path)]
SteroidPathways <- path[grep("Steroid", path)]
CarboPathways1 <- path[grep("Glyco", path)]
CarboPathways2 <- path[grep("Carbo", path)]
InsulinPathways <- path[grep("Insulin", path)]
CortisolPathways <- path[grep("Cortisol", path)]
ArgininPathways <- path[grep("Arginine", path)]
AlaninePathways <- path[grep("Alanine", path)]
AminoAcidPathways <- path[grep("Amino", path)]
pathways <- c(fattyAcidsPathways1, fattyAcidsPathways2, CholersterolPathways, LipidPathways, SteroidPathways, CarboPathways1, CarboPathways2, InsulinPathways, CortisolPathways,ArgininPathways, AlaninePathways, AminoAcidPathways)
pathway_ALLGENE<-GetPathData(pathways)
pathway_ALLGENE<-ConvertedIDgenes(pathways)
pathway_ALLGENE <- unique(unlist(pathway_ALLGENE, recursive = TRUE, use.names = FALSE))

# Score genes based on mutations (Decision tree)
RankCandidates <- matrix(NA, nrow = length(Candidates), ncol = 10, dimnames = list(Candidates,c("SIFT_del", "SIFT_tol","DOMAIN", "stop_gained", "stop_lost", "splice_donor", "splice_acceptor",  "UTRs", "Promoter", "CTCF B-site", "Enhancer", "Expression", "Annotation", "SCORE")))
for (gene in Candidates) {
  Score <- 0
  GeneInfo <- GenesInfo[which(GenesInfo[, "mgi_symbol"] == gene),]
  geneVar <- FullList[which((FullList[,"POS"] > GeneInfo[, "start_position"]) & (FullList[,"POS"] < GeneInfo[, "end_position"])),] 
  if ((("5_prime" %in% geneVar[, "TYPE"]) || ("3_prime" %in% geneVar[, "TYPE"]))){ # if the gene contain mutation only in a regulatory region
  Score = 1
  RankCandidates[gene, "UTRs"] = 1
  }
  if ((("regulatory" %in% geneVar[, "TYPE"]) && (!("missense_variant" %in% geneVar[, "TYPE"])))){
    if (length(grep("promoter", geneVar[,"GENE"])) > 0){
	  RankCandidates[gene, "Promoter"] = 3
	  Score = Score + 3
	}
	if (length(grep("enhancer", geneVar[,"GENE"])) > 0){
      RankCandidates[gene, "Enhancer"] = 1
	  Score = Score + 1
	}
	if (length(grep("enhancer", geneVar[,"GENE"])) > 0){
	  RankCandidates[gene, "CTCF B-site"] = 1
	  Score = Score + 1
	}
  }
  if (length(grep("stop_gained", geneVar[, "TYPE"]) > 0){
    RankCandidates[gene, "stop_gained"] = 3
    Score = Score + 3 
  }
  if (length(grep("stop_lost", geneVar[, "TYPE"]) > 0){
    RankCandidates[gene, "stop_lost"] = 3
    Score = Score + 3 
  }
  if ((length(grep("missense_variant", geneVar[, "TYPE"]) > 0)) || ("splice_donor_variant" %in% geneVar[, "TYPE"]) ||  ("stop_gained" %in% geneVar[, "TYPE"]) || ("stop_lost" %in% geneVar[, "TYPE"])) { # if the gene contains a mutation in the coding sequence
    if (length(grep("deleterious", geneVar[, "TYPE"]) > 0)){
	  Score <- Score + 3 
	  RankCandidates[gene, "SIFT_del"] = 3
	  }else{
	  Score <- Score + 1
	  RankCandidates[gene, "SIFT_tol"] = 1
	  }
	if (any(!(is.na(geneVar[, "DOMAIN"])))){ 
      Score = Score + 3
	  RankCandidates[gene, "DOMAIN"] = 3
	  }
  }
  if (length(grep("splice_donor", geneVar[, "TYPE"]) > 0){
    RankCandidates[gene, "splice_donor"] = 3
	Score = Score + 3
  }
  if (length(grep("splice_acceptor", geneVar[, "TYPE"]) > 0){
    RankCandidates[gene, "splice_acceptor"] = 3
	Score = Score + 3
  }
  if (length(geneVar[, "GENE"]) > 0) {
    if ((gene %in% DiffexprGonadalfat[, "mgi_symbol"]) || (gene %in% Diffexprliver[, "mgi_symbol"])){ # Check the expressions
	  RankCandidates[gene, "Expression"] = 2
	  Score = Score + 2
	}
  }
  if (gene %in% pathway_ALLGENE){	# Check the annotation
    RankCandidates[gene, "Annotation"] = 1
	Score = Score + 1
  }
  RankCandidates[gene,"SCORE"] <- Score  
}
RankCandidates <- data.frame(RankCandidates[order(as.numeric(RankCandidates[,"SCORE"], decreasing = TRUE)),])


library(biomaRt)
bio.mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
biomart.RankCandidates <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "mgi_symbol"), 
                          filters = c("external_gene_name"), values = rownames(RankCandidates), mart = bio.mart)
rownames(biomart.RankCandidates) <- biomart.RankCandidates[,3]
RankCandidates <- RankCandidates[which(rownames(RankCandidates) %in% rownames(biomart.RankCandidates)),]
biomart.RankCandidates <- biomart.RankCandidates[rownames(RankCandidates),]
RankCandidates <- cbind(biomart.RankCandidates, RankCandidates)

RankCandidate <- RankCandidates[order(-RankCandidates[, "SCORE"], decreasing = FALSE),]

#chromosomes <- c(3,12,15,16,17)
chromosomes <- c(7)
RankCandidates <- c()
for(chr in chromosomes){
  RankCandidates <- rbind(RankCandidates, RankCandidate[RankCandidate[,"chromosome_name"] == chr,])
}


write.table(RankCandidates, file = "CandidatesScoresTrig.txt", sep = "\t", quote = FALSE, row.names = FALSE)