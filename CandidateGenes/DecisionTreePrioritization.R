# Decision tree for prioritization of the candidate genes in the associated regions
#
# copyright (c) 2018-2021 - Brockmann group - HU Berlin Manuel Delpero
# 

setwd("C:/Users/Manuel/Desktop/AIL_S1xS2/RAWDATA/SNPsGenesGonLiver/SNPsGenesMetS")

FullList <- read.table("ForPrioritization.txt", sep = "\t", check.names = FALSE, header=TRUE)
setwd("C:/Users/Manuel/Desktop/AIL_S1xS2/RAWDATA")
Diffexprliver <- read.csv("DiffExprLiver.txt", sep = "\t", header = TRUE, check.names = FALSE)
DiffexprGonadalfat <- read.csv("DiffExprGon.txt", sep = "\t", header = TRUE, check.names = FALSE)
GenesInfo <- read.csv("genesInfo.txt", sep = "\t", header = TRUE, check.names = FALSE)
DiffexprGonadalfat <- read.csv("DiffExprGon.txt", sep = "\t", header = TRUE, check.names = FALSE)
mmu00010 <- read.table("kegg/mmu00010_Glycolysis Gluconeogenesis.txt",sep='\t', row.names=2, header = FALSE, fill = TRUE)
mmu04910 <- read.table("kegg/mmu04910_Insulin signaling pathway.txt",sep='\t', row.names=2, header = FALSE, fill = TRUE)
mmu04940 <- read.table("kegg/mmu04940_Type I diabetes mellitus.txt",sep='\t', row.names=2, header = FALSE, fill = TRUE)
mmu04930 <- read.table("kegg/mmu04930_Type II diabetes mellitus.txt",sep='\t', row.names=2, header = FALSE, fill = TRUE)
mmu04950 <- read.table("kegg/mmu04950_MODY.txt",sep='\t', row.names=2, header = FALSE, fill = TRUE)
mmu04931 <- read.table("kegg/mmu04931_Insulin resistance.txt",sep='\t', row.names=2, header = FALSE, fill = TRUE)
mmu04975 <- read.table("kegg/mmu04975_Fat digestion and absorption.txt",sep='\t', row.names=2, header = FALSE, fill = TRUE)
mmu00061 <- read.table("kegg/mmu00061_Fatty acid biosynthesis.txt",sep='\t', row.names=2, header = FALSE, fill = TRUE)
mmu01040 <- read.table("kegg/mmu01040_Biosynthesis of unsaturated fatty acids.txt",sep='\t', row.names=2, header = FALSE, fill = TRUE)
mmu04920 <- read.table("kegg/mmu04920_Adipocytokine signaling pathway.txt",sep='\t', row.names=2, header = FALSE, fill = TRUE)
mmu04923 <- read.table("kegg/mmu04923_Regulation of lipolysis in adipocytes.txt",sep='\t', row.names=2, header = FALSE, fill = TRUE)
mmu04973 <- read.table("kegg/mmu04973_Carbohydrate digestion and absorption.txt",sep='\t', row.names=2, header = FALSE, fill = TRUE)
mmu04974 <- read.table("kegg/mmu04974_Protein digestion and absorption.txt",sep='\t', row.names=2, header = FALSE, fill = TRUE)
GenesInPathways <- c(rownames(mmu00010),rownames(mmu04910),rownames(mmu04940),rownames(mmu04930),rownames(mmu04950),rownames(mmu04931),rownames(mmu04975),rownames(mmu00061),rownames(mmu01040),rownames(mmu04920),rownames(mmu04923),rownames(mmu04973),rownames(mmu04974))
GenesInPathways <- unique(GenesInPathways)
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


# Score genes based on mutations (Decision tree)
RankCandidates <- matrix(NA, nrow = length(Candidates), ncol = 9, dimnames = list(Candidates,c("AAchange", "UTRs", "Promoter", "CTCF B-site", "Enhancer", "DOMAIN", "Expression", "Annotation", "SCORE")))
Score <- 0
for (gene in Candidates) {
  GeneInfo <- GenesInfo[which(GenesInfo[, "mgi_symbol"] == gene),]
  geneVar <- FullList[which((FullList[,"POS"] > GeneInfo[, "start_position"]) & (FullList[,"POS"] < GeneInfo[, "end_position"])),] 
  if ((("5_prime" %in% geneVar[, "TYPE"]) || ("3_prime" %in% geneVar[, "TYPE"])) && (!("missense_variant" %in% geneVar[, "TYPE"]))){ # if the gene contain mutation only in a regulatory region
  Score = 3
  RankCandidates[gene, "UTRs"] = "+"
  }
  if ((("regulatory" %in% geneVar[, "TYPE"]) && (!("missense_variant" %in% geneVar[, "TYPE"])))){
  Score = 3
    if (length(grep("promoter", geneVar[,"GENE"])) > 0){
	  RankCandidates[gene, "Promoter"] = "+"
	}
	if (length(grep("enhancer", geneVar[,"GENE"])) > 0){
	RankCandidates[gene, "Enhancer"] = "+"
	}
	if (length(grep("enhancer", geneVar[,"GENE"])) > 0){
	RankCandidates[gene, "CTCF B-site"] = "+"
	}
  }
  if (("missense_variant" %in% geneVar[, "TYPE"]) || ("splice_donor_variant" %in% geneVar[, "TYPE"]) ||  ("stop_gained" %in% geneVar[, "TYPE"]) || ("stop_lost" %in% geneVar[, "TYPE"])) { # if the gene contains a mutation in the coding sequence
    Score <- 3
	if ((("5_prime" %in% geneVar[, "TYPE"]) || ("3_prime" %in% geneVar[, "TYPE"]))){ 
    RankCandidates[gene, "UTRs"] = "+"
	}
    if ((("regulatory" %in% geneVar[, "TYPE"]) && (!("missense_variant" %in% geneVar[, "TYPE"])))){
      if (length(grep("promoter", geneVar[,"GENE"])) > 0){
	    RankCandidates[gene, "Promoter"] = "+"
	  }
	  if (length(grep("enhancer", geneVar[,"GENE"])) > 0){
	    RankCandidates[gene, "Enhancer"] = "+"
	  }
      if (length(grep("enhancer", geneVar[,"GENE"])) > 0){
	  RankCandidates[gene, "CTCF B-site"] = "+"
	  }
    } 
	if (any(!(is.na(geneVar[, "DOMAIN"])))){ # if the gene contains a mutation in a coding sequence located in a domain
      Score = Score + 2
	  RankCandidates[gene, c("AAchange", "DOMAIN")] = "+"
	}else{
	  score = Score
	  RankCandidates[gene, "AAchange"] = "+"
	  }
  }
  if (length(geneVar[, "GENE"]) > 0) {
    if ((gene %in% DiffexprGonadalfat[, "mgi_symbol"]) || (gene %in% Diffexprliver[, "mgi_symbol"])){ # Check the expressions
	  RankCandidates[gene, "Expression"] = "+"
	  Score = Score +2
	}
  }
  if (gene %in% GenesInPathways){	# Check the annotation
    RankCandidates[gene, "Annotation"] = "+"
	Score = Score + 1
  }
  RankCandidates[gene,"SCORE"] <- Score  
}
RankCandidates <- data.frame(RankCandidates[order(as.numeric(RankCandidates[,"SCORE"], decreasing = TRUE)),])