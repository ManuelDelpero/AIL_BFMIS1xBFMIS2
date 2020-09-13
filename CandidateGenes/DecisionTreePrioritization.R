# Decision tree for prioritization of the candidate genes in the associated regions
#
# copyright (c) 2018-2021 - Brockmann group - HU Berlin Manuel Delpero & Danny Arends
# 

setwd("C:/Users/Manuel/Desktop/AIL_S1xS2/RAWDATA/SNPsGenesGonLiver/SNPsGenesMetS")

FullList <- read.table("ForPrioritization.txt", sep = "\t", check.names = FALSE, header=TRUE)
setwd("C:/Users/Manuel/Desktop/AIL_S1xS2/RAWDATA")
Diffexprliver <- read.csv("DiffExprLiver.txt", sep = "\t", header = TRUE, check.names = FALSE)
DiffexprGonadalfat <- read.csv("DiffExprGon.txt", sep = "\t", header = TRUE, check.names = FALSE)
#DiffexprSkeletalmuscle <- read.csv("DiffExprMuscle.txt", sep = "\t", header = TRUE, check.names = FALSE)
#DiffexprPankreas <- read.csv("DiffExprPankreas.txt", sep = "\t", header = TRUE, check.names = FALSE)
FullList[, "GENE"] <- as.character(FullList[, "GENE"])
FullList[, "GENE"] <- gsub(" ", "", FullList[, "GENE"])
Candidates <- as.character(unique(FullList[, "GENE"]))
Candidates <- Candidates[-grep("ENSM", Candidates)]
Candidates <- Candidates[-1]
FullList <- FullList[which(FullList[, "GENE"] %in% Candidates),]
FullList[, "DOMAIN"] <- as.character(FullList[, "DOMAIN"])
FullList[which(FullList[, "DOMAIN"] == ""), "DOMAIN"] <- NA
write.table(FullList, file="fulllist.txt", sep = "\t", quote = FALSE, row.names=FALSE, na = "")

# Score genes based on mutations 
RankCandidates <- matrix(NA, nrow = length(Candidates), ncol = 6, dimnames = list(Candidates,c("AAchange", "Regulatory", "DOMAIN", "Expression", "Annotation", "SCORE")))
Score <- 0
for (gene in Candidates) {
  geneVar <- FullList[which(FullList[,"GENE"] == gene),] # first branch, if the gene contains a mutation (regulatory regions (UTR) or AA change mutation)
  if ((("5_prime" %in% geneVar[, "TYPE"]) || ("3_prime" %in% geneVar[, "TYPE"])) && (!("missense_variant" %in% geneVar[, "TYPE"]))){ # if the gene contain mutation only in a regulatory region
  Score = 3
  RankCandidates[gene, "Regulatory"] = "+"
  }
  if (("missense_variant" %in% geneVar[, "TYPE"]) || ("splice_donor_variant" %in% geneVar[, "TYPE"]) ||  ("stop_gained" %in% geneVar[, "TYPE"]) || ("stop_lost" %in% geneVar[, "TYPE"])) { # if the gene contains a mutation in the coding sequence
    Score <- 3
	if (("5_prime" %in% geneVar[, "TYPE"]) || ("3_prime" %in% geneVar[, "TYPE"])){
	  RankCandidates[gene, "Regulatory"] = "+"
	  }
	if (any(!(is.na(geneVar[, "DOMAIN"])))){ # if the gene contains a mutation in a coding sequence located in a domain
      Score = Score + 2
	  RankCandidates[gene, c("AAchange", "DOMAIN")] = "+"
	}else{
	  score = Score
	  RankCandidates[gene, "AAchange"] = "+"
	  }
  }
  if ((geneVar[, "GENE"] %in% DiffexprGonadalfat[, "mgi_symbol"]) || (geneVar[, "GENE"] %in% Diffexprliver[, "mgi_symbol"])){
	RankCandidates[gene, "Expression"] = "+"
	Score = Score +2
  }
  RankCandidates[gene,"SCORE"] <- Score  
}
RankCandidates <- data.frame(RankCandidates[order(as.numeric(RankCandidates[,"SCORE"], decreasing = TRUE)),])
# Score the genes based on annotation  