# Get a VCF file with INDELs in genes in QTL regions  
#
# copyright (c) 2018-2021 - Brockmann group - HU Berlin Manuel Delpero & Danny Arends
# 
# first written june, 2020

setwd("/home/manuel/AIL_S1xS2/RAWDATA/")


regions <- read.table("QTLregions2212020.txt", sep = "\t", header = TRUE)

# Keep the regions that overlap between the traits that show correlation (f.i. gonadal adipose tissue weight and liver) 
regions <- regions[c(32,34,35,44,45,47,56),]

# Get genes in regions
genes <- vector("list", nrow(regions))
for(x in 1:nrow(regions)){																			
  if(!is.na(regions[x, "Chr"])){
    genes[[x]] <- getregion(bio.mart, regions[x, "Chr"], regions[x, "StartPos"], regions[x, "StopPos"])
	cat(x, " has ", nrow(genes[[x]]), "genes\n")
    fname <- paste0("genes_in_", regions[x, "Chr"],"-", regions[x, "StartPos"], ":", regions[x, "StopPos"], as.character(regions[x, "Phenotype"]) , ".txt") 
    #write.table(genes[[x]], file = fname, sep="\t", quote = FALSE, row.names = FALSE)
  }else{
 	cat(x, " has NA region\n")
  }
}

# figure out all the unique genes, result: matrix with 4 columns: name, chromosome, start position, end position
uniquegenes <- NULL
for(x in genes){ 
  if(!is.null(x)){
    subset <- x[ ,c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand")]
	uniquegenes <- rbind(uniquegenes, subset) 
  }
}
uniquegenes <- uniquegenes[!duplicated(uniquegenes),] 
table(uniquegenes[ ,"chromosome_name"])

bamfileS1 <- c("/home/danny/NAS/Mouse/DNA/Sequencing/Alignment2020/7/BFMI861-S1P_trimmed.aligned.sorted.dedup.recalibrated.bam",
	           "/home/danny/NAS/Mouse/DNA/Sequencing/Alignment2020/6/BFMI861-S1P_trimmed.aligned.sorted.dedup.recalibrated.bam") 

bamfileS2 <- c("/home/danny/NAS/Mouse/DNA/Sequencing/Alignment2020/9/BFMI861-S2P_trimmed.aligned.sorted.dedup.recalibrated.bam",  
               "/home/danny/NAS/Mouse/DNA/Sequencing/Alignment2020/8/BFMI861-S2P_trimmed.aligned.sorted.dedup.recalibrated.bam")			  

execute <- function(x, intern = FALSE, execute = TRUE){
  cat("----", x, "\n")
  if (execute) {
    res <- system(x, intern = intern)
    cat(">>>>", res[1], "\n")
    if(res[1] >= 1) q("no")
  }
}
		 
## INDELs in genes

# Create the BED file with the regions
bedfile <- c()
for(x in 1:nrow(uniquegenes)){ 
  startpos <- uniquegenes[x, 3]
  endpos <- uniquegenes[x, 4]
  if(uniquegenes[x, 5] == 1){
    startpos <- startpos-500
  }else{
    endpos <- endpos +500
  }
  bedfile <- rbind(bedfile, c(uniquegenes[x, 2], format(startpos, scientific = FALSE), format(endpos, scientific = FALSE)))
}


setwd("/home/manuel/AIL_S1xS2/RAWDATA/INDELsGenesGonLiverS1")
bedfile <- bedfile[-1,]
#bedfile[,1] <- gsub("3", "chr3", bedfile[,1])
#bedfile[,1] <- gsub("12", "chr12", bedfile[,1])
#bedfile[,1] <- gsub("15", "chr15", bedfile[,1])
#bedfile[,1] <- gsub("16", "chr16", bedfile[,1])
#bedfile[,1] <- gsub("17", "chr17", bedfile[,1])
#bedfile[,2] <- as.numeric(bedfile[,2])
#bedfile[,3] <- as.numeric(bedfile[,3])
bedfile <- cbind(as.numeric(bedfile[,1]), as.numeric(bedfile[,2]), as.numeric(bedfile[,3]))
write.table(bedfile, file = "bedfile.bed.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Call INDELs in multiple regions using the BED file
callINDELs <- function(bamfiles) {
  bamstr = paste0(bamfiles, collapse = " ") # Collapse all the bam files in a single string
  scalpel = "/home/florian/Downloads/scalpel-0.5.4/scalpel-discovery --single" # Location of scalpel executable
  reference = "/home/danny/References/Mouse/GRCm38_95/Mus_musculus.GRCm38.dna.toplevel.fa" #Reference genome
  bedfile = "/home/manuel/AIL_S1xS2/RAWDATA/INDELsGenesGonLiverS1/bedfile.bed.txt"
  cmd <- paste0(scalpel, " --bam ", bamstr, " --bed ", bedfile, " --ref ", reference)
  execute(cmd)
  invisible("")
}

setwd("/home/manuel/AIL_S1xS2/RAWDATA/INDELsGenesGonLiverS1")
callINDELs(bamfileS1)
setwd("/home/manuel/AIL_S1xS2/RAWDATA/INDELsGenesGonLiverS2") 
callINDELs(bamfileS2)