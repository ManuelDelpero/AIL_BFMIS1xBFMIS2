# Get a VCF file with INDELs in genes in QTL regions  
#
# copyright (c) 2018-2021 - Brockmann group - HU Berlin Manuel Delpero & Danny Arends
# 
# first written december, 2019

setwd("/home/manuel/AIL_S1xS2_old/RAWDATA/")


regions <- read.table("QTLregions2212020.txt", sep = "\t", header = TRUE)

# Keep the regions that overlap between the traits that show correlation (f.i. gonadal adipose tissue weight and liver) 
regions <- regions[c(39, 40, 41, 49, 50, 52 , 55, 56),]

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

bamfiles <- c("/home/danny/NAS/Mouse/DNA/Sequencing/Alignment2020/9/BFMI861-S2P_trimmed.aligned.sorted.dedup.recalibrated.bam",  
              "/home/danny/NAS/Mouse/DNA/Sequencing/Alignment2020/8/BFMI861-S2P_trimmed.aligned.sorted.dedup.recalibrated.bam",    
              "/home/danny/NAS/Mouse/DNA/Sequencing/Alignment2020/7/BFMI861-S1P_trimmed.aligned.sorted.dedup.recalibrated.bam",
			  "/home/danny/NAS/Mouse/DNA/Sequencing/Alignment2020/6/BFMI861-S1P_trimmed.aligned.sorted.dedup.recalibrated.bam")   
			 
# INDELs in genes
for(x in 1:nrow(uniquegenes)){ 
  startpos <- uniquegenes[x, 3]
  endpos <- uniquegenes[x, 4]
  if(uniquegenes[x, 5] == 1){
    startpos <- startpos-500
  }else{
    endpos <- endpos +500
  }
  callINDELs(bamfiles, uniquegenes[x, 2], startpos, endpos, uniquegenes[x, 1]) 
}