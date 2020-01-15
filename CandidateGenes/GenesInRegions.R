# Get SNPs in genes in regions
#
# copyright (c) 2018-2021 - Brockmann group - HU Berlin Manuel Delpero & Danny Arends
# 
# first written december, 2019

setwd("/home/manuel/AIL_S1xS2/DATA")

regions <- read.table("QTL_regions11420.txt", sep = "\t", header = TRUE)

# Get genes in regions
genesBW <- vector("list", nrow(regions[1:59,]))
for(x in 1:nrow(regions[1:59,])){																				#for loop number of rows in the region
  # check if the region is valid (not NA)
  if(!is.na(regions[x, "Chr"])){
    # get the genes in the region
    genesBW[[x]] <- getregion(bio.mart, regions[x, "Chr"], regions[x, "StartPos"], regions[x, "StopPos"])
	cat(x, " has ", nrow(genesBW[[x]]), "genes\n")
    fname <- paste0("genes_in_", regions[x, "Chr"],"-", regions[x, "StartPos"], ":", regions[x, "StopPos"], regions[x, "Phenotype"] , ".txt")
    write.table(genesBW[[x]], file = fname, sep="\t", quote = FALSE, row.names = FALSE)
  }else{
 	cat(x, " has NA region\n")
  }
}

# figure out all the unique genes, result: matrix with 4 columns: name, chromosome, start position, end position
uniquegenesBW <- NULL
for(x in genesBW){ 
  if(!is.null(x)){
    subset <- x[ ,c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand")]
	uniquegenesBW <- rbind(uniquegenesBW, subset) 
  }
}
uniquegenesBW <- uniquegenesBW[!duplicated(uniquegenesBW),] 
table(uniquegenesBW[ ,"chromosome_name"])

bamfiles <- c("/halde/BFMI_Alignment_Mar19/merged_sorted_860-S12.bam",  # 860-S12  (high coverage)
             "/halde/BFMI_Alignment_Mar19/merged_sorted_861-S1.bam",    # 861-S1 (medium coverage)
             "/halde/BFMI_Alignment_Mar19/merged_sorted_861-S2.bam")    # 861-S2 (medium coverage)

# Snps in genes
for(x in 1:nrow(uniquegenesBW)){ 
  startpos <- uniquegenesBW[x, 3]
  endpos <- uniquegenesBW[x, 4]
  if(uniquegenesBW[x, 5] == 1){
    startpos <- startpos-500
  }else{
    endpos <- endpos +500
  }
  callSNPs(bamfiles, uniquegenesBW[x, 2], startpos, endpos, uniquegenesBW[x, 1]) 
}

setwd("/home/manuel/AIL_B6xBFMI/RAWDATA/SNPsMQM")

filelist <- list.files(".") #we removed all the .txt files

allSNPs <- NULL
for(file in filelist){
  if(length(readLines(file)) > 169){
    mcontent <- read.csv(file, sep = "\t", skip = 169, header=FALSE, colClasses=c("character"))
    gene = gsub(".snps-filtered.vcf", "", file, fixed=T)
    allSNPs <- rbind(allSNPs, cbind(gene, mcontent))
  }
}

header = readLines(filelist[1], n = 169)
cat(paste0(header, collapse = "\n"), "\n", file = "all_combined.vcf")
# File containing all SNPs in the genes
write.table(allSNPs[,-1], file = "all_combined.vcf", sep = "\t", quote=FALSE, append = TRUE, col.names=FALSE, row.names= FALSE)

q("no")