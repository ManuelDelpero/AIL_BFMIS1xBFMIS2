# Get SNPs in genes in regions
#
# copyright (c) 2018-2021 - Brockmann group - HU Berlin Manuel Delpero & Danny Arends
# 
# first written december, 2019

setwd("C:/Users/Manuel/Desktop/AIL_S1xS2/RAWDATA")

regions <- read.table("QTLregions2212020.txt", sep = "\t", header = TRUE)

# Keep the regions that overlap between the traits that show correlation (f.i. gonadal adipose tissue weight and liver) 
regions <- regions[c(39, 40, 41, 49, 50, 51, 52 , 55, 56),]

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

bamfiles <- c("/halde/BFMI_Alignment_Mar19/merged_sorted_860-S12.bam",  # 860-S12  (high coverage)
             "/halde/BFMI_Alignment_Mar19/merged_sorted_861-S1.bam",    # 861-S1 (medium coverage)
             "/halde/BFMI_Alignment_Mar19/merged_sorted_861-S2.bam")    # 861-S2 (medium coverage)


setwd("C:/Users/Manuel/Desktop/AIL_S1xS2/RAWDATA/SNPsGenesGonLiver")

# Snps in genes
for(x in 1:nrow(uniquegenes)){ 
  startpos <- uniquegenes[x, 3]
  endpos <- uniquegenes[x, 4]
  if(uniquegenes[x, 5] == 1){
    startpos <- startpos-500
  }else{
    endpos <- endpos +500
  }
  callSNPs(bamfiles, uniquegenes[x, 2], startpos, endpos, uniquegenes[x, 1]) 
}


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
write.table(allSNPs[,-1], file = "all_combinedBW.vcf", sep = "\t", quote=FALSE, append = TRUE, col.names=FALSE, row.names= FALSE)