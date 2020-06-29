# Get a VCF file with SNPs in genes in QTL regions  
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

bamfiles <- c("/home/manuel/AIL_S1xS2/DNA/chr17153_S12.bam",  # 860-S12  (high coverage) #create subset
              "/home/manuel/AIL_S1xS2/DNA/chr17153_S1.bam",    # 861-S1 (medium coverage)
              "/home/manuel/AIL_S1xS2/DNA/chr17153_S2.bam")    # 861-S2 (medium coverage
			 

setwd("/home/manuel/AIL_S1xS2/RAWDATA/SNPsGenesGonLiver/SNPsGenesGonLiver/")

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

# Sort by chromosomes and position otherwise VEP is complaining
allSNPs <- allSNPs[order(as.numeric(allSNPs[,3])),]
chromosomes <- c(3, 15, 17)

annotation <- c()
for(chr in chromosomes){
  annotation <- rbind(annotation, allSNPs[as.numeric(allSNPs[,"V1"]) == chr,])
}

allSNPs <- annotation

header = readLines(filelist[1], n = 169)
cat(paste0(header, collapse = "\n"), "\n", file = "all_combined.vcf")

# File containing all SNPs in the genes
write.table(allSNPs[,-1], file = "all_combined.vcf", sep = "\t", quote=FALSE, append = TRUE, col.names=FALSE, row.names= FALSE)

# Combine these results with the gene expression analysis, look first at genes that are differentially expressed in the QTL regions, then look at VEP results
write.table(allSNPs[,-1], file = "all_combined.vcf", sep = "\t", quote=FALSE, append = TRUE, col.names=FALSE, row.names= FALSE)

# Snps in regions containing the genes (needed to figure out the genotypes of S1 and S2 for the topmarkers)
regions <- read.table("QTLregions2212020.txt", sep = "\t", header = TRUE)
regions <- regions[c(39, 40, 41, 48, 49, 50, 51, 52 , 55, 56),]
for(x in 1:nrow(regions)){ 
  startpos <- regions[x, 4] -10
  endpos <- regions[x, 4] +10
  callSNPs(bamfiles, regions[x, 2], startpos, endpos, paste0(regions[x, 1], "Chr", regions[x,2])) 
}
