# Call SNPs in QTL regions using seq data from the MGP to identify novel SNPs


setwd("/home/manuel/AIL_S1xS2/RAWDATA")

execute <- function(x, intern = FALSE, execute = TRUE){
  cat("----", x, "\n")
  if (execute) {
    res <- system(x, intern = intern)
    cat(">>>>", res[1], "\n")
    if(res[1] >= 1) q("no")
  }
}

callSNPs <- function(bamfiles, chr = 1, startpos = 1, endpos = 2, outname = "mySNPs") {
  bamstr = paste0(bamfiles, collapse = " ") # Collapse all the bam files in a single string
  bcftools = "/home/arends/Github/bcftools/bcftools" # Location of BCFtools executable
  reference = "/home/arends/NAS/Mouse/Reference_Genomes/GRCm38_68/GRCm38_68.fa" #Reference genome
  region = paste0(chr, ":", format(startpos, scientific = FALSE), "-", format(endpos, scientific = FALSE)) # Region requested

  cmd1 <- paste0(bcftools, " mpileup -q 30 -Ou -r ", region, " -f ", reference, " ", bamstr)
  cmd2 <- paste0(bcftools, " call -mv -Ov ")
  cmd3 <- paste0(bcftools, " view -v snps -i '%QUAL>=30 && DP>10' - -o ", outname, ".snps-filtered.vcf")
  execute(paste0(cmd1, " | ", cmd2, " | ", cmd3))
  invisible("")
}

bamfiles <- c("/home/arends/NAS/Mouse/DNA/Sequencing/Alignment2020/1/SJLP_trimmed.aligned.sorted.dedup.recalibrated.bam",
           "/home/arends/NAS/Mouse/DNA/Sequencing/Alignment2020/2/NZOP_trimmed.aligned.sorted.dedup.recalibrated.bam",
           "/home/arends/NAS/Mouse/DNA/Sequencing/Alignment2020/3/NZOP_trimmed.aligned.sorted.dedup.recalibrated.bam",
           "/home/arends/NAS/Mouse/DNA/Sequencing/Alignment2020/4/KHP_trimmed.aligned.sorted.dedup.recalibrated.bam",
           "/home/arends/NAS/Mouse/DNA/Sequencing/Alignment2020/5/KHP_trimmed.aligned.sorted.dedup.recalibrated.bam",
           "/home/arends/NAS/Mouse/DNA/Sequencing/Alignment2020/6/BFMI861-S1P_trimmed.aligned.sorted.dedup.recalibrated.bam",
           "/home/arends/NAS/Mouse/DNA/Sequencing/Alignment2020/7/BFMI861-S1P_trimmed.aligned.sorted.dedup.recalibrated.bam",
           "/home/arends/NAS/Mouse/DNA/Sequencing/Alignment2020/8/BFMI861-S2P_trimmed.aligned.sorted.dedup.recalibrated.bam",
           "/home/arends/NAS/Mouse/DNA/Sequencing/Alignment2020/9/BFMI861-S2P_trimmed.aligned.sorted.dedup.recalibrated.bam",
           "/home/arends/NAS/Mouse/DNA/Sequencing/Alignment2020/10/BFMI860-12P_trimmed.aligned.sorted.dedup.recalibrated.bam",
           "/home/arends/NAS/Mouse/DNA/Sequencing/Alignment2020/11/BFMI860-12P_trimmed.aligned.sorted.dedup.recalibrated.bam",
           "/home/arends/NAS/Mouse/DNA/Sequencing/Alignment2020/12/BFMI860-12P_trimmed.aligned.sorted.dedup.recalibrated.bam",
           "/home/arends/NAS/Mouse/DNA/Sequencing/Alignment2020/13/BFMI860-12P_trimmed.aligned.sorted.dedup.recalibrated.bam",
           #"/home/arends/NAS/Mouse/DNA/Sequencing/Alignment2020/14/",
           "/home/arends/NAS/Mouse/DNA/Sequencing/Alignment2020/15/B6-3P_trimmed.aligned.sorted.dedup.recalibrated.bam",
           "/home/arends/NAS/Mouse/DNA/Sequencing/Alignment2020/16/B6-4P_trimmed.aligned.sorted.dedup.recalibrated.bam",
           "/home/arends/NAS/Mouse/DNA/Sequencing/Alignment2020/17/B6-5P_trimmed.aligned.sorted.dedup.recalibrated.bam",
           "/home/arends/NAS/Mouse/DNA/Sequencing/Alignment2020/18/B6-5P_trimmed.aligned.sorted.dedup.recalibrated.bam",
           "/home/arends/NAS/Mouse/DNA/Sequencing/Alignment2020/19/B6-5P_trimmed.aligned.sorted.dedup.recalibrated.bam",
           "/home/arends/NAS/Mouse/DNA/Sequencing/Alignment2020/20/AKRP_trimmed.aligned.sorted.dedup.recalibrated.bam",
           #"/home/arends/NAS/Mouse/DNA/Sequencing/Alignment2020/21/"
           #"/home/arends/NAS/Mouse/DNA/Sequencing/Alignment2020/22/"
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/129P2_OlaHsd.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/129S1_SvImJ.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/129S5SvEvBrd.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/A_J.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/AKR_J.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/BALB_cJ.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/BTBR_T__Itpr3tf_J.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/BUB_BnJ.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/C3H_HeH.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/C3H_HeJ.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/C57BL_10J.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/C57BL_6NJ.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/C57BR_cdJ.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/C57L_J.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/C58_J.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/CAST_EiJ.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/CBA_J.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/DBA_1J.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/DBA_2J.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/FVB_NJ.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/I_LnJ.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/JF1_MsJ.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/KK_HiJ.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/LEWES_EiJ.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/LG_J.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/LP_J.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/MOLF_EiJ.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/NOD_ShiLtJ.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/NZB_B1NJ.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/NZO_HlLtJ.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/NZW_LacJ.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/PWK_PhJ.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/RF_J.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/SEA_GnJ.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/SJL_J.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/SM_J.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/SPRET_EiJ.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/ST_bJ.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/WSB_EiJ.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/ZALENDE_EiJ.bam"
           )
		   
regions <- read.table("QTLregions2212020.txt", sep = "\t", header = TRUE)

# Keep the three interesting regions 
regions <- regions[c(44, 23, 47),]

# Get genes in regions
library(biomaRt)

bio.mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")

getregion <- function(bio.mart, chr, startpos, endpos) {
  region <- paste0(chr, ":",startpos, ":", endpos)
  cat("function: ", " has region: ", region, "\n")
  res.biomart <- getBM(attributes = c("ensembl_gene_id",                                                    # Things that we want to get from biomart
                                      "chromosome_name", "start_position", "end_position", "strand", 
                                      "external_gene_name", "mgi_id", "mgi_symbol", "mgi_description"), 
                       filters = c("chromosomal_region", "biotype"),                                        # Things that we will use to query biomart
                       values = list(region, "protein_coding"),                                             # The thing that we are querying
                       mart = bio.mart)

  cat("function: ", " has ", nrow(res.biomart), "\n")
  return (res.biomart) 
}

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

setwd("/home/manuel/AIL_S1xS2/RAWDATA/SNPsGenesGonLiver/SNPsMGP/")

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

# Let´s combine all the SNPs for each gene into one unique vcf file
filelist <- list.files(".") 

allSNPs <- NULL
for(file in filelist){
  if(length(readLines(file)) > 96){
    mcontent <- read.csv(file, sep = "\t", skip = 96, header=FALSE, colClasses=c("character"))
    gene = gsub(".snps-filtered.vcf", "", file, fixed=T)
    allSNPs <- rbind(allSNPs, cbind(gene, mcontent))
  }
}

chromosomes <- c(3, 15, 17)


allSNPsAnnot <- allSNPs[order(allSNPs[,"V2"]),]

annotation <- c()
for (chr in chromosomes){
  annotation <- rbind(annotation, allSNPsAnnot[as.numeric(allSNPsAnnot[,"V1"]) == chr,])
}

allSNPsAnnot <- annotation
header = readLines(filelist[1], n = 96)
cat(paste0(header, collapse = "\n"), "\n", file = "all_combined.vcf")

# File containing all SNPs in the genes
write.table(allSNPsAnnot[,-1], file = "all_combined.vcf", sep = "\t", quote=FALSE, append = TRUE, col.names=FALSE, row.names= FALSE)

# sort the vcf file
cmdsort <- "cat all_combined.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | \"sort -k1,1 -k2,2n\"}' > out_sorted.vcf"

execute(cmdsort)