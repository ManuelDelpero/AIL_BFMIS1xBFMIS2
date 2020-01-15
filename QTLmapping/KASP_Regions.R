# Define regions to genotype with KASP assay for AIL BFMI S1xS2
#
# copyright (c) 2018-2021 - Brockmann group - HU Berlin Manuel Delpero 
# 
# first written december, 2019

#setwd("C:/Users/Manuel/Desktop/AIL_S1xS2/RAWDATA")
setwd("C:/Users/Manuel/Desktop/AIL_S1xS2/RAWDATA")
library(seqinr)

map <- read.table("map.cleaned.txt", sep="\t", check.names = FALSE)
lods <- read.table("lodmatrixADDDOM_nosum.txt", sep = "\t",  header = TRUE)
lods <- t(lods)
chromosomes <- c(1:19, "X", "Y")
map <- map[,c(1,2)]
phenotypes <- rownames(lods)

# Get the main QTLs for each phenotype
res <- c()
for (x in phenotypes){
  Lodscores <- lods[x,]
  ord <- sort(Lodscores, decreasing = TRUE)
  nqtl <- 0
  while ((any(ord > 4)) && (nqtl < 8)){
    info <- cbind(map[names(ord[1]),], names(ord[1]), x, ord[1][[1]])   # get rid of the markers as rownames!!
    markers <- rownames(map[which(map[, "chr"] == info[, "chr"]),])    
    ord <- ord[-which(names(ord) %in% markers)]
    ord <- sort(ord, decreasing = TRUE)
    nqtl <- nqtl + 1
    res <- rbind(res, info)
  }
}
colnames(res) <- c("Chr", "Pos", "Trait", "Lod")

res <- res[-which(is.na(res[, "Chr"] )),]

write.table(res, file = "KASP_regions.txt", sep = "\t", quote = FALSE, row.names = FALSE)
#QTLs <- read.table("KASP_regions.txt", sep = "\t",  header = TRUE)
QTLs <- read.table("QTLs13120.txt", sep = "\t",  header = TRUE)

# Keep just the unique markers
QTLs <- QTLs[-which(duplicated(QTLs[,5])),]
write.table(res, file = "UniqueKASP_regions.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Get the regions for each QTL
regions <- c()
for (x in 1:nrow(QTLs)){
  startpos <- QTLs[, "Pos"] -60
  endpos <- QTLs[, "Pos"] + 60
  chr <- QTLs[, "Chr"]
  regions <- cbind(startpos, endpos, chr)
}
regions <- cbind(regions, as.character(QTLs[,"Marker"]))

colnames(regions) <- c("StartPos", "StopPos", "Chr", "Marker")
rownames(regions) <- QTLs[,"Trait"]
KaspRegions <- regions[c(1,2,7,11,16,23,25,27,28,29,30),]

# Get sequence in the regions
mfasta <- read.fasta("Mus_musculus.GRCm38.dna.toplevel.fa.gz")
for (x in 1:nrow(regions)){
  seq <- mfasta[[regions[x,"Chr"]]][(regions[x,"StartPos"]):(regions[x,"StopPos"])]
  write.table(seq, file = as.character(QTLs[x,"Trait"]), sep="\t", quote = FALSE)
}
  
# SNP calling in the regions
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
  bcftools = "/home/danny/Github/bcftools/bcftools" # Location of BCFtools executable
  reference = "/home/danny/References/Mouse/GRCm38_95/Mus_musculus.GRCm38.dna.toplevel.fa" #Reference genome
  region = paste0(chr, ":", format(startpos, scientific = FALSE), "-", format(endpos, scientific = FALSE)) # Region requested

  cmd1 <- paste0(bcftools, " mpileup -q 30 -Ou -r ", region, " -f ", reference, " ", bamstr)
  cmd2 <- paste0(bcftools, " call -mv -Ov ")
  cmd3 <- paste0(bcftools, " view -v snps -i '%QUAL>=30 && DP>10' - -o ", outname, ".snps-filtered.vcf")
  execute(paste0(cmd1, " | ", cmd2, " | ", cmd3))
  invisible("")
}

bamfiles <- c("/halde/BFMI_Alignment_Mar19/merged_sorted_860-S12.bam",  # 860-S12  (high coverage)
             "/halde/BFMI_Alignment_Mar19/merged_sorted_861-S1.bam",    # 861-S1 (medium coverage)
             "/halde/BFMI_Alignment_Mar19/merged_sorted_861-S2.bam")    # 861-S2 (medium coverage)
			 
			 
name <- c("D21" ,"D35", "D84" , "D150" ,"Gluc172" ,"Gon1"  ,"Gon2"  ,"Longissimus" ,"Length" ,"Triglycerides1" ,"Triglycerides2")
for(x in 1:nrow(KaspRegions)){ 
  startpos <- KaspRegions[x,"StartPos"]
  endpos <- KaspRegions[x,"StopPos"]
  callSNPs(bamfiles, KaspRegions[x,"Chr"], startpos, endpos, name[x]) 
}

