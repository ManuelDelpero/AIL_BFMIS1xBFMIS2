# Analysis in diff. express genes
#
# copyright (c) 2015-2020 - Brockmann group - HU Berlin, Manuel Delpero
# last modified May, 2019
# first written May, 2019

# Get a matrix with gene name, chr, start pos. and end pos.

setwd("/home/manuel/AIL_B6xBFMI/RAWDATA/")
genes <- read.table("gonadalfat_significant_ann.txt", sep = "\t", header = TRUE)
uniquegenes <-  genes[, c("external_gene_name", "chromosome_name", "start_position", "end_position", "strand")] 

# Bam files to use for the SNPs calling function
bamfiles <- c("/halde/BFMI_Alignment_Mar19/merged_sorted_860-S12.bam",  # 860-S12  (high coverage)
             "/halde/BFMI_Alignment_Mar19/merged_sorted_861-S1.bam",    # 861-S1 (medium coverage)
             "/halde/BFMI_Alignment_Mar19/merged_sorted_861-S2.bam")    # 861-S2 (medium coverage)

# Snps in diff. exp. genes
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
