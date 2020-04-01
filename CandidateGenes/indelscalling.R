#
# Use Scalpel to call SNPs from bam files
#

execute <- function(x, intern = FALSE, execute = TRUE){
  cat("----", x, "\n")
  if (execute) {
    res <- system(x, intern = intern)
    cat(">>>>", res[1], "\n")
    if(res[1] >= 1) q("no")
  }
}


callINDELs <- function(bamfiles, chr = 1, bedfile, outname = "myINDELs") {
  bamstr = paste0(bamfiles, collapse = " ") # Collapse all the bam files in a single string
  scalpel = "/home/florian/Downloads/scalpel-0.5.4/scalpel-discovery --single" # Location of scalpel executable
  reference = "/home/danny/References/Mouse/GRCm38_95/Mus_musculus.GRCm38.dna.toplevel.fa" #Reference genome
  region = bedfile # Region requested in bed file
  outdir = # output directory
  cmd <- paste0("nohup", scalpel, " --bam ~", bamstr, " --bed ~", region, " --ref ~", reference, " --window 3500 --mapscore 10 --step 750", " --dir ~", outdir)
  execute(cmd1)
  invisible("")
}

#callSNPs(c("/halde/BFMI_Alignment_Mar19/merged_sorted_860-S12.bam",  # 860-S12  (high coverage)
#           "/halde/BFMI_Alignment_Mar19/merged_sorted_861-S1.bam",   # 861-S1 (medium coverage)
#           "/halde/BFMI_Alignment_Mar19/merged_sorted_861-S2.bam"    # 861-S2 (medium coverage)
#           ), 19, 3000000, 3500000, "test")


# nohup /home/florian/Downloads/scalpel-0.5.4/scalpel-discovery --single\
# --bam ~/NAS/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/129S1_SvImJ.bam\
# --bed ~/scalpelTest/region.bed --ref ~/NAS/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/38_68_ref/GRCm38_68.fa\
# --window 3500 --mapscore 10 --step 750 --dir ~/BBS7_20MB_region_SNPs_Indels/20MB_Indels_all_strains/129S1_SvImJ > 129S1_SvImJ.out