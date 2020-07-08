#
# Use Scalpel to call INDELs from bam files
#

execute <- function(x, intern = FALSE, execute = TRUE){
  cat("----", x, "\n")
  if (execute) {
    res <- system(x, intern = intern)
    cat(">>>>", res[1], "\n")
    if(res[1] >= 1) q("no")
  }
}


callINDELs <- function(bamfiles) {
  bamstr = paste0(bamfiles, collapse = " ") # Collapse all the bam files in a single string
  scalpel = "/home/florian/Downloads/scalpel-0.5.4/scalpel-discovery --single" # Location of scalpel executable
  reference = "/home/danny/References/Mouse/GRCm38_95/Mus_musculus.GRCm38.dna.toplevel.fa" #Reference genome
  bedfile = "/home/manuel/AIL_S1xS2/RAWDATA/INDELsGenesGonLiverS1/bedfile.bed.txt"
  cmd <- paste0(scalpel, " --bam ", bamstr, " --bed ", bedfile, " --ref ", reference)
  execute(cmd)
  invisible("")
}




