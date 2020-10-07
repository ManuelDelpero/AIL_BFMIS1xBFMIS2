# Creating a coverage plot from Samtools output for S1 and S2 seq data
#
# copyright (c) - Manuel Delpero

setwd("/home/manuel/AIL_S1xS2/CoverageS1S2_Regions/")

#chr 3 coverage
CovS1_3 <- read.table("chr3_S1.coverage", header=FALSE)
covS1_3 <- cbind(CovS1_3[,2], CovS1_3[,3] + CovS1_3[,4])

colnames(covS1_3) <- c("position", "coverage")
chr3_region <- covS1_3[which((covS1_3[, "position"] > 95763020) & (covS1_3[, "position"] < 100780367)),]

max(chr3_region[,"coverage"])
 
min(chr3_region[,"coverage"])

PlotCoverage <- function(CoverageRegion, n_Kb = 10000){
  medCoverages <- c()
  for (bp in seq(0,nrow(CoverageRegion),n_Kb)){
    if (bp != 0){ 
      start <- bp - n_Kb
      end <- bp
    }else{
      start <- 0
      end <- bp
    }
    region <- seq(start,end)
    medCoverage <- mean(as.numeric(CoverageRegion[region, "coverage"]))
    medCoverages <- c(medCoverages, medCoverage) 
  } 
plot(seq(0,nrow(CoverageRegion),n_Kb), medCoverages, type = "l", lwd=0.5, main = "Coverage chromosome 3 region", xlab = "Position", ylab = "Coverage")
}
PlotCoverage(chr3_region, n_Kb = 10000)