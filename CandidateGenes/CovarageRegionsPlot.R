# Creating a coverage plot from Samtools output for S1 and S2 seq data
#
# copyright (c) - Manuel Delpero

setwd("/home/manuel/AIL_S1xS2/CoverageS1S2_Regions/")

#chr 3 coverage
CovS1_3 <- read.table("chr3_S1.coverage", header=FALSE)
covS1_3 <- cbind(CovS1_3[,2], CovS1_3[,3] + CovS1_3[,4])

colnames(covS1_3) <- c("position", "coverage")
chr3_regionS1 <- covS1_3[which((covS1_3[, "position"] > 95763020) & (covS1_3[, "position"] < 100780367)),]

max(chr3_regionS1[,"coverage"])
 
min(chr3_regionS1[,"coverage"])

CovS2_3 <- read.table("chr3_S2.coverage", header=FALSE)
covS2_3 <- cbind(CovS2_3[,2], CovS2_3[,3] + CovS2_3[,4])

colnames(covS2_3) <- c("position", "coverage")
chr3_regionS2 <- covS2_3[which((covS2_3[, "position"] > 95763020) & (covS2_3[, "position"] < 100780367)),]

max(chr3_regionS2[,"coverage"])
 
min(chr3_regionS2[,"coverage"])

# function to plot the coverage across the region, input the subset with position and coverage
PlotCoverage <- function(CoverageRegion1, CoverageRegion2, n_Kb = 10000){
  medCoverages1 <- c()
  medCoverages2 <- c()
  pos1 <- c()
  pos2 <- c()
  for (bp in seq(0,nrow(CoverageRegion1),n_Kb)){
    if (bp != 0){ 
      start <- bp - n_Kb
      end <- bp
    }else{
      start <- 0
      end <- bp
    }
    region <- seq(start,end)
    medCoverage1 <- mean(as.numeric(CoverageRegion1[region, "coverage"]))
    medCoverages1 <- c(medCoverages1, medCoverage1) 
	position1 <- mean(c(start,end))
	pos1 <- c(pos1, position1)
  }
  for (bp in seq(0,nrow(CoverageRegion2),n_Kb)){
    if (bp != 0){ 
      start <- bp - n_Kb
      end <- bp
    }else{
      start <- 0
      end <- bp
    }
    region <- seq(start,end)
    medCoverage2 <- mean(as.numeric(CoverageRegion2[region, "coverage"]))
    medCoverages2 <- c(medCoverages2, medCoverage2) 
	position2 <- mean(c(start,end))
	pos2 <- c(pos2, position2)
  }
  title <- readline(prompt="Enter name of the chromosome: ")  
  plot(c(0,nrow(CoverageRegion1)), c(0,100) , main = title, xlab = "Position", ylab = "Coverage", xaxt = "n", ylim = c(0,100), t = "n")
    points(pos1, medCoverages1,  type = "l", lwd=0.8)
	points(pos2, medCoverages2, type = "l", col = "azure4", lwd=0.8)
    axis(1, at = seq(0,nrow(CoverageRegion1),1000000), as.character(seq(0,nrow(CoverageRegion1),1000000) + as.numeric(CoverageRegion1[1,1]))) 
    legend("topleft",
    legend = c("BFMI-S1", "BFMI-S2"),
    col = c("black", "azure4"),
    pch = c(15,15,15),
    bty = "n",
    pt.cex = 1.8,
    cex = 1.5,)	
}
PlotCoverage(chr3_regionS1, chr3_regionS2, n_Kb = 10000)

# chr 7 coverage
CovS1_7 <- read.table("chr7_S1.coverage", header=FALSE)
covS1_7 <- cbind(CovS1_7[,2], CovS1_7[,3] + CovS1_7[,4])

colnames(covS1_7) <- c("position", "coverage")
chr7_region <- covS1_7[which((covS1_7[, "position"] > 16692359) & (covS1_7[, "position"] < 17022025)),]

max(chr7_regionS1[,"coverage"])
 
min(chr7_regionS1[,"coverage"])

CovS2_7 <- read.table("chr7_S2.coverage", header=FALSE)
covS2_7 <- cbind(CovS2_7[,2], CovS2_7[,3] + CovS2_7[,4])

colnames(covS2_7) <- c("position", "coverage")
chr7_regionS2 <- covS2_7[which((covS2_7[, "position"] > 16692359) & (covS2_7[, "position"] < 17022025)),]

max(chr7_regionS2[,"coverage"])
 
min(chr7_regionS2[,"coverage"])

PlotCoverage(chr7_regionS1, chr7_regionS2, n_Kb = 10000)

# chr 12 coverage
CovS1_12 <- read.table("chr12_S1.coverage", header=FALSE)
covS1_12 <- cbind(CovS1_12[,2], CovS1_12[,3] + CovS1_12[,4])

colnames(covS1_12) <- c("position", "coverage")
chr12_region <- covS1_12[which((covS1_12[, "position"] > 3569173) & (covS1_12[, "position"] < 9870025)),]

max(chr12_region[,"coverage"])
 
min(chr12_region[,"coverage"])

PlotCoverage(chr12_region, n_Kb = 10000)

#chr 15 coverage
CovS1_15 <- read.table("chr15_S1.coverage", header=FALSE)
covS1_15 <- cbind(CovS1_15[,2], CovS1_15[,3] + CovS1_15[,4])

colnames(covS1_15) <- c("position", "coverage")
chr15_region <- covS1_15[which((covS1_15[, "position"] > 3569173) & (covS1_15[, "position"] < 9870025)),]

max(chr15_region[,"coverage"])
 
min(chr15_region[,"coverage"])

PlotCoverage(chr15_region, n_Kb = 10000)

#chr 16 coverage
CovS1_16 <- read.table("chr16_S1.coverage", header=FALSE)
covS1_16 <- cbind(CovS1_16[,2], CovS1_16[,3] + CovS1_16[,4])

colnames(covS1_16) <- c("position", "coverage")
chr16_region <- covS1_16[which((covS1_16[, "position"] > 3569173) & (covS1_16[, "position"] < 9870025)),]

max(chr16_region[,"coverage"])
 
min(chr16_region[,"coverage"])

PlotCoverage(chr16_region, n_Kb = 10000)

#chr17 coverage
CovS1_17 <- read.table("chr17_S1.coverage", header=FALSE)
covS1_17 <- cbind(CovS1_17[,2], CovS1_17[,3] + CovS1_17[,4])

colnames(covS1_17) <- c("position", "coverage")
chr17_region <- covS1_17[which((covS1_17[, "position"] > 9483181) & (covS1_17[, "position"] < 25391933)),]

max(chr17_region[,"coverage"])
 
min(chr17_region[,"coverage"])

PlotCoverage(chr17_region, n_Kb = 10000)
