# Plots for AIL BFMI S1xS2
#
# copyright (c) 2018-2021 - Brockmann group - HU Berlin Manuel Delpero 
# 
# first written december, 2019

setwd("C:/Users/Manuel/Desktop/AIL_S1xS2/RAWDATA")

lodmatrixDOM <- read.table("lodmatrixDOM.txt", header = TRUE, sep = "\t", check.names = FALSE)
lodmatrixADD <- read.table("lodmatrixADD.txt", header = TRUE, sep = "\t", check.names = FALSE)
lodmatrixADDDOM <- read.table("lodmatrixADDDOM_nosum.txt", header = TRUE, sep = "\t", check.names = FALSE)
markerannot <- read.csv("map.cleaned.txt", header=TRUE, sep="\t", check.names=FALSE)
markerannot <- markerannot[, c(1,2)]
markerannot$Index <- seq.int(nrow(markerannot))
colnames(markerannot) <- c("Chromosome", "Position", "Index")
markerannot <- markerannot[-which(is.na(markerannot[,2])),]

chromosomes <- c(1:19, "X", "Y")

annotation <- c()
for(chr in chromosomes){
  annotation <- rbind(annotation, map[map[,"chr"] == chr,])
}

lodmatrixDOM <- lodmatrixDOM[rownames(annotation),]
lodmatrixADD <- lodmatrixADD[rownames(annotation),]
lodmatrixADDDOM <- lodmatrixADDDOM[rownames(annotation),]

# Preliminary visual check bodyweight 
bw <- lodmatrixADDDOM[,c(1:32)]
rotate <- function(x) t(apply(x, 2, rev))
image(rotate(bw))
lodannotmatrix <- cbind(annotation[rownames(lodmatrixADDDOM), ], lodmatrixADDDOM)
chr17 <- lodannotmatrix[which(lodannotmatrix[,"chr"] == 17),]
chr17ord <- chr17[order(chr17[,2], decreasing = FALSE),]
plot(main = "QTL Gon weight [Chr 17]", c(min(as.numeric(chr17ord[, "bp_mm10"])), max(as.numeric(chr17ord[, "bp_mm10"]))), c(0,9), ylab = "-log10 [pvalue]", xlab = "Position [mb]", las = 2, t = "n", xaxt = "n")
points(x = as.numeric(chr17ord[,"bp_mm10"]), y = chr17ord[,"Gon"] , type = "l", col="dodgerblue", lwd = 1)

## Manhattan plots
chrs <- c(1:19,"X")
gap <- 40000000
map.sorted <- NULL
chr.lengths <- c()
chr.starts <- c(0)
chrmids <- c()
i <- 1
for(chr in chrs){
  onChr <- which(markerannot[,"Chromosome"] == chr)
  map.sorted <- rbind(map.sorted, markerannot[onChr,])
  chr.lengths <- c(chr.lengths, max(markerannot[onChr, "Position"]))
  chr.starts <- c(chr.starts, chr.starts[i] + max(markerannot[onChr, "Position"]) + gap)
  i <- i + 1
}

names(chr.starts) <- chrs
names(chr.lengths) <- chrs

for (x in chrs){
  chrmid <- as.numeric(chr.lengths[x]/2) + as.numeric(chr.starts[x])
  chrmids <- c(chrmids, chrmid)
}

# Gonadal adipose tissue weight
phenotype <- "Gon"
plot(x = c(-gap, tail(chr.starts,1)), y = c(0,8), t = 'n', xlab="Chromosome", ylab="-log10(P)",xaxt='n', xaxs="i", yaxs="i",las=2,main=paste0("Manhattan plot - ", phenotype, ", Dom + Add effect"))
for(chr in chrs){
  onChr <- rownames(map.sorted[map.sorted[,"Chromosome"] == chr,])
  points(x=chr.starts[chr] + map.sorted[onChr,"Position"], y = mprofiles[onChr, phenotype],t ='p', pch = 16, col=c("black", "cornflowerblue")[(i %% 2 == 1) + 1])
  i <- i + 1
}
axis(1, chrs, at = chrmids)
abline(h= 3.8, col="green",lty=3)
abline(h=4.2, col="orange",lty=3)

# liver weight
phenotype <- "Leber"
plot(x = c(-gap, tail(chr.starts,1)), y = c(0,8), t = 'n', xlab="Chromosome", ylab="-log10(P)",xaxt='n', xaxs="i", yaxs="i",las=2,main=paste0("Manhattan plot - ", phenotype, ", Dom + Add effect"))
for(chr in chrs){
  onChr <- rownames(map.sorted[map.sorted[,"Chromosome"] == chr,])
  points(x=chr.starts[chr] + map.sorted[onChr,"Position"], y = mprofiles[onChr, phenotype],t ='p', pch = 16, col=c("black", "cornflowerblue")[(i %% 2 == 1) + 1])
  i <- i+1
}
axis(1, chrs, at = chrmids)
abline(h= 3.8, col="green",lty=3)
abline(h=4.2, col="orange",lty=3)

# Triglycerides
phenotype <- "Triglycerides"
plot(x = c(-gap, tail(chr.starts,1)), y = c(0,8), t = 'n', xlab="", ylab="-log10(P)",xaxt='n', xaxs="i", yaxs="i",las=2,main=paste0("Manhattan plot - ", phenotype, ", Dom + Add effect"))
for(chr in chrs){
  onChr <- rownames(map.sorted[map.sorted[,"Chromosome"] == chr,])
  points(x=chr.starts[chr] + map.sorted[onChr,"Position"], y = mprofiles[onChr, phenotype],t ='p', pch = 16,  col=c("black", "cornflowerblue")[(i %% 2 == 1) + 1])
  i <- i + 1
}
axis(1, chrs, at = chrmids)
abline(h= 3.8, col="green",lty=3)
abline(h=4.2, col="orange",lty=3)

# D140
phenotype <- "D140"
plot(x = c(-gap, tail(chr.starts,1)), y = c(0,8), t = 'n', xlab="Chromosome", ylab="-log10(P)",xaxt='n', xaxs="i", yaxs="i",las=2,main=paste0("Manhattan plot - ", phenotype, ", Dom + Add effect"))
for(chr in chrs){
  onChr <- rownames(map.sorted[map.sorted[,"Chromosome"] == chr,])
  points(x=chr.starts[chr] + map.sorted[onChr,"Position"], y = mprofiles[onChr, phenotype],t ='p', pch = 16, col=c("black", "cornflowerblue")[(i %% 2 == 1) + 1])
  i <- i + 1
}
axis(1, chrs, at = chrmids)
abline(h= 3.8, col="green",lty=3)
abline(h=4.2, col="orange",lty=3)

# Gluc172
phenotype <- "Gluc172"
plot(x = c(-gap, tail(chr.starts,1)), y = c(0,8), t = 'n', xlab="Chromosome", ylab="-log10(P)",xaxt='n', xaxs="i", yaxs="i",las=2,main=paste0("Manhattan plot - ", phenotype, ", Dom + Add effect"))
for(chr in chrs){
  onChr <- rownames(map.sorted[map.sorted[,"Chromosome"] == chr,])
  points(x=chr.starts[chr] + map.sorted[onChr,"Position"], y = mprofiles[onChr, phenotype],t ='p', pch = 16, col=c("black", "cornflowerblue")[(i %% 2 == 1) + 1])
  i <- i + 1
}
axis(1, chrs, at = chrmids)
abline(h= 3.8, col="green",lty=3)
abline(h=4.2, col="orange",lty=3)

# Gluc172
phenotype <- "Gluc172"
plot(x = c(-gap, tail(chr.starts,1)), y = c(0,8), t = 'n', xlab="Chromosome", ylab="-log10(P)",xaxt='n', xaxs="i", yaxs="i",las=2,main=paste0("Manhattan plot - ", phenotype, ", Dom + Add effect"))
for(chr in chrs){
  onChr <- rownames(map.sorted[map.sorted[,"Chromosome"] == chr,])
  points(x=chr.starts[chr] + map.sorted[onChr,"Position"], y = mprofiles[onChr, phenotype],t ='p', pch = 16, col=c("black", "cornflowerblue")[(i %% 2 == 1) + 1])
  i <- i + 1
}
axis(1, chrs, at = chrmids)
abline(h= 3.8, col="green",lty=3)
abline(h=4.2, col="orange",lty=3)

# ITTauc
phenotype <- "ITTauc"
plot(x = c(-gap, tail(chr.starts,1)), y = c(0,8), t = 'n', xlab="Chromosome", ylab="-log10(P)",xaxt='n', xaxs="i", yaxs="i",las=2,main=paste0("Manhattan plot - ", phenotype, ", Dom + Add effect"))
for(chr in chrs){
  onChr <- rownames(map.sorted[map.sorted[,"Chromosome"] == chr,])
  points(x=chr.starts[chr] + map.sorted[onChr,"Position"], y = mprofiles[onChr, phenotype],t ='p', pch = 16, col=c("black", "cornflowerblue")[(i %% 2 == 1) + 1])
  i <- i + 1
}
axis(1, chrs, at = chrmids)
abline(h= 3.8, col="green",lty=3)
abline(h=4.2, col="orange",lty=3)

# length
phenotype <- "LÃ.nge"
plot(x = c(-gap, tail(chr.starts,1)), y = c(0,8), t = 'n', xlab="Chromosome", ylab="-log10(P)",xaxt='n', xaxs="i", yaxs="i",las=2,main=paste0("Manhattan plot - ", phenotype, ", Dom + Add effect"))
for(chr in chrs){
  onChr <- rownames(map.sorted[map.sorted[,"Chromosome"] == chr,])
  points(x=chr.starts[chr] + map.sorted[onChr,"Position"], y = mprofiles[onChr, phenotype],t ='p', pch = 16, col=c("black", "cornflowerblue")[(i %% 2 == 1) + 1])
  i <- i + 1
}
axis(1, chrs, at = chrmids)
abline(h= 3.8, col="green",lty=3)
abline(h=4.2, col="orange",lty=3)