# Plots for AIL BFMI S1xS2
#
# copyright (c) 2018-2021 - Brockmann group - HU Berlin Manuel Delpero 
# 
# first written december, 2019

setwd("C:/Users/Manuel/Desktop/AIL_S1xS2/RAWDATA")

lodmatrixDOM <- read.table("lodmatrixDOM.txt", header = TRUE, sep = "\t", check.names = FALSE)
lodmatrixADD <- read.table("lodmatrixADD.txt", header = TRUE, sep = "\t", check.names = FALSE)
lodmatrixADDDOM <- read.table("lodmatrixADDDOM_nosum.txt", header = TRUE, sep = "\t", check.names = FALSE)
map <- read.table("map.cleaned.txt", sep="\t")

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

# Manhattan plots
chrs <- c(1:19,"X")
gap <- 20000000
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

for (x in chrs){
  chrmid <- chr.lengths[x]/2 + chr.starts[x]
  chrmids <- c(chrmids, chrmid)
}

names(chr.starts) <- chrs
names(chr.lengths) <- chrs
phenotype <- "Gon"
plot(x = c(-gap, tail(chr.starts,1)), y = c(0,8), t = 'n', xlab="", ylab="-log10(P)",xaxt='n', xaxs="i", yaxs="i",las=2,main=paste0("Manhattan plot - ", phenotype, ", Dom + Add effect"))
for(chr in chrs){
  onChr <- rownames(map.sorted[map.sorted[,"Chromosome"] == chr,])
  points(x=chr.starts[chr] + map.sorted[onChr,"Position"], y = mprofiles[onChr, phenotype],t ='p', col=c("black", "cornflowerblue")[(i %% 2 == 1) + 1])
  if(chr != "X") text(x=chr.starts[chr] + chr.lengths[chr] / 2, y = 0.25, chr, col=c("white", "black")[(i %% 2 == 1) + 1],cex=0.8)
  if(chr == "X") text(x=chr.starts[chr] + 2*gap, y = 0.25, chr, col=c("white", "black")[(i %% 2 == 1) + 1],cex=0.8)
  i <- i + 1
}
axis(1, chrs, at = chrmids)
abline(h= 3.8, col="green",lty=3)
abline(h=4.2, col="orange",lty=3)
