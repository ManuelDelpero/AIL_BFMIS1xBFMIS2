# Various plots to represent the results for AIL BFMI S1xS2
#
# copyright (c) 2018-2021 - Brockmann group - HU Berlin Manuel Delpero 
# 
# first written december, 2019

setwd("C:/Users/Manuel/Desktop/AIL_S1xS2/RAWDATA")

phenotypes <- read.csv("allPhenotypes.txt", header = TRUE, check.names = FALSE, sep = "\t", colClasses = "character")
genotypes <- read.csv("genomatrix.clean.txt", header = TRUE, check.names = FALSE, sep="\t", colClasses="character")
lodmatrixDOM <- read.table("lodmatrixDOM.txt", header = TRUE, sep = "\t", check.names = FALSE)
lodmatrixADD <- read.table("lodmatrixADD.txt", header = TRUE, sep = "\t", check.names = FALSE)
mprofiles <- read.table("lodmatrixADDDOM_nosum.txt", header = TRUE, sep = "\t", check.names = FALSE)
markerannot <- read.csv("map.cleaned.txt", header=TRUE, sep="\t", check.names=FALSE)
markerannot <- markerannot[, c(1,2)]
markerannot$Index <- seq.int(nrow(markerannot))
colnames(markerannot) <- c("Chromosome", "Position", "Index")
markerannot <- markerannot[-which(is.na(markerannot[,2])),]

chromosomes <- c(1:19, "X", "Y")

annotation <- c()
for(chr in chromosomes){
  annotation <- rbind(annotation, markerannot[markerannot[,"Chromosome"] == chr,])
}

lodmatrixDOM <- lodmatrixDOM[rownames(annotation),]
lodmatrixADD <- lodmatrixADD[rownames(annotation),]
lodmatrixADDDOM <- mprofiles[rownames(annotation),]

# Preliminary visual check bodyweight 
bw <- mprofiles[,c(1:32)]
rotate <- function(x) t(apply(x, 2, rev))
image(rotate(bw))
lodannotmatrix <- cbind(annotation[rownames(lodmatrixADDDOM), ], lodmatrixADDDOM)
chr17 <- lodannotmatrix[which(lodannotmatrix[,"Chromosome"] == 17),]
chr17ord <- chr17[order(chr17[,2], decreasing = FALSE),]
# Gon chr 17
plot(main = "QTL Gon weight [Chr 17]", c(min(as.numeric(chr17ord[, "Position"])), max(as.numeric(chr17ord[, "Position"]))), c(0,9), ylab = "-log10 [pvalue]", xlab = "Position [mb]", las = 2, t = "n", xaxt = "n")
points(x = as.numeric(chr17ord[,"Position"]), y = chr17ord[,"Gon"] , type = "l", col="dodgerblue", lwd = 1)

# Liver chr 17
plot(main = "QTL liver weight [Chr 17]", c(min(as.numeric(chr17ord[, "Position"])), max(as.numeric(chr17ord[, "Position"]))), c(0,9), ylab = "-log10 [pvalue]", xlab = "Position [mb]", las = 2, t = "n", xaxt = "n")
points(x = as.numeric(chr17ord[,"Position"]), y = chr17ord[,"Leber"] , type = "l", col="dodgerblue", lwd = 1)

# Body weight different time points for chr 15 (main QTL)
lodannotmatrix <- cbind(annotation[rownames(lodmatrixADDDOM), ], lodmatrixADDDOM)
dataset <- lodannotmatrix[, c("Chromosome", "Position", "D28", "D49", "D35", "D42", "D49", "D56", "D63", "D70", "D77", "D84", "D91", "D98", "D105", "D112", "D119", "D126", "D133", "D140", "D160", "D172")]
chr15 <- dataset[which(dataset[,"Chromosome"] == 15),]
chr15 <- chr15[order(chr15[,"Position"]),]
plot(main = "QTL profile bodyweight [Chr 15]", c(min(as.numeric(chr15[, "Position"])), max(as.numeric(chr15[, "Position"]))), c(0,9), ylab = "-log10 [pvalue]", xlab = "Position [mb]", las = 2, t = "n", xaxt = "n")
  points(x = as.numeric(chr15[,"Position"]), y = chr15[,"D70"] , type = "l", col="dodgerblue", lwd = 1)
  points(x = as.numeric(chr15[,"Position"]), y = chr15[,"D140"] , type = "l", col="blue", lwd = 1)
  points(x = as.numeric(chr15[,"Position"]), y = chr15[,"D63"] , type = "l", col="deepskyblue", lwd = 1)
  points(x = as.numeric(chr15[,"Position"]), y = chr15[,"D105"] , type = "l", col="purple", lwd = 1)
  points(x = as.numeric(chr15[,"Position"]), y = chr15[,"D98"] , type = "l", col="dodgerblue4", lwd = 1)
  abline(h=4.5, col="green")
  abline(h=4, col="orange")
  axis(1, at = c(0,25000000, 50000000, 75000000, 100000000), c("0", "25", "50", "75", "100"))
  legend("topright", bg="gray",
  legend = c("Week 9", "Week 10", "Week 14", "Week 15","Week 20"),
    col = c("deepskyblue", "dodgerblue", "dodgerblue4", "purple", "blue"),
    pch = 15,
    pt.cex = 1.7,
    pt.bg = "lightsteelblue1",
    cex = 1,
    text.col = "black")

## Manhattan plots (Plot the effect with the highest  lod score and use three different colors for each one) 
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
plot(x = c(-gap, tail(chr.starts,1)), y = c(0,8), t = 'n', xlab="Chromosome", ylab="-log10(P)",xaxt='n', xaxs="i", yaxs="i",las=2,main=paste0("Manhattan plot - ", phenotype))
for(chr in chrs){
  onChr <- rownames(map.sorted[map.sorted[,"Chromosome"] == chr,])
  currentADDDOM <- mprofiles[onChr, phenotype]
  currentDOM <- lodmatrixDOM[onChr, phenotype]
  currentADD <- lodmatrixADD[onChr, phenotype]
  if (chr == "X")
    points(x=chr.starts[chr] + map.sorted[onChr,"Position"], y = mprofiles[onChr, phenotype], t ='p', pch = 16, col= "black")
  for (p in 1:length(currentADDDOM)){
    pos <- chr.starts[chr] + map.sorted[onChr,"Position"]
    if ((currentADDDOM[p] >  currentDOM[p]) && (currentADDDOM[p] > currentADD[p])){
      if (chr %in% seq(1,20,2))
        points(x=pos[p], y = currentADDDOM[p], t ='p', pch = 16, col= "cornflowerblue")
      else (points(x=pos[p], y = currentADDDOM[p], t ='p', pch = 16, col= "black"))
    }
    if ((currentDOM[p] >  currentADDDOM[p]) && (currentDOM[p] > currentADD[p])){
      if (chr %in% seq(1,20,2))
        points(x=pos[p], y = currentDOM[p], t ='p', pch = 17, col= "cornflowerblue")
      else (points(x=pos[p], y = currentDOM[p], t ='p', pch = 17, col= "black"))
    }
      if ((currentADD[p] >  currentADDDOM[p]) && (currentADD[p] > currentDOM[p])){
        if (chr %in% seq(1,20,2))
          points(x=pos[p], y = currentADD[p], t ='p', pch = 18, col= "cornflowerblue")
	else (points(x=pos[p], y = currentADD[p], t ='p', pch = 18, col= "black"))
    }
  }
}
axis(1, chrs, at = chrmids)
abline(h= 3.8, col="green",lty=3)
abline(h=4.2, col="orange",lty=3)
axis(1, chrs, at = chrmids)
abline(h= 3.8, col="green",lty=3)
abline(h=4.2, col="orange",lty=3)
legend("topright", bg="gray",
  legend = c("ADD + DOM dev.", "ADD", "DOM dev."),
  pch = c(16, 18, 17))
  
# Effect plot for the top marker on chr17
colnames(genotypes) <- gsub("V 888-", "", colnames(genotypes))
phenotypes <- phenotypes[colnames(genotypes),]
boxplot(as.numeric(phenotypes[, "Gon"]) ~ unlist(genotypes["JAX00432128",]), col = c("lightskyblue1", "cyan3", "dodgerblue4"), main = "Effect plot for marker JAX00432128 [Gonadal weight]", ylab = "Weight [Gr.]", xlab = "Genotypes")
  legend("topright", bg="gray",
  legend = c( "BFMI-S2", "HET", "BFMI-S1"),
    col = c("lightskyblue1", "cyan3", "dodgerblue4"),
    pch = 15,
    pt.cex = 1.7,
    pt.bg = "lightsteelblue1",
    cex = 1,
    text.col = "black")
  

# liver weight
phenotype <- "Leber"
plot(x = c(-gap, tail(chr.starts,1)), y = c(0,8), t = 'n', xlab="Chromosome", ylab="-log10(P)",xaxt='n', xaxs="i", yaxs="i",las=2,main=paste0("Manhattan plot - ", phenotype))
for(chr in chrs){
  onChr <- rownames(map.sorted[map.sorted[,"Chromosome"] == chr,])
  currentADDDOM <- mprofiles[onChr, phenotype]
  currentDOM <- lodmatrixDOM[onChr, phenotype]
  currentADD <- lodmatrixADD[onChr, phenotype]
  if (chr == "X")
    points(x=chr.starts[chr] + map.sorted[onChr,"Position"], y = mprofiles[onChr, phenotype], t ='p', pch = 16, col= "black")
  for (p in 1:length(currentADDDOM)){
    pos <- chr.starts[chr] + map.sorted[onChr,"Position"]
    if ((currentADDDOM[p] >  currentDOM[p]) && (currentADDDOM[p] > currentADD[p])){
      if (chr %in% seq(1,20,2))
        points(x=pos[p], y = currentADDDOM[p], t ='p', pch = 16, col= "cornflowerblue")
      else (points(x=pos[p], y = currentADDDOM[p], t ='p', pch = 16, col= "black"))
    }
    if ((currentDOM[p] >  currentADDDOM[p]) && (currentDOM[p] > currentADD[p])){
      if (chr %in% seq(1,20,2))
        points(x=pos[p], y = currentDOM[p], t ='p', pch = 17, col= "cornflowerblue")
      else (points(x=pos[p], y = currentDOM[p], t ='p', pch = 17, col= "black"))
    }
      if ((currentADD[p] >  currentADDDOM[p]) && (currentADD[p] > currentDOM[p])){
        if (chr %in% seq(1,20,2))
          points(x=pos[p], y = currentADD[p], t ='p', pch = 18, col= "cornflowerblue")
	else (points(x=pos[p], y = currentADD[p], t ='p', pch = 18, col= "black"))
    }
  }
}
axis(1, chrs, at = chrmids)
abline(h= 3.8, col="green",lty=3)
abline(h=4.2, col="orange",lty=3)
axis(1, chrs, at = chrmids)
abline(h= 3.8, col="green",lty=3)
abline(h=4.2, col="orange",lty=3)
legend("topright", bg="gray",
  legend = c("ADD + DOM dev.", "ADD", "DOM dev."),
  pch = c(16, 18, 17))

# Effect plot for top marker for Gon weight using liver weight as phenotypes, there is a switch!!  responsile for ectopic fat storage in the liver
boxplot(as.numeric(phenotypes[, "Leber"]) ~ unlist(genotypes["JAX00432128",]), col = c("lightskyblue1", "cyan3", "dodgerblue4"), main = "Effect plot for marker JAX00432128 [liver weight]", ylab = "Weight [Gr.]", xlab = "Genotypes")
  legend("topright", bg="gray",
  legend = c( "BFMI-S2", "HET", "BFMI-S1"),
    col = c("lightskyblue1", "cyan3", "dodgerblue4"),
    pch = 15,
    pt.cex = 1.7,
    pt.bg = "lightsteelblue1",
    cex = 1,
    text.col = "black")

# Triglycerides
phenotype <- "Triglycerides"
plot(x = c(-gap, tail(chr.starts,1)), y = c(0,8), t = 'n', xlab="Chromosome", ylab="-log10(P)",xaxt='n', xaxs="i", yaxs="i",las=2,main=paste0("Manhattan plot - ", phenotype))
for(chr in chrs){
  onChr <- rownames(map.sorted[map.sorted[,"Chromosome"] == chr,])
  currentADDDOM <- mprofiles[onChr, phenotype]
  currentDOM <- lodmatrixDOM[onChr, phenotype]
  currentADD <- lodmatrixADD[onChr, phenotype]
  if (chr == "X")
    points(x=chr.starts[chr] + map.sorted[onChr,"Position"], y = mprofiles[onChr, phenotype], t ='p', pch = 16, col= "black")
  for (p in 1:length(currentADDDOM)){
    pos <- chr.starts[chr] + map.sorted[onChr,"Position"]
    if ((currentADDDOM[p] >  currentDOM[p]) && (currentADDDOM[p] > currentADD[p])){
      if (chr %in% seq(1,20,2))
        points(x=pos[p], y = currentADDDOM[p], t ='p', pch = 16, col= "cornflowerblue")
      else (points(x=pos[p], y = currentADDDOM[p], t ='p', pch = 16, col= "black"))
    }
    if ((currentDOM[p] >  currentADDDOM[p]) && (currentDOM[p] > currentADD[p])){
      if (chr %in% seq(1,20,2))
        points(x=pos[p], y = currentDOM[p], t ='p', pch = 17, col= "cornflowerblue")
      else (points(x=pos[p], y = currentDOM[p], t ='p', pch = 17, col= "black"))
    }
      if ((currentADD[p] >  currentADDDOM[p]) && (currentADD[p] > currentDOM[p])){
        if (chr %in% seq(1,20,2))
          points(x=pos[p], y = currentADD[p], t ='p', pch = 18, col= "cornflowerblue")
	else (points(x=pos[p], y = currentADD[p], t ='p', pch = 18, col= "black"))
    }
  }
}
axis(1, chrs, at = chrmids)
abline(h= 3.8, col="green",lty=3)
abline(h=4.2, col="orange",lty=3)
axis(1, chrs, at = chrmids)
abline(h= 3.8, col="green",lty=3)
abline(h=4.2, col="orange",lty=3)
legend("topright", bg="gray",
  legend = c("ADD + DOM dev.", "ADD", "DOM dev."),
  pch = c(16, 18, 17))

# D140
phenotype <- "D140"
plot(x = c(-gap, tail(chr.starts,1)), y = c(0,8), t = 'n', xlab="Chromosome", ylab="-log10(P)",xaxt='n', xaxs="i", yaxs="i",las=2,main=paste0("Manhattan plot - ", phenotype))
for(chr in chrs){
  onChr <- rownames(map.sorted[map.sorted[,"Chromosome"] == chr,])
  currentADDDOM <- mprofiles[onChr, phenotype]
  currentDOM <- lodmatrixDOM[onChr, phenotype]
  currentADD <- lodmatrixADD[onChr, phenotype]
  if (chr == "X")
    points(x=chr.starts[chr] + map.sorted[onChr,"Position"], y = mprofiles[onChr, phenotype], t ='p', pch = 16, col= "black")
  for (p in 1:length(currentADDDOM)){
    pos <- chr.starts[chr] + map.sorted[onChr,"Position"]
    if ((currentADDDOM[p] >  currentDOM[p]) && (currentADDDOM[p] > currentADD[p])){
      if (chr %in% seq(1,20,2))
        points(x=pos[p], y = currentADDDOM[p], t ='p', pch = 16, col= "cornflowerblue")
      else (points(x=pos[p], y = currentADDDOM[p], t ='p', pch = 16, col= "black"))
    }
    if ((currentDOM[p] >  currentADDDOM[p]) && (currentDOM[p] > currentADD[p])){
      if (chr %in% seq(1,20,2))
        points(x=pos[p], y = currentDOM[p], t ='p', pch = 17, col= "cornflowerblue")
      else (points(x=pos[p], y = currentDOM[p], t ='p', pch = 17, col= "black"))
    }
      if ((currentADD[p] >  currentADDDOM[p]) && (currentADD[p] > currentDOM[p])){
        if (chr %in% seq(1,20,2))
          points(x=pos[p], y = currentADD[p], t ='p', pch = 18, col= "cornflowerblue")
	else (points(x=pos[p], y = currentADD[p], t ='p', pch = 18, col= "black"))
    }
  }
}
axis(1, chrs, at = chrmids)
abline(h= 3.8, col="green",lty=3)
abline(h=4.2, col="orange",lty=3)
axis(1, chrs, at = chrmids)
abline(h= 3.8, col="green",lty=3)
abline(h=4.2, col="orange",lty=3)
legend("topright", bg="gray",
  legend = c("ADD + DOM dev.", "ADD", "DOM dev."),
  pch = c(16, 18, 17))

# Gluc172
phenotype <- "Gluc172"
plot(x = c(-gap, tail(chr.starts,1)), y = c(0,8), t = 'n', xlab="Chromosome", ylab="-log10(P)",xaxt='n', xaxs="i", yaxs="i",las=2,main=paste0("Manhattan plot - ", phenotype))
for(chr in chrs){
  onChr <- rownames(map.sorted[map.sorted[,"Chromosome"] == chr,])
  currentADDDOM <- mprofiles[onChr, phenotype]
  currentDOM <- lodmatrixDOM[onChr, phenotype]
  currentADD <- lodmatrixADD[onChr, phenotype]
  if (chr == "X")
    points(x=chr.starts[chr] + map.sorted[onChr,"Position"], y = mprofiles[onChr, phenotype], t ='p', pch = 16, col= "black")
  for (p in 1:length(currentADDDOM)){
    pos <- chr.starts[chr] + map.sorted[onChr,"Position"]
    if ((currentADDDOM[p] >  currentDOM[p]) && (currentADDDOM[p] > currentADD[p])){
      if (chr %in% seq(1,20,2))
        points(x=pos[p], y = currentADDDOM[p], t ='p', pch = 16, col= "cornflowerblue")
      else (points(x=pos[p], y = currentADDDOM[p], t ='p', pch = 16, col= "black"))
    }
    if ((currentDOM[p] >  currentADDDOM[p]) && (currentDOM[p] > currentADD[p])){
      if (chr %in% seq(1,20,2))
        points(x=pos[p], y = currentDOM[p], t ='p', pch = 17, col= "cornflowerblue")
      else (points(x=pos[p], y = currentDOM[p], t ='p', pch = 17, col= "black"))
    }
      if ((currentADD[p] >  currentADDDOM[p]) && (currentADD[p] > currentDOM[p])){
        if (chr %in% seq(1,20,2))
          points(x=pos[p], y = currentADD[p], t ='p', pch = 18, col= "cornflowerblue")
	else (points(x=pos[p], y = currentADD[p], t ='p', pch = 18, col= "black"))
    }
  }
}
axis(1, chrs, at = chrmids)
abline(h= 3.8, col="green",lty=3)
abline(h=4.2, col="orange",lty=3)
axis(1, chrs, at = chrmids)
abline(h= 3.8, col="green",lty=3)
abline(h=4.2, col="orange",lty=3)
legend("topright", bg="gray",
  legend = c("ADD + DOM dev.", "ADD", "DOM dev."),
  pch = c(16, 18, 17))


# Effect plot for top marker for Gon weight using final glucose as phenotypes, this marker is responsile for ectopic fat storage and also for the glucose level 
boxplot(as.numeric(phenotypes[, "Gluc172"]) ~ unlist(genotypes["JAX00432128",]), col = c("lightskyblue1", "cyan3", "dodgerblue4"), main = "Effect plot for marker JAX00432128 [Final glucose]", ylab = "Weight [Gr.]", xlab = "Genotypes")
  legend("topright", bg="gray",
  legend = c( "BFMI-S2", "HET", "BFMI-S1"),
    col = c("lightskyblue1", "cyan3", "dodgerblue4"),
    pch = 15,
    pt.cex = 1.7,
    pt.bg = "lightsteelblue1",
    cex = 1,
    text.col = "black")

# ITTauc
phenotype <- "ITTauc"
plot(x = c(-gap, tail(chr.starts,1)), y = c(0,8), t = 'n', xlab="Chromosome", ylab="-log10(P)",xaxt='n', xaxs="i", yaxs="i",las=2,main=paste0("Manhattan plot - ", phenotype))
for(chr in chrs){
  onChr <- rownames(map.sorted[map.sorted[,"Chromosome"] == chr,])
  currentADDDOM <- mprofiles[onChr, phenotype]
  currentDOM <- lodmatrixDOM[onChr, phenotype]
  currentADD <- lodmatrixADD[onChr, phenotype]
  if (chr == "X")
    points(x=chr.starts[chr] + map.sorted[onChr,"Position"], y = mprofiles[onChr, phenotype], t ='p', pch = 16, col= "black")
  for (p in 1:length(currentADDDOM)){
    pos <- chr.starts[chr] + map.sorted[onChr,"Position"]
    if ((currentADDDOM[p] >  currentDOM[p]) && (currentADDDOM[p] > currentADD[p])){
      if (chr %in% seq(1,20,2))
        points(x=pos[p], y = currentADDDOM[p], t ='p', pch = 16, col= "cornflowerblue")
      else (points(x=pos[p], y = currentADDDOM[p], t ='p', pch = 16, col= "black"))
    }
    if ((currentDOM[p] >  currentADDDOM[p]) && (currentDOM[p] > currentADD[p])){
      if (chr %in% seq(1,20,2))
        points(x=pos[p], y = currentDOM[p], t ='p', pch = 17, col= "cornflowerblue")
      else (points(x=pos[p], y = currentDOM[p], t ='p', pch = 17, col= "black"))
    }
      if ((currentADD[p] >  currentADDDOM[p]) && (currentADD[p] > currentDOM[p])){
        if (chr %in% seq(1,20,2))
          points(x=pos[p], y = currentADD[p], t ='p', pch = 18, col= "cornflowerblue")
	else (points(x=pos[p], y = currentADD[p], t ='p', pch = 18, col= "black"))
    }
  }
}
axis(1, chrs, at = chrmids)
abline(h= 3.8, col="green",lty=3)
abline(h=4.2, col="orange",lty=3)
axis(1, chrs, at = chrmids)
abline(h= 3.8, col="green",lty=3)
abline(h=4.2, col="orange",lty=3)
legend("topright", bg="gray",
  legend = c("ADD + DOM dev.", "ADD", "DOM dev."),
  pch = c(16, 18, 17))

# length
phenotype <- "LÃ.nge"
plot(x = c(-gap, tail(chr.starts,1)), y = c(0,8), t = 'n', xlab="Chromosome", ylab="-log10(P)",xaxt='n', xaxs="i", yaxs="i",las=2,main=paste0("Manhattan plot - ", phenotype))
for(chr in chrs){
  onChr <- rownames(map.sorted[map.sorted[,"Chromosome"] == chr,])
  currentADDDOM <- mprofiles[onChr, phenotype]
  currentDOM <- lodmatrixDOM[onChr, phenotype]
  currentADD <- lodmatrixADD[onChr, phenotype]
  if (chr == "X")
    points(x=chr.starts[chr] + map.sorted[onChr,"Position"], y = mprofiles[onChr, phenotype], t ='p', pch = 16, col= "black")
  for (p in 1:length(currentADDDOM)){
    pos <- chr.starts[chr] + map.sorted[onChr,"Position"]
    if ((currentADDDOM[p] >  currentDOM[p]) && (currentADDDOM[p] > currentADD[p])){
      if (chr %in% seq(1,20,2))
        points(x=pos[p], y = currentADDDOM[p], t ='p', pch = 16, col= "cornflowerblue")
      else (points(x=pos[p], y = currentADDDOM[p], t ='p', pch = 16, col= "black"))
    }
    if ((currentDOM[p] >  currentADDDOM[p]) && (currentDOM[p] > currentADD[p])){
      if (chr %in% seq(1,20,2))
        points(x=pos[p], y = currentDOM[p], t ='p', pch = 17, col= "cornflowerblue")
      else (points(x=pos[p], y = currentDOM[p], t ='p', pch = 17, col= "black"))
    }
      if ((currentADD[p] >  currentADDDOM[p]) && (currentADD[p] > currentDOM[p])){
        if (chr %in% seq(1,20,2))
          points(x=pos[p], y = currentADD[p], t ='p', pch = 18, col= "cornflowerblue")
	else (points(x=pos[p], y = currentADD[p], t ='p', pch = 18, col= "black"))
    }
  }
}
axis(1, chrs, at = chrmids)
abline(h= 3.8, col="green",lty=3)
abline(h=4.2, col="orange",lty=3)
axis(1, chrs, at = chrmids)
abline(h= 3.8, col="green",lty=3)
abline(h=4.2, col="orange",lty=3)
legend("topright", bg="gray",
  legend = c("ADD + DOM dev.", "ADD", "DOM dev."),
  pch = c(16, 18, 17))
