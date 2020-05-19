# Various plots to represent the results for AIL BFMI S1xS2
#
# copyright (c) 2018-2021 - Brockmann group - HU Berlin Manuel Delpero 
# 
# first written december, 2019

setwd("C:/Users/Manuel/Desktop/AIL_S1xS2/RAWDATA")

phenotypes <- read.csv("allPhenotypes.txt", header = TRUE, check.names = FALSE, sep = "\t", colClasses = "character")
genotypes <- read.csv("genotypesComplete.txt", header = TRUE, check.names = FALSE, sep="\t", colClasses="character")
lodmatrixDOM <- read.table("lodmatrixDOM_nosum.txt", header = TRUE, sep = "\t", check.names = FALSE)
lodmatrixADD <- read.table("lodmatrixADD.txt", header = TRUE, sep = "\t", check.names = FALSE)
mprofiles <- read.table("lodmatrixADDDOM_nosum.txt", header = TRUE, sep = "\t", check.names = FALSE)
markerannot <- read.csv("map.cleaned.txt", header=TRUE, sep="\t", check.names=FALSE)
markerannot <- markerannot[, c(1,2)]
markerannot$Index <- seq.int(nrow(markerannot))
colnames(markerannot) <- c("Chromosome", "Position", "Index")
markerannot <- markerannot[-which(is.na(markerannot[,2])),]
colnames(genotypes) <- gsub("AIL", "", colnames(genotypes))
phenotypes <- phenotypes[colnames(genotypes),]
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
axis(1, at = c(0,25000000, 50000000, 75000000, 100000000), c("0", "25", "50", "75", "100"))


# Liver chr 17
plot(main = "QTL liver weight [Chr 17]", c(min(as.numeric(chr17ord[, "Position"])), max(as.numeric(chr17ord[, "Position"]))), c(0,9), ylab = "-log10 [pvalue]", xlab = "Position [mb]", las = 2, t = "n", xaxt = "n")
points(x = as.numeric(chr17ord[,"Position"]), y = chr17ord[,"Leber"] , type = "l", col="dodgerblue", lwd = 1)
axis(1, at = c(0,25000000, 50000000, 75000000, 100000000), c("0", "25", "50", "75", "100"))

# Body weight different time points for chr 15 (main QTL)
lodannotmatrix <- cbind(annotation[rownames(lodmatrixADD), ], lodmatrixADD)
dataset <- lodannotmatrix[, c("Chromosome", "Position", "D28", "D49", "D35", "D42", "D49", "D56", "D63", "D70", "D77", "D84", "D91", "D98", "D105", "D112", "D119", "D126", "D133", "D140", "D160", "D172")]
chr15 <- dataset[which(dataset[,"Chromosome"] == 15),]
chr15 <- chr15[order(chr15[,"Position"]),]
plot(main = "QTL profile bodyweight [Chr 15]", c(min(as.numeric(chr15[, "Position"])), max(as.numeric(chr15[, "Position"]))), c(0,9), ylab = "-log10 [pvalue]", xlab = "Position [mb]", las = 2, t = "n", xaxt = "n")
  points(x = as.numeric(chr15[,"Position"]), y = chr15[,"D70"] , type = "l", col="dodgerblue", lwd = 1)
  points(x = as.numeric(chr15[,"Position"]), y = chr15[,"D140"] , type = "l", col="blue", lwd = 1)
  points(x = as.numeric(chr15[,"Position"]), y = chr15[,"D63"] , type = "l", col="deepskyblue", lwd = 1)
  points(x = as.numeric(chr15[,"Position"]), y = chr15[,"D126"] , type = "l", col="purple", lwd = 1)
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


## Manhattan plots (Plot the effect with the highest  lod score and use three different symbols for each one) 
par(cex.lab=1.4, cex.main = 1.8, cex.axis = 1.5)
mat <- matrix(c(1,1,2,3), 2, 2, byrow = TRUE)
layout(mat, widths = rep.int(3, ncol(mat)))

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
plot(x = c(-gap, tail(chr.starts,1)), y = c(0,10), t = 'n', xlab="Chromosome", ylab="-log10[P]",xaxt='n', xaxs="i", yaxs="i", las=2, main=paste0("Manhattan plot - Gonadal adipose tissue weight"))
for(chr in chrs){
  onChr <- rownames(map.sorted[map.sorted[,"Chromosome"] == chr,])
  currentADDDOM <- lodmatrixADDDOM[onChr, phenotype]
  currentDOM <- lodmatrixDOM[onChr, phenotype]
  currentADD <- lodmatrixADD[onChr, phenotype]
  if (chr == "X")
    points(x=chr.starts[chr] + map.sorted[onChr,"Position"], y = mprofiles[onChr, phenotype], t ='p', pch = 16, cex = 1.5, col= "black")
  for (p in 1:length(currentADDDOM)){
    pos <- chr.starts[chr] + map.sorted[onChr,"Position"]
    if ((currentADDDOM[p] >  currentDOM[p]) && (currentADDDOM[p] > currentADD[p])){
      if (chr %in% seq(1,20,2))
        points(x=pos[p], y = currentADDDOM[p], t ='p', pch = 0, cex = 1.5, col= "cornflowerblue")
      else (points(x=pos[p], y = currentADDDOM[p], t ='p', pch = 0, cex = 1.5, col= "black"))
    }
    if ((currentDOM[p] >  currentADDDOM[p]) && (currentDOM[p] > currentADD[p])){
      if (chr %in% seq(1,20,2))
        points(x=pos[p], y = currentDOM[p], t ='p', pch = 17, cex = 1,5, col= "cornflowerblue")
      else (points(x=pos[p], y = currentDOM[p], t ='p', pch = 17, cex = 1.5, col= "black"))
    }
      if ((currentADD[p] >  currentADDDOM[p]) && (currentADD[p] > currentDOM[p])){
        if (chr %in% seq(1,20,2))
          points(x=pos[p], y = currentADD[p], t ='p', pch = 18, cex = 1.5, col= "cornflowerblue")
	else (points(x=pos[p], y = currentADD[p], t ='p', pch = 18, cex = 1.5,col= "black"))
    }
  }
}
axis(1, chrs, at = chrmids)
abline(h= 4, col="orange",lty=3)
abline(h= 4.5, col="green",lty=3)
axis(1, chrs, at = chrmids)
legend("topright", #bg="gray"
  bty = "n",
  legend = c("Dominance", "Additive", "Dominance dev."),
  pch = c(0, 18, 17))
  

 # Effect plot for the top marker on chr 3
UNCHS043909 <- cbind(phenotypes[, "Gon"], t(genotypes["UNCHS043909",]))
UNC5791802 <- cbind(phenotypes[, "Gon"], t(genotypes["UNC5791802",]))
boxplot(as.numeric(UNC5791802[which(UNC5791802[,2] == "B"),1]), as.numeric(UNC5791802[which(UNC5791802[,2] == "H"),1]), as.numeric(UNC5791802[which(UNC5791802[,2] == "A"),1]), col = c("lightskyblue1", "cyan3", "dodgerblue4"), main = "Gonadal adipose tissue weight [Chr 3]", ylab = "Weight [Gr.]", xlab = "Genotypes [UNC5791802]", las = 2, t = "n", xaxt = "n", ylim = c(0, 7))
  axis(1, at = 1:3 , c("TT", "TC", "CC"))
  legend("topright", #bg="gray",
  legend = c( "BFMI-S1", "HET", "BFMI-S2"),
    col = c("lightskyblue1", "cyan3", "dodgerblue4"),
    pch = 15,
    pt.cex = 1.7,
    #pt.bg = "lightsteelblue1",
    cex = 1,
	bty = "n",
    text.col = "black")
   
# Effect plot for the top marker on chr 17
boxplot(as.numeric(UNCHS043909[which(UNCHS043909[,2] == "A"),1]), as.numeric(UNCHS043909[which(UNCHS043909[,2] == "H"),1]), as.numeric(UNCHS043909[which(UNCHS043909[,2] == "B"),1]), col = c("lightskyblue1", "cyan3", "dodgerblue4"), main = "Gonadal adipose tissue weight [Chr 17]", ylab = "Weight [Gr.]", xlab = "Genotypes [UNCHS043909]" , las = 2, t = "n", xaxt = "n",  ylim = c(0, 7))
  axis(1, at = 1:3 , c("TT", "TC", "CC"))
  legend("topright", #bg="gray",
  legend = c( "BFMI-S1", "HET", "BFMI-S2"),
    col = c("lightskyblue1", "cyan3", "dodgerblue4"),
    pch = 15,
    pt.cex = 1.7,
    #pt.bg = "lightsteelblue1",
    cex = 1,
	bty = "n",
    text.col = "black")
  
# Body weight week 18
phenotype <- "D126"
plot(x = c(-gap, tail(chr.starts,1)), y = c(0,15), t = 'n', xlab="Chromosome", ylab="-log10[P]",xaxt='n', xaxs="i", yaxs="i",las=2,main=paste0("Manhattan plot - Bodyweight week 18"))
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
        points(x=pos[p], y = currentADDDOM[p], t ='p', pch = 0, col= "cornflowerblue")
      else (points(x=pos[p], y = currentADDDOM[p], t ='p', pch = 0, col= "black"))
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
abline(h= 4, col="orange",lty=3)
abline(h= 4.5, col="green",lty=3)
axis(1, chrs, at = chrmids)
legend("topright", #bg="gray"
  bty = "n",
  legend = c("Dominance", "Additive", "Dominance dev."),
  pch = c(0, 18, 17))
  



# Effect plot for the top marker on chr 15
UNCHS041907 <- cbind(phenotypes[, "D126"], t(genotypes["UNC25922623",]))
boxplot(as.numeric(UNCHS041907[which(UNCHS041907[,2] == "B"),1]), as.numeric(UNCHS041907[which(UNCHS041907[,2] == "H"),1]), as.numeric(UNCHS041907[which(UNCHS041907[,2] == "A"),1]), col = c("lightskyblue1", "cyan3", "dodgerblue4"), main = "Body weight week 18 [Chr 15]", ylab = "Weight [Gr.]", xlab = "Genotypes [UNC25922623]" , las = 2, t = "n", xaxt = "n",  ylim = c(20, 60))
  axis(1, at = 1:3 , c("TT", "TC", "CC"))
  legend("topright", bg="gray",
  legend = c( "BFMI-S1", "HET", "BFMI-S2"),
    col = c("lightskyblue1", "cyan3", "dodgerblue4"),
    pch = 15,
    pt.cex = 1.7,
    pt.bg = "lightsteelblue1",
    cex = 1,
	bty = "n",
    text.col = "black")

# Effect plot for the top marker on chr 16 
UNCHS041907 <- cbind(phenotypes[, "D126"], t(genotypes["UNCHS041907",]))
boxplot(as.numeric(UNCHS041907[which(UNCHS041907[,2] == "A"),1]), as.numeric(UNCHS041907[which(UNCHS041907[,2] == "H"),1]), as.numeric(UNCHS041907[which(UNCHS041907[,2] == "B"),1]), col = c("lightskyblue1", "cyan3", "dodgerblue4"), main = "Body weight week 18 [Chr 16]", ylab = "Weight [Gr.]", xlab = "Genotypes [UNCHS041907]" , las = 2, t = "n", xaxt = "n",  ylim = c(20, 60))
  axis(1, at = 1:3 , c("TT", "TC", "CC"))
  legend("topright", bg="gray",
  legend = c( "BFMI-S1", "HET", "BFMI-S2"),
    col = c("lightskyblue1", "cyan3", "dodgerblue4"),
    pch = 15,
    pt.cex = 1.7,
    pt.bg = "lightsteelblue1",
    cex = 1,
	bty = "n",
    text.col = "black")
	

  
# liver weight
phenotype <- "Leber"
plot(x = c(-gap, tail(chr.starts,1)), y = c(0,8), t = 'n', xlab="Chromosome", ylab="-log10[P]",xaxt='n', xaxs="i", yaxs="i",las=2,main=paste0("Manhattan plot - liver weight"))
for(chr in chrs){
  onChr <- rownames(map.sorted[map.sorted[,"Chromosome"] == chr,])
  currentADDDOM <- mprofiles[onChr, phenotype]
  currentDOM <- lodmatrixDOM[onChr, phenotype]
  currentADD <- lodmatrixADD[onChr, phenotype]
  if (chr == "X")
    points(x=chr.starts[chr] + map.sorted[onChr,"Position"], y = mprofiles[onChr, phenotype], t ='p', pch = 16, cex =1.5, col= "black")
  for (p in 1:length(currentADDDOM)){
    pos <- chr.starts[chr] + map.sorted[onChr,"Position"]
    if ((currentADDDOM[p] >  currentDOM[p]) && (currentADDDOM[p] > currentADD[p])){
      if (chr %in% seq(1,20,2))
        points(x=pos[p], y = currentADDDOM[p], t ='p', pch = 0, cex = 1.5, col= "cornflowerblue")
      else (points(x=pos[p], y = currentADDDOM[p], t ='p', pch = 0, cex = 1.5, col= "black"))
    }
    if ((currentDOM[p] >  currentADDDOM[p]) && (currentDOM[p] > currentADD[p])){
      if (chr %in% seq(1,20,2))
        points(x=pos[p], y = currentDOM[p], t ='p', pch = 17, cex = 1.5, col= "cornflowerblue")
      else (points(x=pos[p], y = currentDOM[p], t ='p', pch = 17, cex = 1.5, col= "black"))
    }
      if ((currentADD[p] >  currentADDDOM[p]) && (currentADD[p] > currentDOM[p])){
        if (chr %in% seq(1,20,2))
          points(x=pos[p], y = currentADD[p], t ='p', pch = 18, cex = 1.5, col= "cornflowerblue")
	else (points(x=pos[p], y = currentADD[p], t ='p', pch = 18, cex = 1.5, col= "black"))
    }
  }
}
axis(1, chrs, at = chrmids)
abline(h= 4, col="orange",lty=3)
abline(h= 4.5, col="green",lty=3)
axis(1, chrs, at = chrmids)
legend("topright", #bg="gray"
  bty = "n",
  legend = c("Dominance", "Additive", "Dominance dev."),
  pch = c(0, 18, 17))

# Effect plot for the top marker on chr 17 for the liver
UNCHS043909 <- cbind(phenotypes[, "Leber"], t(genotypes["UNCHS043909",]))
UNC27568354 <- cbind(phenotypes[, "Leber"], t(genotypes["UNC27568354",]))
boxplot(as.numeric(UNCHS043909[which(UNCHS043909[,2] == "A"),1]), as.numeric(UNCHS043909[which(UNCHS043909[,2] == "H"),1]), as.numeric(UNCHS043909[which(UNCHS043909[,2] == "B"),1]), col = c("lightskyblue1", "cyan3", "dodgerblue4"), main = "Liver weight [Chr 17]", ylab = "Weight [Gr.]", xlab = "Genotypes [UNCHS043909]" , las = 2, t = "n", xaxt = "n",  ylim = c(0, 8))
  axis(1, at = 1:3 , c("TT", "TC", "CC"))
  legend("topright", #bg="gray",
  legend = c( "BFMI-S1", "HET", "BFMI-S2"),
    col = c("lightskyblue1", "cyan3", "dodgerblue4"),
    pch = 15,
    pt.cex = 1.7,
    pt.bg = "lightsteelblue1",
    cex = 1,
	bty = "n",
    text.col = "black")

boxplot(as.numeric(UNC27568354[which(UNC27568354[,2] == "A"),1]), as.numeric(UNC27568354[which(UNC27568354[,2] == "H"),1]), as.numeric(UNC27568354[which(UNC27568354[,2] == "B"),1]), col = c("lightskyblue1", "cyan3", "dodgerblue4"), main = "Liver weight [Chr 17]", ylab = "Weight [Gr.]", xlab = "Genotypes [UNC27568354]" , las = 2, t = "n", xaxt = "n",  ylim = c(0, 8))
  axis(1, at = 1:3 , c("TT", "TC", "CC"))
  legend("topright", #bg="gray",
  legend = c( "BFMI-S1", "HET", "BFMI-S2"),
    col = c("lightskyblue1", "cyan3", "dodgerblue4"),
    pch = 15,
    pt.cex = 1.7,
    pt.bg = "lightsteelblue1",
    cex = 1,
	bty = "n",
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
abline(h= 4, col="orange",lty=3)
abline(h= 4.5, col="green",lty=3)
legend("topright", bg="gray",
  legend = c("ADD + DOM dev.", "ADD", "DOM dev."),
  pch = c(16, 18, 17))

# D140
phenotype <- "D140"
plot(x = c(-gap, tail(chr.starts,1)), y = c(0,10), t = 'n', xlab="Chromosome", ylab="-log10[P]" xaxt='n', xaxs="i", yaxs="i",las=2,main=paste0("Manhattan plot - ", phenotype))
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
abline(h= 4, col="orange",lty=3)
abline(h= 4.5, col="green",lty=3)
axis(1, chrs, at = chrmids)
legend("topright", bg="gray",
  legend = c("ADD + DOM dev.", "ADD", "DOM dev."),
  pch = c(16, 18, 17))

# Gluc172
phenotype <- "Gluc172"
plot(x = c(-gap, tail(chr.starts,1)), y = c(0,8), t = 'n', xlab="Chromosome", ylab="-log10[P]",xaxt='n', xaxs="i", yaxs="i",las=2,main=paste0("Manhattan plot - Glucose day 172"))
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
        points(x=pos[p], y = currentADDDOM[p], t ='p', pch = 0, cex = 1.5, col= "cornflowerblue")
      else (points(x=pos[p], y = currentADDDOM[p], t ='p', pch = 0, cex = 1.5, col= "black"))
    }
    if ((currentDOM[p] >  currentADDDOM[p]) && (currentDOM[p] > currentADD[p])){
      if (chr %in% seq(1,20,2))
        points(x=pos[p], y = currentDOM[p], t ='p', pch = 17, cex = 1.5, col= "cornflowerblue")
      else (points(x=pos[p], y = currentDOM[p], t ='p', pch = 17, cex = 1.5, col= "black"))
    }
      if ((currentADD[p] >  currentADDDOM[p]) && (currentADD[p] > currentDOM[p])){
        if (chr %in% seq(1,20,2))
          points(x=pos[p], y = currentADD[p], t ='p', pch = 18, cex = 1.5, col= "cornflowerblue")
	else (points(x=pos[p], y = currentADD[p], t ='p', pch = 18, cex = 1.5, col= "black"))
    }
  }
}
axis(1, chrs, at = chrmids)
abline(h= 4, col="orange",lty=3)
abline(h= 4.5, col="green",lty=3)
axis(1, chrs, at = chrmids)
legend("topright", #bg="gray"
  bty = "n",
  legend = c("Dominance", "Additive", "Dominance dev."),
  pch = c(0, 18, 17))

# Effect plot for top marker for Gon weight (chr 17) using glucose as phenotypes, this marker is responsile for ectopic fat storage and also for the glucose level 
JAX00432128 <- cbind(phenotypes[, "Gluc172"], t(genotypes["JAX00432128",]))
UNC5791802 <- cbind(phenotypes[, "Gluc172"], t(genotypes["UNC5791802",]))
boxplot(as.numeric(JAX00432128[which(JAX00432128[,2] == "A"),1]), as.numeric(JAX00432128[which(JAX00432128[,2] == "H"),1]), as.numeric(JAX00432128[which(JAX00432128[,2] == "B"),1]), col = c("lightskyblue1", "cyan3", "dodgerblue4"), main = "Effect plot marker JAX00432128 [Chr 17]", ylab = "Glucose [mg/dL]", xlab = "Genotypes[JAX00432128]" , ylim = c(0, 600), las = 2, t = "n", xaxt = "n")
  axis(1, at = 1:3 , c("TT", "TC", "CC"))
  legend("topright", #bg="gray",
  legend = c( "BFMI-S1", "HET", "BFMI-S2"),
    col = c("lightskyblue1", "cyan3", "dodgerblue4"),
    pch = 15,
    pt.cex = 1.7,
	bty = "n",
    pt.bg = "lightsteelblue1",
    cex = 1,
    text.col = "black")
	
# Effect plot for top marker for Gon weight (chr3) using glucose as phenotypes, this marker is responsile for ectopic fat storage and also for the glucose level 
boxplot(as.numeric(UNC5791802[which(UNC5791802[,2] == "B"),1]), as.numeric(UNC5791802[which(UNC5791802[,2] == "H"),1]), as.numeric(UNC5791802[which(UNC5791802[,2] == "A"),1]), col = c("lightskyblue1", "cyan3", "dodgerblue4"), main = "Effect plot marker UNC5791802 [Chr 3]", ylab = "Glucose [mg/dL]", xlab = "Genotypes[UNC5791802]",  ylim = c(0, 600), las = 2, t = "n", xaxt = "n")
  axis(1, at = 1:3 , c("TT", "TC", "CC"))
  legend("topright", #bg="gray",
  legend = c( "BFMI-S1", "HET", "BFMI-S2"),
    col = c("lightskyblue1", "cyan3", "dodgerblue4"),
    pch = 15,
    pt.cex = 1.7,
	bty = "n",
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
  
  
# Effect plot for top marker for Gon weight (chr 17) using ITT auc as phenotypes, this marker is responsile for ectopic fat storage and also for the glucose level 
JAX00432128 <- cbind(phenotypes[, "ITTauc"], t(genotypes["JAX00432128",]))
UNC5791802 <- cbind(phenotypes[, "ITTauc"], t(genotypes["UNC5791802",]))
boxplot(as.numeric(JAX00432128[which(JAX00432128[,2] == "A"),1]), as.numeric(JAX00432128[which(JAX00432128[,2] == "H"),1]), as.numeric(JAX00432128[which(JAX00432128[,2] == "B"),1]), col = c("lightskyblue1", "cyan3", "dodgerblue4"), main = "Effect plot marker JAX00432128 [Chr 17]", ylab = "Weight [Gr.]", xlab = "Genotypes" , las = 2, t = "n", xaxt = "n")
  axis(1, at = 1:3 , c("TT", "TC", "CC"))
  legend("topright", bg="gray",
  legend = c( "BFMI-S1", "HET", "BFMI-S2"),
    col = c("lightskyblue1", "cyan3", "dodgerblue4"),
    pch = 15,
    pt.cex = 1.7,
    pt.bg = "lightsteelblue1",
    cex = 1,
    text.col = "black")
	
# Effect plot for top marker for Gon weight (chr3) using ITT auc as phenotypes, this marker is responsile for ectopic fat storage and also for the glucose level 
boxplot(as.numeric(UNC5791802[which(UNC5791802[,2] == "A"),1]), as.numeric(UNC5791802[which(UNC5791802[,2] == "H"),1]), as.numeric(UNC5791802[which(UNC5791802[,2] == "B"),1]), col = c("lightskyblue1", "cyan3", "dodgerblue4"), main = "Effect plot marker UNC5791802 [Chr 3]", ylab = "Weight [Gr.]", xlab = "Genotypes" , las = 2, t = "n", xaxt = "n")
  axis(1, at = 1:3 , c("TT", "TC", "CC"))
  legend("topright", bg="gray",
  legend = c( "BFMI-S1", "HET", "BFMI-S2"),
    col = c("lightskyblue1", "cyan3", "dodgerblue4"),
    pch = 15,
    pt.cex = 1.7,
    pt.bg = "lightsteelblue1",
    cex = 1,
    text.col = "black")


# length
phenotype <- "LÃ¤nge"
plot(x = c(-gap, tail(chr.starts,1)), y = c(0,12), t = 'n', xlab="Chromosome", ylab="-log10(P)",xaxt='n', xaxs="i", yaxs="i",las=2,main="Manhattan plot - Length")
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
        points(x=pos[p], y = currentADDDOM[p], t ='p', pch = 0, col= "cornflowerblue")
      else (points(x=pos[p], y = currentADDDOM[p], t ='p', pch = 0, col= "black"))
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
abline(h= 4, col="orange",lty=3)
abline(h= 4.5, col="green",lty=3)
axis(1, chrs, at = chrmids)
legend("topright", #bg="gray"
  bty = "n",
  legend = c("Dominance", "Additive", "Dominance dev."),
  pch = c(0, 18, 17))
  
# Effect plot for top marker for the length (chr16) 
UNCHS041714 <- cbind(phenotypes[, "LÃ¤nge"], t(genotypes["UNCHS041714",]))
boxplot(as.numeric(UNCHS041714[which(UNCHS041714[,2] == "A"),1]), as.numeric(UNCHS041714[which(UNCHS041714[,2] == "H"),1]), as.numeric(UNCHS041714[which(UNCHS041714[,2] == "B"),1]), col = c("lightskyblue1", "cyan3", "dodgerblue4"), main = "Effect plot marker UNCHS041714 [Chr 16]", ylab = "Length [cm]", ylim = c(9.5,13), xlab = "Genotypes" , las = 2, t = "n", xaxt = "n")
  axis(1, at = 1:3 , c("TT", "TC", "CC"))
  legend("topright", #bg="gray",
  legend = c( "BFMI-S1", "HET", "BFMI-S2"),
    col = c("lightskyblue1", "cyan3", "dodgerblue4"),
    pch = 15,
    pt.cex = 1.7,
    pt.bg = "lightsteelblue1",
    cex = 1,
    text.col = "black")
  
