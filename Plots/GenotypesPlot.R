setwd("C:/Users/Manuel/Desktop/AIL_S1xS2/RAWDATA")

mdata <- read.table("genotypes.cleaned.txt")
map <- read.table("map.cleaned.txt", colClasses="character")

ord <- order(as.numeric(as.character(map[, "bp_mm10"])))

map <- map[ord,]


nmap <- NULL
for(chr in c(1:20, "X")){
  nmap <- rbind(nmap, map[which(map[, "chr"] == chr),])
}
map <- nmap

mdata <- mdata[rownames(map),]

mndata <- apply(mdata,2,function(x){
  as.numeric(factor(x, levels = c("A", "H", "B")))
})

#mcor <- cor(t(mndata), use="pair")
#image(abs(mcor), xaxt= "n", yaxt= "n")

chrs <- c(1:19, "X", "Y")
getChrIdx <- function(map){
  chrstart <- 0
  nmap <- c()
  for(chr in chrs){
    nmar <- chrstart + length(which(map[, "chr"] == chr))
    nmap <- c(nmap, nmar)
    chrstart <- nmar
  }
  return(nmap)
}

linepos <- c(getChrIdx(map) / max(getChrIdx(map)))
abline(v = linepos)
linepos <- c(0, linepos)
ll <- c()
for (x in 1 : length(linepos)) {
  ll <- c(ll, linepos[x] + ((linepos[x + 1] - linepos[x] ) / 2))
}
ll <- na.omit(ll)

library(RColorBrewer)
colz <- brewer.pal(n = 9, name = "YlOrRd")

# Genetic map
pal <- colorRampPalette(colz)

corM <- cor(t(mndata), use = "pair")
image(corM, xaxt = 'n', yaxt = 'n')
abline(v = linepos, col = "white")
abline(h = linepos, col = "white")
axis(1, at = ll, chrs, las=1)
axis(2, at = ll, chrs, las=1)
box()


#OldRange = (nrow(map) - 1)  
#NewRange = (1 - 0)
#values <- c() 
#for (x in lengths){
#  NewValue = (((x - 1) * NewRange) / OldRange) + 0
#  values <- c(values, NewValue)
#}

# calculate het % for each marker and plot
mdata_table <- apply(mdata,1,table)
het_percent <- c()
for (x in 1:length(mdata_table)){
  if ("H" %in% names(mdata_table[[x]])){
    het <- mdata_table[[x]][["H"]]
  }else{het <- 0}
  het_percent <- c(het_percent, het)
}
het_percent <- (het_percent/200) * 100
het_percent <- cbind(het_percent, map[,"chr"])
rownames(het_percent) <- rownames(map)

chrs <- c(1:19,"X")
gap <- 80000000
map.sorted <- NULL
chr.lengths <- c()
chr.starts <- c(0)
chrmids <- c()
i <- 1
for(chr in chrs){
  onChr <- which(map[,"chr"] == chr)
  map.sorted <- rbind(map.sorted, map[onChr,])
  chr.lengths <- c(chr.lengths, max(as.numeric(map[onChr, "bp_mm10"])))
  chr.starts <- c(chr.starts, chr.starts[i] + max(as.numeric(map[onChr, "bp_mm10"])) + gap)
  i <- i + 1
}

chr.start <- chr.starts[-5]
chr.ends <- chr.start + chr.lengths
names(chr.starts) <- chrs
names(chr.lengths) <- chrs

for (x in chrs){
  chrmid <- as.numeric(chr.lengths[x]/2) + as.numeric(chr.starts[x])
  chrmids <- c(chrmids, chrmid)
}


plot(x = c(-gap, tail(chr.starts,1)), y = c(0,100), t = 'n', xlab="Chromosome", ylab="% of Het",xaxt='n', xaxs="i", yaxs="i", las=2)
for(chr in chrs){
  onChr <- rownames(map[map[,"chr"] == chr,])
  current <- as.numeric(het_percent[onChr, 1])
  if (chr == "X"){
    points(x=chr.starts[chr] + as.numeric(map[onChr,"bp_mm10"]), y = as.numeric(het_percent[onChr, 1]), t ='p', pch = 16, cex = 1, col= "orange")
  }else{
    pos <- chr.starts[chr] + as.numeric(map[onChr,"bp_mm10"])
    for (p in 1:length(current)){
      if (as.numeric(chr) %in% seq(1,20,2)){
        points(x=pos[p], y = current[p], type ='o', pch = 16, cex = 1, col= "gray61")
      }else{ points(x=pos[p], y = current[p], type ='o', pch = 16, cex = 1, col= "orange")}
    }
  }
}
axis(1, chrs, at = chrmids)