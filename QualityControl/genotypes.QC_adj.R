# Quality control on genotypes
#
# copyright (c) 2018-2021 - Brockmann group - HU Berlin, Manuel Delpero & Danny Arends
# last modified September, 2019
# first written August, 2019

setwd("C:/Users/Manuel/Desktop/AIL_S1xS2/RAWDATA")

locusxdnaheader <- unlist(strsplit(readLines("Humboldt_Univ_Zu_Berlin_MURGIGV01_20191205_LocusXDNA.csv", n=16)[16],","))
locusxdnasnps <- unlist(strsplit(readLines("Humboldt_Univ_Zu_Berlin_MURGIGV01_20191205_LocusXDNA.csv", n=18)[18],","))

locusxdna <- readLines("Humboldt_Univ_Zu_Berlin_MURGIGV01_20191205_LocusXDNA.csv")[-c(1:22)]
splitted <- strsplit(locusxdna, ",")

calls <- matrix(NA, length(locusxdna) / 2, length(splitted[[1]]))
scores <- matrix(NA, length(locusxdna) / 2, length(splitted[[1]]))
for(x in 1:length(splitted)){
  if(x %% 2 == 1) calls[x/2,] <- splitted[[x]]
  if(x %% 2 == 0) scores[x/2,] <- splitted[[x]]
}

markers <- locusxdnaheader[4:length(locusxdnaheader)]

colnames(calls) <- c("Label", "plateWell", "Date","oligoPoolId","bundleId", "status", "Type", "Nas", markers)
colnames(scores) <- c("Label", "plateWell", "Date","oligoPoolId","bundleId", "status", "Type", "Nas", markers)

gts <- calls[,markers]
rownames(gts) <- gsub("V 888-", "AIL", calls[, "Label"])
qual <- apply(scores[,markers],2,as.numeric)
rownames(qual) <- gsub("V 888-", "AIL", calls[, "Label"])

# Filter by quality score
gts[qual < 0.7] <- NA
gts[gts == "U"] <- NA

# missing markers
gts <- t(gts)
idx <- which(apply(gts,1, function(x){sum(is.na(x)) == length(x)}))
gts <- gts[-idx,]

# non segregating markers
idx <- which(apply(gts,1,function(x){length(table(x)) == 1}))
gts <- gts[-idx,]

genotypes <- gts 

# At least 2 groups with 10 observations
good <- which(unlist(lapply(apply(genotypes,1,table), function(x){
  length(which(x > 10)) >= 2
})))
genotypes <- genotypes[good, ]
dim(genotypes)

# Prevent small groups, set the groups with < 10 individuals to NA
mtab <- apply(genotypes,1,table)
for(x in 1:nrow(genotypes)){
  if(any(mtab[[x]] < 6)){
    gt <- names(mtab[[x]])[which(mtab[[x]] < 6)]
    genotypes[x, genotypes[x, ] == gt] <- NA
  }
}
dim(genotypes)

# No duplicated markers, keep the first one we see
genotypes <- genotypes[-which(duplicated(genotypes)),]
dim(genotypes)
write.table(gts, "genotypes.raw.txt", sep="\t", quote=FALSE)

# load in the map file from Karl Broman
map <- read.table("snp_map.karl.txt", sep = ",", header = TRUE, row.names=1)
map <- map[rownames(gts),]

chrs <- 1:21 
names(chrs) <- c(1:19, "X", "Y")

plot(c(1,21), c(0,200000000), t = 'n', xaxt = "n", las= 2, ylab = "Pos", xlab = "Chr")
aa <- apply(map, 1, function(r) { points(x = chrs[r[1]], y = r[2], pch = "-"); })
axis(1, at = chrs, names(chrs))
