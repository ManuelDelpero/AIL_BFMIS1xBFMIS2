# AIL_S1xS2 Analysis on Phenotypes
#
# copyright (c) - Manuel Delpero
# first written April, 2019
#
# Script for selection of animals to genotype with the gigaMUGA array

setwd("C:/Users/Manuel/Desktop/AIL_S1xS2/RAWDATA")
mriLEAN <- read.csv("MRIlean.txt", sep = "\t", header=TRUE, check.names=FALSE)
mriFAT <- read.csv("MRIfat.txt", sep = "\t", header=TRUE, check.names=FALSE)
allPhenotypes <- read.csv("allPhenotypes.txt",sep = "\t", header=TRUE, check.names=FALSE)

# extremes for MRI
lowfat <- c()
highfat <- c()
for (x in 1:nrow(mriFAT)){
  if (!is.na(mriFAT[x,"174"])){
    if (mriFAT[x,"174"] > 18){
      highfat <- c(highfat, rownames(mriFAT[x,]))
      }else if (mriFAT[x,"174"] < 12){
      lowfat <- c(lowfat, rownames(mriFAT[x,]))
    }
  }
}

# extremes for bodyweight
weight <- allPhenotypes[, "D174"]
weight <- data.frame(weight)
weight <- cbind(rownames(allPhenotypes), weight)
weightordered <- weight[order(as.numeric(weight[,2]) ,decreasing = TRUE),]
highweight <- weightordered[1:100,]
highweight <- as.character(highweight[,1])
lowweight <- weightordered[(nrow(weightordered) - 100):nrow(weightordered),]

# extremes for triglycerides
Triglyceridesordered <- Triglycerides[order(Triglycerides[,2] ,decreasing = TRUE),]
HighTrig <- Triglyceridesordered[1:150,]
HighTrig <- as.character(HighTrig[,1])

# Get rid of NAs and sort to select the Low extremes
Trigs <- c()
for (x in 1:nrow(Triglyceridesordered)){
  if (!is.na(Triglyceridesordered[x, 2])){
    trig <- Triglyceridesordered[x, c(1,2)]
	Trigs <- rbind(trig,Trigs)
  }
}

LowTrig <- Trigs[1:150,1]

# extremes for gluc172
Gluc <- allPhenotypes[,"Gluc172"]
Gluc <- data.frame(Gluc)
Gluc <- cbind(rownames(allPhenotypes), Gluc)
Glucordered <- Gluc[order(Gluc[,2] ,decreasing = TRUE),]

highgluc <- Glucordered[1:100,]
highgluc <- highgluc[,1]
highgluc <- as.character(highgluc)
lowwgluc <- Glucordered[(nrow(Glucordered) - 100):nrow(Glucordered),]

# extremes for Gon weight
gonWeight <- allPhenotypes[,"Gon"]
gonWeight <- data.frame(gonWeight)
gonWeight <- cbind(rownames(allPhenotypes), gonWeight)
gonWeightordered <- gonWeight[order(gonWeight[,2] ,decreasing = TRUE),]
highgGonWeight <- gonWeightordered[1:150,]
highgGonWeight <- highgGonWeight[,1]
highgGonWeight <- as.character(highgGonWeight)
lowGonWeight <- gonWeightordered[(nrow(gonWeightordered) - 150):nrow(gonWeightordered),]
lowGonWeight <- lowGonWeight[,1]
lowGonWeight <- as.character(lowGonWeight)
boxplot(main = "extremes Gonadal Fat weight", gonWeightordered[1:100,2], gonWeightordered[(nrow(gonWeightordered) - 150):nrow(gonWeightordered),2], notch = TRUE)

# extremes to select (we mainly want to focus on th Gon fat and Liver due to the evident phenotypes)
Group1 <- highgGonWeight[which(highgGonWeight %in% as.character(LowTrig))]
Group2 <- lowGonWeight[which(lowGonWeight %in% as.character(HighTrig))]
boxplot(main = "extremes Gonadal Fat weight and triglycerides", gonWeightordered[which(as.numeric(as.character(gonWeightordered[,1])) %in% as.numeric(Group1)), 2], gonWeightordered[which(as.numeric(as.character(gonWeightordered[,1])) %in% as.numeric(Group2)), 2], notch = TRUE)

Group1 <- gonWeightordered[which(as.numeric(as.character(gonWeightordered[,1])) %in% as.numeric(Group1)),]
Group2 <- gonWeightordered[which(as.numeric(as.character(gonWeightordered[,1])) %in% as.numeric(Group2)),]
colnames(Group1) <- c("ID", "GonFatWeight")
colnames(Group2) <- c("ID", "GonFatWeight")
Group1 <- Group1[order(as.numeric(Group1[,1]) ,decreasing = FALSE),]
Group2 <- Group2[order(as.numeric(Group2[,1]) ,decreasing = FALSE),]
write.table(Group1, file = "Group1.txt", sep ="\t", row.names = FALSE)
write.table(Group2, file = "Group2.txt", sep ="\t", row.names = FALSE)

if (i == 0){
highextremes4 <- c()
highextreme4 <- c()
highextremes3 <- c() 
highextreme3 <- c()
highextremes2 <- c()
highextreme2 <- c()
lowextremes4 <- c()
lowextreme4 <- c()
lowextremes3 <- c() 
lowextreme3 <- c()
lowextremes2 <- c()
lowextreme2 <- c()

for (x in 1:nrow(allPhenotypes)){
  if ((rownames(allPhenotypes[x,]) %in% highfat) && (rownames(allPhenotypes[x,]) %in% highweight) && (rownames(allPhenotypes[x,]) %in% HighTrig) && (rownames(allPhenotypes[x,]) %in% highgluc)) 
    highextreme4 <- rownames(allPhenotypes[x,])
    highextremes4 <- unique(c(highextreme4, highextremes4))
  if ((rownames(allPhenotypes[x,]) %in% highfat) && (rownames(allPhenotypes[x,]) %in% highweight) && (rownames(allPhenotypes[x,]) %in% HighTrig) || (rownames(allPhenotypes[x,]) %in% highfat) && (rownames(allPhenotypes[x,]) %in% highweight) && (rownames(allPhenotypes[x,]) %in% highgluc) || (rownames(allPhenotypes[x,]) %in% highfat) && (rownames(allPhenotypes[x,]) %in% highgluc) && (rownames(allPhenotypes[x,]) %in% HighTrig) || (rownames(allPhenotypes[x,]) %in% highgluc) && (rownames(allPhenotypes[x,]) %in% highweight) && (rownames(allPhenotypes[x,]) %in% HighTrig)) 
    highextreme3 <- rownames(allPhenotypes[x,])
    highextremes3 <- unique(c(highextreme3, highextremes3))
  if ((rownames(allPhenotypes[x,]) %in% highfat) && (rownames(allPhenotypes[x,]) %in% highweight) || (rownames(allPhenotypes[x,]) %in% highfat) && (rownames(allPhenotypes[x,]) %in% HighTrig) || (rownames(allPhenotypes[x,]) %in% highfat) && (rownames(allPhenotypes[x,]) %in% highgluc) || (rownames(allPhenotypes[x,]) %in% highweight) && (rownames(allPhenotypes[x,]) %in% HighTrig) || (rownames(allPhenotypes[x,]) %in% highweight) && (rownames(allPhenotypes[x,]) %in% highgluc) || (rownames(allPhenotypes[x,]) %in% HighTrig) && (rownames(allPhenotypes[x,]) %in% highgluc))
    highextreme2 <- rownames(allPhenotypes[x,])
    highextremes2 <- unique(c(highextreme2, highextremes2))
}
}