# AIL_S1xS2 Analysis on Phenotypes
#
# copyright (c) - Manuel Delpero
# first written April, 2019
#
# Script for analysis on Phenotypes

setwd("C:/Users/Manuel/Desktop/AIL_S1xS2/RAWDATA")

# Reading the files that contain the phenotypes
allPhenotypes <- read.csv("growthWeight",sep = "\t", header=TRUE, check.names=FALSE)
PlasmaGluc <- read.csv("Gluc147.txt",sep = "\t", header= FALSE, check.names=FALSE)
Triglycerides <- read.csv("Triglycerides.txt",sep = "\t", header= TRUE, check.names=FALSE)
Sex <- read.csv("sex.txt",sep = "\t", header= TRUE, check.names=FALSE)
Sex <- Sex[1:530,]
mriLEAN <- read.csv("MRIlean.txt", sep = "\t", header=TRUE, check.names=FALSE)
mriFAT <- read.csv("MRIfat.txt", sep = "\t", header=TRUE, check.names=FALSE)
oralGTTDATA <- read.csv("oralGTT.txt", header=TRUE, check.names=FALSE, sep="\t")
insulinDATA <- read.csv("insulinTTDATA.txt", header=TRUE, check.names=FALSE, sep="\t")
allTissues <- read.csv("All_Tissues.txt", header=TRUE, check.names=FALSE, sep="\t")
allTissuesIDs <- allTissues[,1] 
rownames(allTissues) <- allTissues[,1]		
allTissues <- allTissues[,-1]															# Update the allTissue file and keep the data for the animals that survived until the end (IDs contained in the allTissue file)

# Combine all the data into one data frame
allPhenotypes <- allPhenotypes[which(allPhenotypes[, "ID-Nr."] %in% Sex[,"ID"]),]
Sex <- Sex[which(Sex[,"ID"] %in% allPhenotypes[, "ID-Nr."]),]
allPhenotypes <- cbind(Sex[,2], allPhenotypes)
PlasmaGluc <- PlasmaGluc[-(1:7),]
PlasmaGluc <- PlasmaGluc[which(PlasmaGluc[,1] %in% allPhenotypes[, "ID-Nr."]),]
rownames(allPhenotypes) <- PlasmaGluc[,1]
PlasmaGluc <- PlasmaGluc[,2]
PlasmaGluc <- data.frame(PlasmaGluc)
colnames(PlasmaGluc ) <- "Gluc172"
allPhenotypes <- cbind(allPhenotypes, PlasmaGluc)
allPhenotype <- gsub("V 888-", "", rownames(allPhenotypes))
rownames(allPhenotypes) <- allPhenotype
out <- which(!(insulinDATA[,1] %in% oralGTTDATA[,1]))
oralGTTDATA <- oralGTTDATA[-out,]
rownames(oralGTTDATA) <- oralGTTDATA[,1]
oralGTTDATA <- oralGTTDATA[,-1]
testData <- cbind(oralGTTDATA, insulinDATA[,2:5])
testDatarows <- gsub("V 888-", "", rownames(testData))
rownames(testData) <- testDatarows
allPhenotypes <- allPhenotypes[which(rownames(allPhenotypes) %in% rownames(testData)),]
testData <- testData[which(rownames(testData) %in% rownames(allPhenotypes)),]
allPhenotypes <- cbind(allPhenotypes, testData)
allTissues <- allTissues[which(rownames(allTissues) %in% rownames(allPhenotypes)), ]
allPhenotypes <- allPhenotypes[which(rownames(allPhenotypes) %in% rownames(allTissues)), ]
allPhenotypes <- cbind(allPhenotypes, allTissues)
Triglycerides <- Triglycerides[which(Triglycerides[,1] %in% rownames(allPhenotypes)),]
allPhenotypes <- allPhenotypes[which(rownames(allPhenotypes) %in% Triglycerides[,1]),]
allPhenotypes <- cbind(allPhenotypes, Triglycerides[,2])
colallPhenotypes <- c("Sex", "ID" , "D21", "D28", "D35", "D42", "D49", "D56", "D63", "D70", "D77", "D84", "D91", "D98", "D105", "D112", "D119", "D125", "D126", "D133", "D139", "D140", "D142", "D144", "D147", "D150", "D154", "D157", "D160", "D163", "D166", "D169", "D172", "D174", "Gluc172"  ,   "0 min"    ,   "15 min",  "30 min"     , "60 min", "120 min"   ,  "0 min"   ,    "15 min", "30 min"    ,  "60 min",     "Gewicht"   ,  "Hypotalamus" ,"Pankreas",  "Gehirn"     , "Gon"      ,   "SCF"     ,    "Leber"      , "Quadrizeps" ,"Longissimus" ,"BAT"       ,  "Herz" ,"LÃ¤nge", "Triglycerides")      
colnames(allPhenotypes) <- colallPhenotypes
write.table(allPhenotypes, "allPhenotypes.txt", sep = "\t", quote = FALSE, row.names = TRUE)
allPhenotypes <- read.csv("allPhenotypes.txt",sep = "\t", header=TRUE, check.names=FALSE)

# MRI analysis
pdf("MRICurve.pdf")
par(mfrow = c(2,1))
timepoints <- as.numeric(colnames(mriLEAN))
plot(c(min(timepoints), max(timepoints)), c(0, max(mriLEAN,na.rm=TRUE)), main= "Lean mass curve", t = 'n', xlab="Time (days)", ylab="Lean mass (grams)")
for(row in 1:nrow(mriLEAN)){
  color <- "lightgreen"
  points(timepoints, mriLEAN[row,], t = 'l', col=color)
}

plot(c(min(timepoints), max(timepoints)), c(0, max(mriLEAN,na.rm=TRUE)), main= "Fat mass curve", t = 'n', xlab="Time (days)", ylab="Fat mass (grams)")
for(row in 1:nrow(mriFAT)){
  color <- "lightblue"
  points(timepoints, mriFAT[row,], t = 'l', col=color)
}
dev.off()
 
## QC by bodyweight
# Bodyweight plot
plot(c(min(days), max(days)), c(0, max(allPhenotypes[,1:31],na.rm=TRUE)), main= "Growth curve", t = 'n', xlab="Time (days)", ylab="Weight (gramms)")
for(row in 1:nrow(allPhenotypes)){
  color <- "lightgreen"
  points(days, allPhenotypes[row,], t = 'l', col=color)
}

# Divide the 2 diet after week 20
highFATnoCARB <- allPhenotypes[,names(allPhenotypes) %in% c("140", "142", "144", "147", "150","154")]
highFAThighCARB <- allPhenotypes[,names(allPhenotypes) %in% c("157", "160", "163", "166", "169", "172", "174")]

op <- par(mfrow=c(2,1))

days_HighFatNoCarb <- as.numeric(colnames(highFATnoCARB))
plot(c(min(days_HighFatNoCarb), max(days_HighFatNoCarb)), c(0, max(highFATnoCARB,na.rm=TRUE)), t = 'n', xlab="time (days)", ylab="Weight (gramms)")
for(x in 1:nrow(highFATnoCARB)){
  color <- rgb(0.5, 0.5, 0.8, 0.5)
  points(days_HighFatNoCarb, highFATnoCARB[x,], t = 'l', col=color)
}

days_HighFatHighCarb <- as.numeric(colnames(highFAThighCARB))
plot(c(min(days_HighFatHighCarb), max(days_HighFatHighCarb)), c(0, max(highFAThighCARB,na.rm=TRUE)), t = 'n', xlab="time (days)", ylab="Weight (gramms)")
for(x in 1:nrow(highFAThighCARB)){
  color <- rgb(0.5, 0.5, 0.8, 0.5)
  points(days_HighFatHighCarb, highFAThighCARB[x,], t = 'l', col=color)
}

diffhighFATnoCARB = (highFATnoCARB[, c("154")] - highFATnoCARB[, c("140")]) / 14
diffhighFAThighCARB = (highFAThighCARB[, c("174")] - highFAThighCARB[, c("157")]) / 17

# Test if there is a difference
t.test(diffhighFATnoCARB, diffhighFAThighCARB, paired = TRUE, var.equal = TRUE)
boxplot(cbind(diffhighFATnoCARB, diffhighFAThighCARB))

# Test if we were allowed to use a t-test
op <- par(mfrow=c(2,1))
hist(diffhighFATnoCARB)
shapiro.test(diffhighFATnoCARB)

hist(diffhighFAThighCARB)
shapiro.test(diffhighFAThighCARB)

# Test if we were allowed to use a t-test
wilcox.test(diffhighFATnoCARB, diffhighFAThighCARB, paired = TRUE)

# Oral Glucose Tolerance test
glucDATA <- allPhenotypes[,36:40]
colnames(glucDATA) <- c("0","15","30","60","120")
time <- colnames(glucDATA)
x <- c(0,15,30,60,120)
plot(main="Oral Glucose Tolerance Test", c(min(x), max(x)), c(0, max(glucDATA[,time], na.rm=TRUE)), t = 'n', xlab="Time (min)", ylab="Blood Glucose (mg/dl)")
for(n in 1:nrow(glucDATA)){
  color <- rgb(0.5, 0.5, 0.8, 0.5)
  points(x, glucDATA[n,], t = 'l', col=color)
}
pdf("oralGTT.pdf")
bpt <- boxplot(glucDATA, main="Oral Glucose Tolerance Test", ylab="Blood Glucose (mg/dl)", xlab="Time(min)", col=c("yellow"), notch=TRUE)
lines(1:5, bpt$stats[ 3, ], col="blue", lwd=2)
dev.off()

# Insulin Tolerance Test
insulinDATA <- allPhenotypes[,41:44]
colnames(insulinDATA) <- c("0","15","30","60")
time <- colnames(insulinDATA)
x <- c(0,15,30,60)
plot(main="Insulin Tolerance Test", c(min(x), max(x)), c(0, max(insulinDATA[,time], na.rm=TRUE)), t = 'n', xlab="Time (min)", ylab="Blood Glucose (mg/dl)")
for(n in 1:nrow(insulinDATA)){
  color <- rgb(0.5, 0.5, 0.8, 0.5)
  points(x, insulinDATA[n,], t = 'l', col=color)
}
pdf("insulinTest.pdf")
bptt <- boxplot(insulinDATA[,1], as.numeric(as.character(insulinDATA[,2])), as.numeric(as.character(insulinDATA[,3])), as.numeric(as.character(insulinDATA[,4])), main="Insulin Tolerance Test", ylab="Blood Glucose (mg/dl)", xlab="Time(min)", col=c("green"), notch=TRUE)
axis(1, at = 1:4 , c("0", "15", "30", "60"))
lines(1:4, bptt$stats[3, ], col="blue", lwd=2)
dev.off()

# Difference in weight between gonadal adipose tissue and liver
allTissues <- allTissues[-which(apply(apply(allTissues,1,is.na),2,sum) > 0),]
allTissues[33,"Leber"] <- 10 * allTissues[33,"Leber"]
rownames(allTissues) <- 1:nrow(allTissues)

# Order columns by gon weight
ordering <- sort(allTissues[,"Gon"], index.return=TRUE)$ix
sortWeight <- allTissues[ordering,]

# Adjust tissues weight by total weight
sortWeight <- sweep(sortWeight,sortWeight[,"Gewicht"],MARGIN=1,"/")

plot(main="Weight relationship", c(1, nrow(sortWeight)), c(0, max(sortWeight[,"Gon"], na.rm=TRUE)), t = "n", xlab="Individuals", ylab="Weight (adjust)")
lines(sortWeight[,"Gon"], col = "blue" , lwd=2 , pch=19 , type="l")
lines(sortWeight[,"Leber"], col = "orange" , lwd=2 , pch=19 , type="l")
lines(sortWeight[,"SCF"], col = "green" , lwd=2 , pch=19 , type="l")
 legend("topleft",
  legend = c("Liver", "Gon", "SCF"),
  col = c("orange", "blue", "green"),
  pch = c(20,20,20),
  bty = "n",
  pt.cex = 2,
  cex = 1.2,
  text.col = "black",
  #horiz = F ,
  #inset = c(0.1, 0.1, 0.1)
)

mmodel <- lm(allTissues[,"Gon"] ~ allTissues[,"Leber"])
abline(a = mmodel$coefficients["(Intercept)"], b = mmodel$coefficients["allTissues[, \"Leber\"]"])

# Adjust by weight plus SCF, Longissimus
plot(main="Weight relationship between gonadal adipose tissue and liver (adjust)", c(1, nrow(sortWeight)), c(0, max(sortWeight[,"Gon"], na.rm=TRUE)), t = "n", xlab="Individuals", ylab="Weight (grams)")
lines(sortWeight[,"Gon"], col = "red" , lwd=1 , pch=20 , type="l")
lines(sortWeight[,"Leber"], col = "blue" , lwd=1 , pch=20 , type="l")
lines(sortWeight[,"sCF"], col = "green" , lwd=1 , pch=20 , type="l")
lines(sortWeight[,"longiss"], col = "black" , lwd=1 , pch=20 , type="l")
lines(sortWeight[,"brain"], col = "brown" , lwd=1 , pch=20 , type="l")
  legend("topleft",
   legend = c("Liver", "Gon", "SCF", "longissimus", "brain"),
   col = c("blue", "red", "green", "black", "brown"),
   pch = c(20,20,20),
   bty = "n",
   pt.cex = 2,
   cex = 1.2,
   text.col = "black",
   horiz = F ,
   inset = c(0.1, 0.1, 0.1))

#extremes for FTIR
sortWeight <- sortWeight[c(1:10,267:277),]

# Create a big matrix with all the phenotypes
phenotypes <- cbind(allPhenotypes, mriLEAN, mriFAT, insulinDATA, glucDATA, allTissues)
