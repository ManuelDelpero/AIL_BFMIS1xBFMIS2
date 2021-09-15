# 
#
# copyright (c) 2018-2021 - Brockmann group - HU Berlin Manuel Delpero
# 
# 

setwd("C:/Users/Manuel/Desktop/AIL_S1xS2/RAWDATA")

data <- read.csv("qPCR_Complete_Liver.txt", header = TRUE, check.names = FALSE, sep="\t")
#data <- read.csv("qPCR_complete_FatGapdh.txt", header = TRUE, check.names = FALSE, sep="\t")

data[,"dCT"] <- data[,"meanG"] - data[,"meanH"]
data[,"Transf"] <- 2^-(data[,"dCT"])
genes <- unique(data[,"Target Name"])

results <- matrix(NA, length(as.character(unique(data[,"Target Name"]))), 6,)
colnames(results) <- c("meanS1_relExpr", "meanS2_relExpr", "sdS1_relExpr", "sdS2_relExpr", "2^-(ddCT)", "p-value")
rownames(results) <- as.character(unique(data[,"Target Name"]))
toplot <- c()
toplotSDs <- c()
for (gene in as.character(genes)){
  currentG <- data[grep(gene, data[,"Target Name"]),]
  p <- t.test(currentG[grep("S1", currentG[,"Sample Name"]), "Transf"], currentG[grep("S2", currentG[,"Sample Name"]), "Transf"])[[3]]
  meanS1 <- mean(currentG[grep("S1", currentG[,"Sample Name"]), "dCT"], na.rm = TRUE)
  meanS2 <- mean(currentG[grep("S2", currentG[,"Sample Name"]), "dCT"], na.rm = TRUE)
  sdS1 <- sd(currentG[grep("S1", currentG[,"Sample Name"]), "dCT"], na.rm = TRUE)
  sdS2 <- sd(currentG[grep("S2", currentG[,"Sample Name"]), "dCT"], na.rm = TRUE)
  relativeExprS1 <- 2^-(meanS1)
  relativeExprS2 <- 2^-(meanS2)
  relativeExprS1SD <- sd(currentG[grep("S1", currentG[,"Sample Name"]), "Transf"]/ mean(currentG[grep("S2", currentG[,"Sample Name"]), "Transf"], na.rm = TRUE))
  relativeExprS2SD <- sd(currentG[grep("S2", currentG[,"Sample Name"]), "Transf"]/ mean(currentG[grep("S2", currentG[,"Sample Name"]), "Transf"], na.rm = TRUE))
  fold <- 2^-(meanS1-meanS2)
  results[gene,"2^-(ddCT)"] <- fold
  results[gene,"meanS2_relExpr"] <- 1
  results[gene,"sdS1_relExpr"] <- relativeExprS1SD
  results[gene,"sdS2_relExpr"] <- relativeExprS2SD 
  results[gene,"p-value"] <- p
  toplot <- rbind(toplot, relativeExprS1/relativeExprS2,1) 
  toplotSDs <- rbind(toplotSDs, relativeExprS1SD,relativeExprS2SD) 
}
  
#barplot
par(cex.lab=1.2, cex.main = 1.8, cex.axis = 1.4)
par(mfrow = c(2,2))
#par(mfrow = c(1,2))
#fat
x <- barplot(toplot[,1],beside=TRUE,col= c("gray88", "gray20"), main = "", ylab = expression(paste("Relative expression [2ˆ(- ",Delta,Delta,"CT)]")), ylim = c(0, 4.5), space=c(0.2,0 , 0.5,0, 0.5,0, 0.5,0, 0.5,0 ,0.5,0, 0.5,0))
  axis(1, at = c(1.2, 3.7, 6.2, 8.7, 11.2, 13.7, 16.2), c(expression(italic("Fmo5")), expression(italic("Notch2")), expression(italic("Acat2")), expression(italic("Trappc9")), expression(italic("Zfat")), expression(italic("Rrn3")), expression(italic("Trap1"))))
  y <- 2.5
  # set an offset for tick lengths
  offset <- 0.1
  # draw first horizontal line
  lines(x[1:2],c(y, y))
  lines(x[5:6],c(y, y))
  lines(x[7:8],c(y, y))
  lines(x[9:10],c(y , y ))
  lines(x[13:14],c(y, y))
  lines(x[3:4],c(y, y))
  lines(x[11:12],c(y, y))
  # draw asterics
  text(1.2,y +offset,"***")
  text(6.2,y+offset,"***")
  text(8.7,y+offset,"***")
  text(11.2,y+offset,"*")
  text(16.2,y+offset,"***")
  text(3.7,y+offset,"**")
  text(13.7, y+offset,"*")
  legend("topleft", legend=c("S1", "S2"), fill= c("gray88","gray20"), bty = "n")
  for (row in 1:nrow(toplot)){
    segments(x[row,1], toplot[row,1] - as.numeric(toplotSDs[row,1]), x[row,1], toplot[row,1] + as.numeric(toplotSDs[row,1]), lwd = 1.5)
    arrows(x[row,1], toplot[row,1] - as.numeric(toplotSDs[row,1]), x[row,1], toplot[row,1] + as.numeric(toplotSDs[row,1]), lwd = 1.5, angle = 90, code = 3, length = 0.05)
  }
  
#liver
x <- barplot(toplot[,1],beside=TRUE,col= c("gray88", "gray20"), main = "", ylab = expression(paste("Relative expression [2ˆ(- ",Delta,Delta,"CT)]")), ylim = c(0, 4.5), space=c(0.2,0 , 0.5,0, 0.5,0, 0.5,0, 0.5,0 ,0.5,0, 0.5,0))
  axis(1, at = c(1.2, 3.7, 6.2, 8.7, 11.2, 13.7, 16.2), c(expression(italic("Fmo5")), expression(italic("Notch2")), expression(italic("Acat2")), expression(italic("Trappc9")), expression(italic("Zfat")), expression(italic("Rrn3")), expression(italic("Trap1"))))
  y <- 2.5
  # set an offset for tick lengths
  offset <- 0.1
  # draw first horizontal line
  lines(x[1:2],c(y, y))
  lines(x[5:6],c(y, y))
  lines(x[7:8],c(y, y))
  lines(x[9:10],c(y , y ))
  lines(x[13:14],c(y, y))
  lines(x[3:4],c(y, y))
  lines(x[11:12],c(y, y))
  # draw asterics
  text(1.2,y +offset,"***")
  text(6.2,y+offset+offset,"p = 0.88")
  text(8.7,y+offset,"***")
  text(11.2,y+offset+ offset,"p = 0.06")
  text(16.2,y+offset+ offset,"p = 0.96")
  text(3.7,y+offset,"**")
  text(13.7, y+offset+ offset,"p = 0.07")
  legend("topleft", legend=c("S1", "S2"), fill= c("gray88","gray20"), bty = "n")
  for (row in 1:nrow(toplot)){
    segments(x[row,1], toplot[row,1] - as.numeric(toplotSDs[row,1]), x[row,1], toplot[row,1] + as.numeric(toplotSDs[row,1]), lwd = 1.5)
    arrows(x[row,1], toplot[row,1] - as.numeric(toplotSDs[row,1]), x[row,1], toplot[row,1] + as.numeric(toplotSDs[row,1]), lwd = 1.5, angle = 90, code = 3, length = 0.05)
  }