# Effect plots for each qtl
par(cex.lab=2, cex.main = 2, cex.axis = 2)
mat <- matrix(c(1,2,3), 1, ,byrow = TRUE)
layout(mat, widths = rep.int(3, ncol(mat)))
#chr3
UNC5791802 <- cbind(phenotypes[, "Gon"], t(genotypes["UNC5812781",]))
boxplot(as.numeric(UNC5791802[which(UNC5791802[,2] == "A"),1]), as.numeric(UNC5791802[which(UNC5791802[,2] == "H"),1]), as.numeric(UNC5791802[which(UNC5791802[,2] == "B"),1]), col = c("#2C86CA", "#ADCCF6", "#F9F9F9"), main = "GonAT weight",  ylab = "Weight [g]", xlab = "UNC5812781 [chr 3]" , las = 2, t = "n", xaxt = "n", ylim = c(0,4))
  axis(1, at = 1:3 , c("S1", "HET", "S2"))
  y1 <- 7
  #lines(c(1,3),c(y1,y1))
  #text(2, y1+0.3,"***", cex = 1.4)
  #legend("topright", #bg="gray",
  #legend = c( "S1", "S1/S2", "S2"),
  #  col = c(""#2C86CA", "#ADCCF6", "#F9F9F9""),
  #  pch = 15,
  #  pt.cex = 1.7,
  #  pt.bg = "lightsteelblue1",
  #  cex = 1.2,
  #  text.col = "black") 

UNC5791802 <- cbind(phenotypes[, "Gluc166"], t(genotypes["UNC5812781",]))
boxplot(as.numeric(UNC5791802[which(UNC5791802[,2] == "B"),1]), as.numeric(UNC5791802[which(UNC5791802[,2] == "H"),1]), as.numeric(UNC5791802[which(UNC5791802[,2] == "A"),1]), col = c("#2C86CA", "#ADCCF6", "#F9F9F9"), main = "Blood glucose", ylab = "", xlab = "UNC5812781 [chr 3]" , las = 2, t = "n", xaxt = "n", ylim = c(0,550))
  axis(1, at = 1:3 , c("S1", "HET", "S2"))
  y1 <- 550
  #lines(c(1,3),c(y1,y1))
  #lines(c(1,2),c(y1+30,y1+30))
  #text(2, y1+10,"***", cex = 1.4)
  #text(1.5, y1+40,"***", cex = 1.4)
  #legend("topright", #bg="gray",
  #legend = c( "S1", "S1/S2", "S2"),
  #  col = c(""#2C86CA", "#ADCCF6", "#F9F9F9""),
  #  pch = 15,
  #  pt.cex = 1.7,
  #  #pt.bg = "lightsteelblue1",
  #  cex = 1.2,
  #  text.col = "black") 
	

	
#chr15
UNC25805470 <- cbind(phenotypes[, "Gon"], t(genotypes["UNC25805470",]))
boxplot(as.numeric(UNC25805470[which(UNC25805470[,2] == "A"),1]), as.numeric(UNC25805470[which(UNC25805470[,2] == "H"),1]), as.numeric(UNC25805470[which(UNC25805470[,2] == "B"),1]), col = c("#2C86CA", "#ADCCF6", "#F9F9F9"), main = "gWAT weight", ylab = "Weight [g]", xlab = "UNC25805470 [chr 15]" , las = 2, t = "n", xaxt = "n", ylim = c(0,4))
  axis(1, at = 1:3 , c("S1", "HET", "S2"))
  y1 <- 7
  #lines(c(1,3),c(y1,y1))
  #text(2, y1+0.3,"***", cex = 1.4)
  #legend("topright", #bg="gray",
  #legend = c( "S1", "S1/S2", "S2"),
  #  col = c(""#2C86CA", "#ADCCF6", "#F9F9F9""),
  #  pch = 15,
  #  pt.cex = 1.7,
  #  pt.bg = "lightsteelblue1",
  #  cex = 1.2,
  #  text.col = "black") 
  
  
# chr 15 glucose effect plot with adjusted values
numgeno <- matrix(NA, nrow(genotypes), ncol(genotypes), dimnames=list(rownames(genotypes), colnames(genotypes)))
for(x in 1:nrow(genotypes)){
  h1 <- "A"
  het <- "H"
  h2 <- "B"
  numgeno[x, which(genotypes[x, ] == h1)] <- -1
  numgeno[x, which(genotypes[x, ] == het)] <- 0
  numgeno[x, which(genotypes[x, ] == h2)] <- 1
}
	
phenos <- phenotypes[colnames(genotypes),]
phenos[,"Gluc166"] <- as.numeric(phenos[,"Gluc166"])
phenos[,"Leber"] <- as.numeric(phenos[,"Leber"])
phenos[,"Gon"] <- as.numeric(phenos[,"Gon"])

JAX00063853 <- as.numeric(unlist(numgeno["JAX00063853",]))
mdata <- data.frame(cbind(glucose = phenos[, "Gluc172"], liver = phenos[, "Leber"], gon = phenos[, "Gon"], geno = JAX00063853))
isNA <- which(apply(apply(mdata,1,is.na),2,any))
mdata <- mdata[-isNA,]
lmGluc <- lm(glucose ~ liver + gon, data = mdata)
Gluc166.adj <- mean(mdata[,"glucose"], na.rm = TRUE) + resid(lmGluc)
mdata <- cbind(mdata, Gluc166.adj)

for(x in 1:nrow(mdata)){
  h1 <- -1
  het <- 0
  h2 <- 1
  mdata[x, which(mdata[x, ] == h1)] <- "A"
  mdata[x, which(mdata[x, ] == het)] <- "H"
  mdata[x, which(mdata[x, ] == h2)] <- "B"
}

boxplot(as.numeric(mdata[which(mdata[,"geno"] == "A"),5]), as.numeric(mdata[which(mdata[,"geno"] == "H"),5]), as.numeric(mdata[which(mdata[,"geno"] == "B"),5]), col = c("#2C86CA", "#ADCCF6", "#F9F9F9"), main = "Blood glucose", ylab = "", xlab = "JAX00063853 [chr 15]" , las = 2, t = "n", xaxt = "n", ylim = c(0,420))
  axis(1, at = 1:3 , c("S1", "HET", "S2"))
  y1 <- 200
  #lines(c(1,3),c(y1,y1))
  #lines(c(2,3),c(y1+30,y1+30))
  #text(2, y1+10,"***", cex = 1.4)
  #text(2.5, y1+40,"***", cex = 1.4)
  #legend("topright", #bg="gray",
  #legend = c( "S1", "HET", "S2"),
  #  col = c("#2C86CA", "#ADCCF6", "#F9F9F9"),
  #  pch = 15,
  #  pt.cex = 1.7,
  #  pt.bg = "lightsteelblue1",
  #  cex = 1.2,
  #  text.col = "black")
	
#chr17
UNC27568354 <- cbind(phenotypes[, "Gon"], t(genotypes["UNCHS043909",]))
boxplot(as.numeric(UNC27568354[which(UNC27568354[,2] == "A"),1]), as.numeric(UNC27568354[which(UNC27568354[,2] == "H"),1]), as.numeric(UNC27568354[which(UNC27568354[,2] == "B"),1]), col = c("#2C86CA", "#ADCCF6", "#F9F9F9"), main = "gWAT weight", ylab = "Weight [g]", xlab = "UNC27568354 [chr 17]" , las = 2, t = "n", xaxt = "n", ylim = c(0,4))
  axis(1, at = 1:3 , c("S1", "HET", "S2"))
  y1 <- 7
  #lines(c(1,3),c(y1,y1))
  #lines(c(2,3),c(y1+0.7,y1+0.7))
  #text(2, y1+0.3,"***", cex = 1.4)
 # text(2.5, y1+1,"***", cex = 1.4)
  #legend("topright", #bg="gray",
  #legend = c( "S1", "S1/S2", "S2"),
  #  col = c("#2C86CA", "#ADCCF6", "#F9F9F9"),
  #  pch = 15,
  #  pt.cex = 1.7,
  #  pt.bg = "lightsteelblue1",
  #  cex = 1.2,
  #  text.col = "black") 

UNCHS043909 <- cbind(phenotypes[, "Gluc166"], t(genotypes["UNCHS043909",]))
boxplot(as.numeric(UNCHS043909[which(UNCHS043909[,2] == "A"),1]), as.numeric(UNCHS043909[which(UNCHS043909[,2] == "H"),1]), as.numeric(UNCHS043909[which(UNCHS043909[,2] == "B"),1]), col = c("#2C86CA", "#ADCCF6", "#F9F9F9"), main = "Blood glucose", ylab ="Glucose [mg/dl]", xlab = "UNCHS043909 [chr 17]" , las = 2, t = "n", yaxt = "n", xaxt = "n", ylim = c(0,550))
  axis(1, at = 1:3 , c("S1", "HET", "S2"))
  y1 <- 550
  #lines(c(1,3),c(y1,y1))
  #lines(c(1,2),c(y1+30,y1+30))
 # text(2, y1+10,"***", cex = 1.4)
 # text(1.5, y1+40,"***", cex = 1.4)
 # legend("topright", #bg="gray",
 # legend = c( "S1", "S1/S2", "S2"),
 #   col = c("#2C86CA", "#ADCCF6", "#F9F9F9"),
 #   pch = 15,
 #   pt.cex = 1.7,
 #   pt.bg = "lightsteelblue1",
    #cex = 1.2,
    #text.col = "black") 
	
#chr17
UNCHS043909 <- cbind(phenotypes[, "Leber"], t(genotypes["UNCHS043909",]))
boxplot(as.numeric(UNCHS043909[which(UNCHS043909[,2] == "A"),1]), as.numeric(UNCHS043909[which(UNCHS043909[,2] == "H"),1]), as.numeric(UNCHS043909[which(UNCHS043909[,2] == "B"),1]), col = c("#2C86CA", "#ADCCF6", "#F9F9F9"), main = "Liver weight", ylab = "Weight [g]", xlab = "UNCHS043909 [chr 17]" , las = 2, t = "n", xaxt = "n", ylim = c(0,5))
  axis(1, at = 1:3 , c("S1", "HET", "S2"))
  y1 <- 6
  #lines(c(1,3),c(y1,y1))
  #lines(c(2,3),c(y1+0.7,y1+0.7))
  #text(2, y1+0.3,"***", cex = 1.4)
  #text(2.5, y1+1,"***", cex = 1.4)
  #legend("topright", #bg="gray",
  legend = c( "S1", "S1/S2", "S2"),
    col = c("#2C86CA", "#ADCCF6", "#F9F9F9"),
    pch = 15,
    pt.cex = 1.7,
    pt.bg = "lightsteelblue1",
    cex = 1.2,
    text.col = "black") 
	
#chr17
UNCHS043909 <- cbind(phenotypes[, "Triglycerides"], t(genotypes["JAX00632487",]))
boxplot(as.numeric(UNCHS043909[which(UNCHS043909[,2] == "B"),1]), as.numeric(UNCHS043909[which(UNCHS043909[,2] == "H"),1]), as.numeric(UNCHS043909[which(UNCHS043909[,2] == "A"),1]), col = c("#2C86CA", "#ADCCF6", "#F9F9F9"), main = "Liver triglycerides", ylab = "[ug/ug]", xlab = "JAX00632487 [chr 7]" , las = 2, t = "n", xaxt = "n", ylim = c(0,350))
  axis(1, at = 1:3 , c("S1", "HET", "S2"))
  #lines(c(1,3),c(y1,y1))
  #lines(c(2,3),c(y1+0.7,y1+0.7))
  #text(2, y1+0.3,"***", cex = 1.4)
  #text(2.5, y1+1,"***", cex = 1.4)
  #legend("topright", #bg="gray",
  legend = c( "S1", "S1/S2", "S2"),
    col = c("#2C86CA", "#ADCCF6", "#F9F9F9"),
    pch = 15,
    pt.cex = 1.7,
    pt.bg = "lightsteelblue1",
    cex = 1.2,
    text.col = "black")

# Chr16 
UNC5791802 <- cbind(phenotypes[, "D174"], t(genotypes["UNCHS041714",]))
boxplot(as.numeric(UNC5791802[which(UNC5791802[,2] == "A"),1]), as.numeric(UNC5791802[which(UNC5791802[,2] == "H"),1]), as.numeric(UNC5791802[which(UNC5791802[,2] == "B"),1]), col = c("#2C86CA", "#ADCCF6", "#F9F9F9"), main = "Body weight", ylab = "Body weight [g]", xlab = "UNCHS041714 [chr 16]" , las = 2, t = "n", xaxt = "n", yaxt = "n", ylim = c(30,60))
  axis(1, at = 1:3 , c("S1", "HET", "S2"))
  y1 <- 550
  #lines(c(1,3),c(y1,y1))
  #lines(c(1,2),c(y1+30,y1+30))
  #text(2, y1+10,"***", cex = 1.4)
  #text(1.5, y1+40,"***", cex = 1.4)
  #legend("topright", #bg="gray",
  #legend = c( "S1", "S1/S2", "S2"),
  #  col = c("lightskyblue1", "cyan3", "dodgerblue4"),
  #  pch = 15,
  #  pt.cex = 1.7,
  #  #pt.bg = "lightsteelblue1",
  #  cex = 1.2,
  #  text.col = "black") 	
	
	
plot(main = "Chr 7", c(min(as.numeric(datasetRow[, "Position"])), max(as.numeric(datasetRow[, "Position"]))), c(-0.5,6), ylab = "LOD score", xlab = "Position [mb]", las = 2, t = "n", xaxt = "n")
  points(x = as.numeric(dataset[,"Position"]), y = dataset[,"Gon"], t ='l', col="cornflowerblue", lwd = 2)
  points(x = as.numeric(datasetRow[,"Position"]), y = datasetRow[,"Triglycerides"], t ='l', col="gray32", lwd = 2)
  #lines(x = as.numeric(dataset[,"Position"]), y = dataset[,"Gluc166"], t = "l", lty = 2, col="red", lwd = 0.7)
  lines(x = as.numeric(datasetRow[,"Position"]), y = datasetRow[,"Gluc166"], t = "l", lty = 1, col="red", lwd = 2)
  points(x = as.numeric(datasetRow[,"Position"]), y = rep(-0.5,nrow(datasetRow)), pch = "|", col = "black", lwd = 1)
  abline(h=4.7, col="orange" )
  abline(h=4.2, col="orange", lty = 2)
  axis(1, at = c(10000000, 30000000, 60000000, 90000000, 120000000, 150000000, 180000000 ), c("10", "30", "60", "90", "120", "150", "180" ))
  legend("topright",
  legend = c("GonAT weight = SNP + error", "BGc = SNP + error", "BGc [adj] = liver weight + GonAT weight + SNP + error", "Liver weight = litter size + SNP + error"),
    bty = "n",
    col = c("cornflowerblue", "red", "red", "gray32" ),
    lty=c(1,1,2,1),
    pt.cex = 1.2,
    pt.bg = "lightsteelblue1",
    cex = 1.2,
    text.col = "black")
	
	
