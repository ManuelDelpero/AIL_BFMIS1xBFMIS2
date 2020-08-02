# AIL_S1xS2 Analysis on Phenotypes
#
# copyright (c) - Manuel Delpero
# first written july, 2020
# 

setwd("C:/Users/Manuel/Desktop/AIL_S1xS2/RAWDATA")

pheno <- read.table("allPhenotypes_final.txt", sep = "\t", row.names=1)
# Order columns by gon weight
index <- order(pheno[,"Gon"], decreasing = FALSE)
sortWeight <- pheno[index,]

# Plot to represent the correlation between the tissues weight, add the correlation values!
plot(main="Weight relationship between gonadal adipose tissue and liver (adjust)", c(1, nrow(sortWeight)), c(0, max(sortWeight[,"Gon"], na.rm=TRUE)), t = "n", xlab="Individuals", ylab="Weight (grams)", ylim = c(0,0.13))
lines(sortWeight[,"Gon"], col = "red" , lwd=1 , pch=20 , type="l")
lines(sortWeight[,"Leber"], col = "blue" , lwd=1 , pch=20 , type="l")
lines(sortWeight[,"SCF"]/sortWeight[, "Gewicht"], col = "green" , lwd=1 , pch=20 , type="l")
   legend("topleft",
   legend = c("Liver", "Gon", "SCF"),
   col = c("blue", "red", "green"),
   pch = c(20,20,20),
   bty = "n",
   pt.cex = 2,
   cex = 1.2,)
   
# CTL mapping curve across chromosome 15




