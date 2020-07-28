setwd("C:/Users/Manuel/Desktop/AIL_S1xS2")

mdata <- read.table("outputVEPMGP.vcf", header = TRUE)

gts <- apply(mdata[,10:ncol(mdata)],2,function(x){unlist(lapply(strsplit(x, ":"),"[",1)) })
gts[gts=="./."] <- NA

vepp <- unlist(lapply(strsplit(as.character(mdata[,8]), ";"), "[",15))
veppsplit <- strsplit(gsub("CSQ=", "", vepp), ",")
rsIDs <- unlist(lapply(strsplit(unlist(lapply(veppsplit, "[", 1)), "|",fixed=TRUE), "[",18))
types <- unlist(lapply(lapply(veppsplit, strsplit, "|", fixed=TRUE), function(x){ 
  type <- unlist(lapply(x,"[",2))
  gene <- unlist(lapply(x,"[",4))
  transcript <- unlist(lapply(x,"[",7))
  for(x in 1:length(type)){
    if(type[x] %in% c("missense_variant")) { return(paste0(type[x], "\t", gene[x], "")) }
    if(type[x] %in% c("regulatory_region_variant")) { return(paste0("regulatory\t", transcript[x], "")) }
    if(type[x] %in% c("splice_acceptor_variant")) { return(paste0("splice_acceptor\t", gene[x], "")) }
    if(type[x] %in% c("splice_donor_variant")) { return(paste0("splice_donor\t", gene[x], "")) }
    if(type[x] %in% c("3_prime_UTR_variant")) { return(paste0("3_prime\t", gene[x], "")) }
    if(type[x] %in% c("5_prime_UTR_variant")) { return(paste0("5_prime\t", gene[x], "")) }
    if(type[x] %in% c("synonymous_variant")) { return(paste0("synonymous\t", gene[x], "")) }
    if(type[x] %in% c("stop_gained")) { return(paste0("stop_gained\t", gene[x], "")) }
    if(type[x] %in% c("stop_lost")) { return(paste0("stop_lost\t", gene[x], "")) }
    if(type[x] %in% c("start_lost")) { return(paste0("start_lost\t", gene[x], "")) }
  }
  return("\t")
}))

mexcel <- cbind("RSID" = rsIDs, mdata[, c(1,2,4,5,6)], "TYPE" = types, gts)

# Only SNPs that are found in the BFMI
mexcel <- mexcel[which(!(mexcel[, "BFMI861.S1"] == "0/0" & mexcel[, "BFMI861.S2"] == "0/0" & mexcel[, "BFMI860.12"] == "0/0")),]

write.table(mexcel, file="MGP_Chr3_36000000-37000000.snps-filtered.tsv", sep = "\t", quote = FALSE, row.names=FALSE, na = "")
