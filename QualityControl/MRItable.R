setwd("C:/Users/Manuel/Desktop/AIL_S1xS2/RAWDATA/Data")

MRIdata <- read.csv("mri.txt", header = TRUE, sep= "\t", check.names = FALSE)
description <- read.csv("descript.txt", header = TRUE, sep= "\t", check.names = FALSE)
IDs <- read.csv("IDs.txt", header = TRUE, sep= "\t", check.names = FALSE)
IDs <- gsub("V 888-", "", IDs[,1])
rownames(description) <- IDs
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
months <- data.frame(months)

dates <- strsplit(unlist(lapply(strsplit(as.character(MRIdata[,"TimeDateDura"]),";"),"[",1))," ")                  				  	    # Get the malformed dates
MRIdata[,"TimeDateDura"] <- unlist(lapply(dates, function(x){                                                           			 	# Transform to DD/MM/YY
	monthNumber <- which(months[,1]==x[2])
	paste(gsub(",","",x[3]), monthNumber, x[4],sep="/")
}))




for(x in 1:nrow(MRIdata)){
    bDay <- as.character(description[as.character(MRIdata[x,"Label"]),"W-dat"]) 														# Date of birth
	mDay <- as.character(MRIdata[x, "TimeDateDura"]) 																					# Date of measurement
	mDay <- gsub("2018", "18", mDay)
	mDay <- gsub("2019", "19", mDay)
	daysDiff <- as.numeric(round(difftime(strptime(mDay, format = "%d/%m/%y"), strptime(bDay, format = "%d/%m/%y"), units="days")))
	if(length(daysDiff) == 0) daysDiff<- 666                                                                                            # If one of the dates is missing use 666
	cat(mDay, bDay, daysDiff, "\n")
	MRIdata[x, "Age"] <- daysDiff
}

animals <- unique(as.character(MRIdata[,"Label"]))
timepoints <- unique(as.character(MRIdata[,"Age"]))
timepoints <- timepoints[-c(2)]


fat <- matrix(NA, length(animals), length(timepoints), dimnames=list(animals, timepoints))
lean <- matrix(NA, length(animals), length(timepoints), dimnames=list(animals, timepoints))

 for(tp in timepoints){
   for(animal in animals){
      ii <- which(MRIdata[,"Label"] == animal & MRIdata[,"Age"] == tp)
      if(length(ii) > 0){
        fat[animal, tp] <- mean(MRIdata[ii,"Fat"])
        lean[animal, tp] <- mean(MRIdata[ii,"Lean"])
     }
   }
 }

