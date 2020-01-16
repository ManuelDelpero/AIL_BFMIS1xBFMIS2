# Peak detections version # 2 with 5MB look-ahead
#
# copyright (c) 2015-2020 - Brockmann group - HU Berlin, Danny Arends & Manuel DelPero
# last modified May, 2019
# first written May, 2019

peak.detect <- function(profile, map, cutoff = 4, loddrop = 1.5){
  peakstart <- c()
  peakend <- c()
  peaks <- c()
  while(!all(is.na(profile)) && max(profile, na.rm=TRUE) > cutoff){
    topM <- which.max(profile)
    cat("topMarker:", topM, "\n")
    right <- fwalk(profile, map, loddrop, topM, 1)
    left <- fwalk(profile, map, loddrop, topM, -1)
    peakstart <- c(peakstart, names(profile)[left])
    peakend <- c(peakend, names(profile)[right])
    peaks <- rbind(peaks, c(left, topM, right))
    leftRange <- getMarkersWithinRange(map, left, -1)
    rightRange <- getMarkersWithinRange(map, right, 1)
    profile[min(leftRange):max(rightRange)] <- -1
  }
  return(peaks)
}

fwalk <- function(profile, map, loddrop, idx, dir = 1){
  cat("IDX", idx, "Dir", dir, "\n")
  topChr <- map[idx, "Chromosome"]
  n <- 1
  bottom <- (profile[idx] - loddrop)
  within5mb <- getMarkersWithinRange(map, idx + (n * dir), dir)
  while(profile[idx + (n * dir)] > bottom || any(profile[within5mb] > bottom)){
    if(map[idx + ((n + 1) * dir), "Chromosome"] != topChr) break;
    n <- n + 1
    within5mb <- getMarkersWithinRange(map, idx + (n * dir), dir)
    if(idx + (n * dir) == 1 || idx + (n * dir) == (length(profile))) return(idx + (n * dir));
  }
  return(idx + (n * dir))
}

getMarkersWithinRange <- function(map, idx, dir, window = 5000000){
  curChr <- map[idx, "Chromosome"]
  curPos <- map[idx, "Position"]
  lookTo <- (curPos + dir * window)
  withinrange <- which(map[,"Chromosome"] == curChr & map[,"Position"] %in% curPos:lookTo)
  return(withinrange)
}
