#' @title CCThreshold
#' @description Averages multiple Fullmaps together, by averaging HpM.
#' @param path.name A character string defining folder containing individual FullMaps.
#' @param exp.name Manually specify a name for the outputted pool file. Defaults to a concatenation of individual file names.
#' @param out.mode output in "txt", "rds", or "both" mode.
#' @examples
#' @author Will Gittens
#' @import data.table
#' @export
CCThreshold<-function(path.name=getwd(),threshold=1,readfilename="Fullmap."){
parent.directory <- dirname(path.name)
  setwd(parent.directory)
library("data.table")
files <- list.files(pattern =readfilename )


CCList(path.name) # Calculate Million reads per sample for conveting to HpM

# on a per strand basis but HpMillion mapped reads total
DSBListm <- mapply(function(x,y,z) {
  x$Crick <- x$Crick/y
  x$Watson <- x$Watson/y
    return(x)
}, DSBList, Mreads, SIMPLIFY = F)


DSBList2 <- lapply(DSBListm, function(x){
  x <- x[x$Crick > threshold | x$Watson > threshold,]
  x[, 3:4][x[, 3:4]<threshold] <- 0L # THIS LINE IS EXTREMELY IMPORTANT - WITHOUT IT YOU SEE A STRONG ARTIFACT IN THE OFFSET ANALYSIS AT W-C = 0.
  return(x)
})

# on a per strand basis but HpMillion mapped reads total
DSBList2 <- mapply(function(x,y,z) {
  x$Crick <- x$Crick*y
  x$Watson <- x$Watson*y
  return(x)
}, DSBList2, Mreads, SIMPLIFY = F)

Msites <- sapply(DSBList2, function(x) {nrow(x)/1000000})
Msites.before.threshold <- sapply(DSBListm, function(x) {
  nrow(x)/1000000})
Mreads.before.threshold <- Mreads

Mreads <- sapply(DSBList2, function(x) {sum(as.numeric(x$Watson) + as.numeric(x$Crick))/1000000})
perc.Mreads.removed.by.threshold <- ((Mreads.before.threshold - Mreads)/Mreads.before.threshold)*100
perc.Msites.removed.by.threshold <- ((Msites.before.threshold - Msites)/Msites.before.threshold)*100

filter.metrics <- data.frame(cbind(Mreads.before.threshold, Mreads, perc.Mreads.removed.by.threshold, Msites.before.threshold, Msites, perc.Msites.removed.by.threshold))

assign(paste0("DSBList.", threshold, "HpM"), DSBList2)

ifelse(!dir.exists(file.path(parent.directory, paste0("TH",threshold, "HpM"))), dir.create(file.path(parent.directory, paste0("TH",threshold, "HpM"))), FALSE); setwd(file.path(parent.directory, paste0("TH",threshold, "HpM")))
ifelse(!dir.exists(file.path(parent.directory, paste0("TH",threshold, "HpM"), "FullMaps")), dir.create(file.path(parent.directory, paste0("TH",threshold, "HpM"), "FullMaps")), FALSE); setwd(file.path(parent.directory, paste0("TH",threshold, "HpM"), "FullMaps"))
CCWrite(DSBList2, filenames = paste0(DSBListNames, "_TH", threshold, "HpM.txt" ))
}
