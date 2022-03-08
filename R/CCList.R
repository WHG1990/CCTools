#' @title CCList
#' @description Reads a folder of FullMaps into a DSBList, calculates Mreads and creates the DSBListNames.
#' @param path.name A character string defining the folder containing FullMaps.
#' @examples
#' @author Will Gittens
#' @import data.table
#' @export
CCList <- function(path.name){
  library("data.table")
setwd(path.name)
DSBList=list();DSBListNames=NULL
#Read in all tables with string "Full.Map."
files = list.files(pattern="FullMap.") # import files names with "FullMap." string into variable "files"
DSBListNames <- substr(files, 9, nchar(files)-4) # Shorten filename by 8 characters from beginning and 6 characters form end (i.e. remove "FullMap." and "_c.txt")
for(i in 1:length(files)){
  DSBList[[i]] <- fread(files[i], sep = "\t", header=TRUE, col.names = c("Chr", "Pos", "Watson", "Crick"))
}
names(DSBList) <- DSBListNames
Mreads <- MreadR(DSBList)
DSBList <<- DSBList
DSBListNames <<- DSBListNames
Mreads <<- Mreads
nfiles <<- length(DSBList)
}

