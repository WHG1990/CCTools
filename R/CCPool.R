#' @title CCPool
#' @description Pools multiple Fullmaps together, by summing total reads.
#' @param path.name A character string defining folder containing individual FullMaps.
#' @param exp.name Manually specify a name for the outputted pool file. Defaults to a concatenation of individual file names.
#' @param out.mode output in "txt", "rds", or "both" mode.
#' @examples
#' @author Will Gittens
#' @import data.table
#' @import stringr
#' @export
CCPool <- function(path.name, exp.name, out.mode = "txt") {
  library("stringr")
  library("data.table")
setwd(path.name)
  DSBList.2<- list()
  file.list <- list.files(pattern = ".rds|.txt")
  abbrev.list <- str_match(file.list, "FullMap.(.*?).txt")[,2]
  temp <- str_match(file.list, "(.*?).rds")[,2]
  abbrev.list[is.na(abbrev.list)] <- temp[!is.na(temp)]
  ext.list <- substr(file.list, nchar(file.list)-2, nchar(file.list))
  ################# data loading
for (g in 1:length(file.list)) {
  if(ext.list[g] == "txt"){
    DSBList <- fread(file.list[g], sep = "\t", header=TRUE, col.names = c("Chr", "Pos", "Watson", "Crick"))
    DSBList.2[[g]] <- DSBList
    rm(DSBList)
  }
    if(ext.list[g] == "rds"){
    load(file.list[g])
      if(exists("DSBList.bf")){DSBList <- DSBList.bf; rm(DSBList.bf)}
      if(exists("DSBList.pooled")){DSBList <- DSBList.pooled; rm(DSBList.pooled)}
      if(exists("pooled.DSBList")){DSBList <- pooled.DSBList; rm(pooled.DSBList)}
      if(exists("dup.DSBList")){DSBList <- dup.DSBList; rm(dup.DSBList)}
      if(exists("undup.DSBList")){DSBList <- undup.DSBList; rm(undup.DSBList)}
      if("list" %in% class(DSBList)) {stop("ERROR: one or more RDS files contains more than one Fullmap. CCPool only works on inidividual Fullmaps.")}
    DSBList.2[[g]] <- DSBList
    rm(DSBList)
  }
}
DSBList <- DSBList.2; rm(DSBList.2)
DSBList <- lapply(DSBList, data.table)
DSBList <- lapply(DSBList, setkey, Chr, Pos)
DSBList.pooled <- DSBList[1]
for (i in 2:length(DSBList)) {
  DSBList.pooled[[1]] <- merge(DSBList.pooled[[1]], DSBList[[i]], by = c("Chr", "Pos"), all = T)
  DSBList.pooled <- lapply(DSBList.pooled, function(x){
    x[is.na(x)] <- 0L
    return(x)
  })
  DSBList.pooled <- lapply(DSBList.pooled, function(x){
    x$Watson <- x$Watson.x + x$Watson.y
    x$Crick <- x$Crick.x  + x$Crick.y
    return(x)
  })
  DSBList.pooled <- lapply(DSBList.pooled, function(x){
    x$Watson.x <- NULL;  x$Watson.y <- NULL;  x$Crick.x <- NULL;  x$Crick.y <- NULL
    return(x)
  })
}
DSBListNames <- paste(abbrev.list, collapse = "_")
DSBListNames <- paste0(DSBListNames, "_Pool")
DSBList <- DSBList.pooled
nfiles <- length(DSBList)
Mreads <- sapply(DSBList, function(x){
  sum(x$Watson + x$Crick)/1000000
})
if(!missing(exp.name)){DSBListNames <- exp.name}
ifelse(!dir.exists(file.path(path.name, "Pool")), dir.create(file.path(path.name, "Pool")), FALSE); setwd(file.path(path.name, "Pool"))
if(out.mode == "rds"){save(DSBList, DSBListNames, nfiles, MMreads, Mreads, file = paste0(DSBListNames[1], ".rds"))}
if(out.mode == "fullmap" | out.mode == "txt"){
  CCWrite(DSBList. = DSBList, filenames = DSBListNames)
}
if(out.mode == "both"){
  save(DSBList, DSBListNames, nfiles, MMreads, Mreads, file = paste0(DSBListNames[1], ".rds"))
  CCWrite(DSBList. = DSBList, filenames = DSBListNames)
}
}


