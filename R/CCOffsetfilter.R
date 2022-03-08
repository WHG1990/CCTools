#' @title CCOffset.filter
#' @description Filters Fullmaps by an offset value
#' @param path.name path to fullmaps
#' @param offset choose the range of offsets to analyse (3 for Top2, 1 for Spo11)
#' @examples CCOffset()
#' @author Will Gittens, George Brown
#' @export
CCOffset.filter<-function(path.name=getwd(),offset=3){
  library("data.table")
  mrg_lst <- function(x, y) {
    out <- list()
    for (n in 1:length(x)) {
      out[[n]] <- merge(x[[n]], y[[n]], nomatch = 0L)
    }
    return(out)
  }

  mrg_lst2 <- function(x, y) {
    out <- list()
    for (n in 1:length(x)) {
      out[[n]] <- merge(x[[n]], y[[n]], all = T)
    }
    return(out)
  }


  ###############################################################################################################################################
  # Import histogram FullMap files for each strain in working directory and tally up the total number of Million mapped reads/ or load .rds
  library(data.table)
  setwd(path.name);parent.directory=getwd()
  CCList(getwd())

  osDSBList <- DSBList
  DSBList2 <- DSBList # copy the data
  osDSBList <-
    lapply(osDSBList, function(x) {
      # drop all rows from osDSBList where Crick is 0
      subset(x, x$Crick > 0)
    })
  DSBList2 <-
    lapply(DSBList2, function(x) {
      # drop all rows from DSBList where Watson is 0
      subset(x, x$Watson > 0)
    })

  DSBList2 <- lapply(DSBList2, data.table)
  DSBList2 <- lapply(DSBList2, setkey, Chr, Pos)

  osDSBList <-  lapply(osDSBList, function(x) { x$Watson <- NULL; return(x)})
  DSBList2 <-   lapply(DSBList2, function(x) { x$Crick <- NULL; return(x)}) # drop Crick column from DSBList2 and Watson column from osDSBList
  time1 <- Sys.time()
  osDSBList <- lapply(osDSBList, data.table)

  # subtract offset from the Pos column of osDSBList, to create osDSBList2.
  osDSBList2 <- lapply(osDSBList, function(x) {
    x$Pos <- x$Pos - offset
    return(x)
  })

  osDSBList2 <- lapply(osDSBList2, setkey, Chr, Pos)

  temp <-
    mrg_lst(DSBList2, osDSBList2)# merge DSBlist2 and osDSBList2, sample by sample, keeping only rows which are present in both

  ##################
  # now do the reverse with merge(..., all = T), to reset the crick positions
  osDSBList <- temp
  DSBList2 <- temp # copy the data

  DSBList2 <- lapply(DSBList2, data.table)
  DSBList2 <- lapply(DSBList2, setkey, Chr, Pos)
  osDSBList <-  lapply(osDSBList, function(x) {x$Watson <- NULL; return(x)})
  DSBList2 <-   lapply(DSBList2, function(x) {x$Crick <- NULL; return(x)}) # drop Crick column from DSBList2 and Watson column from osDSBList
  osDSBList <- lapply(osDSBList, data.table)
  osDSBList2 <- lapply(osDSBList, function(x) {
    x$Pos <- x$Pos - (offset-2*offset)
    return(x)
  })

  osDSBList2 <- lapply(osDSBList2, setkey, Chr, Pos)
  temp <-
    mrg_lst2(DSBList2, osDSBList2)# merge DSBlist2 and osDSBList2, sample by sample, keeping only rows which are present in both
  temp <- lapply(temp, function(x) {
    x[is.na(x)] <- 0
    return(x)
  })

  assign(paste0("DSBList_", offset, "_bp_offset"), temp)

  ifelse(!dir.exists(file.path(parent.directory, paste0(offset, "_bp_Offset"))), dir.create(file.path(parent.directory, paste0(offset, "_bp_Offset"))), FALSE); setwd(file.path(parent.directory, paste0(offset, "_bp_Offset"))); directory2 <- getwd()
  foldern="offsetFiltered"
  # save(paste0("DSBList_", offset, "_bp_offset"), DSBListNames, exp.name, nfiles, Mreads, file = paste0(exp.name))
  ifelse(!dir.exists(file.path(directory2, paste0(foldern))), dir.create(file.path(directory2, paste0(foldern))), FALSE); setwd(file.path(directory2, paste0(foldern)))
  for (i in 1:length(temp)){
    write.table(temp[[i]], file = paste0("FullMap.",DSBListNames[i],"_3bpOffset.txt"), sep = "\t", row.names = F, quote = F)
  }


  #generating reciprocal subset
  DSBList.1 <- NULL
  for (k in 1:length(DSBList)){
    temp <- data.table(DSBList[[k]])
    temp2 <- data.table(DSBList_3_bp_offset[[k]])
    setkey(temp, Chr, Pos)
    setkey(temp2, Chr, Pos)
    temp<- merge(temp, temp2, all.x = T)
    temp<-temp[is.na(temp$Watson.y)&is.na(temp$Crick.y)]
    DSBList.1[[k]] <- temp
  }

  DSBList.1 <- lapply(DSBList.1, function(x){
    x$Watson.y <- NULL
    x$Crick.y <- NULL
    colnames(x) <- c("Chr", "Pos", "Watson", "Crick")
    return(x)
  })
  foldern="NOToffsetFiltered"
  ifelse(!dir.exists(file.path(directory2, paste0(foldern))), dir.create(file.path(directory2, paste0(foldern))), FALSE); setwd(file.path(directory2, paste0(foldern)))
  for (i in 1:length(DSBList.1)){
    write.table(DSBList.1[[i]], file = paste0("FullMap.",DSBListNames[i],"_NOT3bpOffset.txt"), sep = "\t", row.names = F, quote = F)
  }
  setwd(path.name)
}
