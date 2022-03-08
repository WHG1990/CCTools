#' @title CCSplit
#' @description Splits fullmaps into seprate chromosome-seperated maps
#' @param path.name A character string defining the RDS file or folder containing FullMaps.
#' @param in.mode Reading in .RDS file or a folder of FullMaps? Choose "rds" or "fullmap".  Defaults to "rds".
#' @examples
#' @author Will Gittens
#' @export
CCSplit <- function(path.name = getwd(), in.mode = "fullmap") {
  if(in.mode == "rds") {
    load(path.name)
    exp.name <- substr(basename(path.name), 0, nchar(basename(path.name))-4 )
  }
  if(in.mode == "fullmap") {
    parent.directory <- dirname(path.name)
    CCList(path.name)
  }
  
  DSBList2 <- lapply(DSBList, function(x){
x <- split(x, x$Chr)
    return(x)
  })
  DSBList2 <- unlist(DSBList2, recursive = F)
  DSBListNames <- rep(DSBListNames,each = 18)
  suff <- sapply(DSBList2, function(x){
    paste0("_Chr",x$Chr[1])
  })
  DSBListNames <- paste0(DSBListNames, suff)
    CCWrite(DSBList. = DSBList2, filenames = paste0(DSBListNames))
}



