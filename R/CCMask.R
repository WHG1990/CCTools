#' @title CCMask
#' @description removes mtDNA,rDNAand 2-micron plasmid from Fullmaps
#' @param path.name A character string defining the RDS file or folder containing FullMaps.
#' @param in.mode Reading in .RDS file or a folder of FullMaps? Choose "rds" or "fullmap".  Defaults to "rds".
#' @examples
#' @author Will Gittens
#' @export
CCMask <- function(path.name = getwd(), in.mode = "fullmap") {
  if(in.mode == "rds") {
    load(path.name)
    exp.name <- substr(basename(path.name), 0, nchar(basename(path.name))-4 )
  }
  if(in.mode == "fullmap") {
    parent.directory <- dirname(path.name)
    CCList(path.name)
  }

  DSBList2 <- lapply(DSBList, function(x){
    x <- x[!(x$Chr == 12 & x$Pos >= 451000 & x$Pos <= 469000), ]
    x <- x[(x$Chr != 17),]
    x <- x[(x$Chr != 18),]
    return(x)
  })
  CCWrite(DSBList. = DSBList2, filenames = paste0(DSBListNames, "_-mtDNArDNA2UP"))
}
