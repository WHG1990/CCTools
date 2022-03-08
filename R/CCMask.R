#' @title CCMask
#' @description removes mtDNA,rDNAand 2-micron plasmid from Fullmaps
#' @param path.name A character string defining the RDS file or folder containing FullMaps.
#' @param in.mode Reading in .RDS file or a folder of FullMaps? Choose "rds" or "fullmap".  Defaults to "rds".
#' @param genome "Cer3H4L2" or "pombase2202028" are implemented at the moment. The pombase220208 reference is chromosomes 1—3 plus the mtDNA (4).
#' @param rDNA Mask out the rDNA by setting rDNA = F.
#' @param mtDNA Mask out the mtDNA by setting mtDNA = F.
#' @param twomicron Mask out the twomicron by setting twomicron = F.
#' @examples
#' @author Will Gittens
#' @export
CCMask <- function(path.name = getwd(), in.mode = "fullmap", genome = "Cer3H4L2", rDNA = F, mtDNA = F, twomicron = F, write = T) {
  if(in.mode == "rds") {
    load(path.name)
    exp.name <- substr(basename(path.name), 0, nchar(basename(path.name))-4 )
  }
  if(in.mode == "fullmap") {
    parent.directory <- dirname(path.name)
    CCList(path.name)
  }
DSBListNames <- paste0(DSBListNames, "-")
  if(genome == "Cer3H4L2"){
    if(!rDNA){
      DSBList2 <- lapply(DSBList, function(x){
        x <- x[!(x$Chr == 12 & x$Pos >= 451000 & x$Pos <= 469000), ]
        return(x)
      })
    DSBListNames <- paste0(DSBListNames, "rDNA")
      }
    if(!mtDNA){
      DSBList2 <- lapply(DSBList, function(x){
        x <- x[(x$Chr != 17),]
         return(x)
      })
      DSBListNames <- paste0(DSBListNames, "mtDNA")
    }

    if(!twomicron){
      DSBList2 <- lapply(DSBList, function(x){
        x <- x[(x$Chr != 18),]
        return(x)
      })
      DSBListNames <- paste0(DSBListNames, "2UP")
    }
  }

  if(genome == "pombase220208"){
    if(!rDNA){
      DSBList2 <- lapply(DSBList, function(x){
        x <- x[!(x$Chr == 3 & x$Pos <= 24600 & x$Pos <= 2439550), ]
        return(x)
      })
      DSBListNames <- paste0(DSBListNames, "rDNA")
    }
    if(!mtDNA){
      DSBList2 <- lapply(DSBList, function(x){
        x <- x[(x$Chr != 4),]
        return(x)
      })
      DSBListNames <- paste0(DSBListNames, "mtDNA")
    }

  }

 if(write){CCWrite(DSBList. = DSBList2, filenames = DSBListNames)}
if(!write){
  DSBList <<- DSBList2
  DSBListNames <- DSBListNames
  }
}
