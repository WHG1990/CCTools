#' @title MreadR
#' @description Calculates the total number of reads in a DSBList, outputted as a vector of length equal to the length of the DSBList.
#' @param x A DSBList object, or a path to a DSBList saved in RDS format, expressed as a character string.
#' @examples
#' @author Will Gittens 2020
#' @export
MreadR <- function(x) {
  if(class(x) == "character") {
    load(x)

    sapply(DSBList, function(x){

      if("Total" %in% colnames(x)){
        sum(x$Total)/1000000
      } else {
      sum(x$Watson + x$Crick)/1000000
      }
    })
  }
  else {
    sapply(x, function(x){
      if("Total" %in% colnames(x)){
        sum(x$Total)/1000000
      } else {
        sum(x$Watson + x$Crick)/1000000
      }
    })
  }
 }
