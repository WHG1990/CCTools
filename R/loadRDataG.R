#' @title loadRDataG
#' @description loads a RData or RDS file and returns it. This is useful for reassigning the name iupon loadin.
#' @param loads a RData or RDS file and returns it. This is useful for reassigning the name iupon loadin.
#' hey
#' @export
loadRDataG <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}
