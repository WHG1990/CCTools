#' @title CCAnnotate
#' @description loads a RData file attached in th4e data folder of the package and returns it. This is useful for reassigning the name iupon loadin.
#' @param ... a RData or RDS file and returns it. This is useful for reassigning the name iupon loadin
#' @export
CCAnnotate <- function(...){
  options(warn=-1)
  options(warn=0)
  rel_mem <- function(nm) {
    # State the environment
    rm(list=nm, envir = .GlobalEnv )
  }
  data(...)
  delayedAssign("output", ...)
  x=deparse(substitute(...))
  rel_mem(x)
  return(output)
}
#' @title loadRData
#' @description loads a RData file attached in th4e data folder of the package and returns it. This is useful for reassigning the name iupon loadin.
#' @param ... a RData or RDS file and returns it. This is useful for reassigning the name iupon loadin
#' @export
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}
#' @title AppendMe
#' @description combines dataframes making a column stating the source
#' @param dfNames a vector a dataframe names
#' @export
AppendMe <- function(dfNames) {
  do.call(rbind, lapply(dfNames, function(x) {
    cbind(get(x), Library = x)
  }))
}
