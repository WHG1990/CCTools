#' @title rounddown
#' @description Rounds down a number to the specified base.
#' @param x Input number.
#' @param base The base number
#' @author Will Gittens
#' @export
rounddown <- function(x,base){
  base*floor(x/base)
}
