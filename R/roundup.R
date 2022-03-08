#' @title roundup
#' @description Rounds up a number to the specified base.
#' @param x Input number.
#' @param base The base number
#' @author Will Gittens
#' @export
roundup <- function(x,base){
  base*ceiling(x/base)
}
