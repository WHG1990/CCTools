#' @title polygon2
#' @description wrapper for polygon which makes it behave like other types of base R plotting
#' @param x x axis coordinate
#' @param y y axis coordinate
#' @examples
#' @author Will Gittens
#' @export
polygon2 <- function(x,y, ...){
xx <- c(x, rev(x))
yy <- c(rep(0, length(x)), rev(y))
polygon(xx, yy, ...)
}


