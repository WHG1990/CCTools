#' @title hannSlide
#' @description Smooth an input vector (or data table column), with a sliding Hann window of defined width.
#' @param vec Input vector or (or data table column) to smooth.
#' @param win Width of hanning window in bp.
#' @examples
#' @author Will Gittens
#' @import e1071
#' @export
hannSlide <- function(vec, win = 50) {
require("e1071")
hw=hanning.window(win) #create hanning window (require package e1071 to be loaded)
scalar=1/sum(hw) #scalar [1]
out <- as.numeric(scalar*(stats::filter(vec,hw)))
out[is.na(out)] <- 0
return(out)
}
