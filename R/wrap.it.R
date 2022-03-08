#' wrap.it
#' @description wraps a character string (x) every n characters
#' @param x A character string (e.g an axis label)
#' @param n The number of characters on each line.
#' @author Will Gittens
#' @export
wrap.it <- function(x, n)
{
  sapply(x, function(y) paste(

    paste(sapply(seq(1, nchar(y), n), function(i) paste0(substring(y, i, min(i + (n-1), nchar(y))), '\n')), collapse=''),
                              collapse = "\n"),
         USE.NAMES = F)
}


