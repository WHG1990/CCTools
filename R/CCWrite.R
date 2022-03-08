#' @title CCWrite
#' @description Write a DSBList into separate FullMap fies.
#' @param DSBList. A DSBList object.
#' @param filenames The names of the files. Defaults to DSBListNames.
#' @examples
#' @author Will Gittens
#' @export
CCWrite <- function(DSBList. = DSBList, filenames = DSBListNames){
  for (i in 1:length(DSBList.)){
    write.table(DSBList.[[i]], file = paste0("FullMap.", filenames[i],".txt"), sep = "\t", row.names = F, quote = F)
  }
}
