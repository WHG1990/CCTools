#' @title CCSearch
#' @description Looks locally around path.name for a file matching file.name, and returns it's path.
#' @param path.name Path to the folder around which you wish to search for the file.
#' @param file.name The file name to search for.
#' @examples
#' @author Will Gittens
#' @export
CCSearch <- function(path.name, file.name){
  if(missing(path.name)){path.name <- getwd()}
  path.vec <- strsplit(path.name, "/")[[1]]
  path.vec <- path.vec[2:length(path.vec)]
  path.vec2 <- path.vec
  for (i in 1:length(path.vec)){
    temp <- c(path.vec[1:i])
    temp <- paste0(temp, sep = "/", collapse = "")
    temp <- paste0("/", temp)
    path.vec2[i] <- temp
  }
  path.vec <- path.vec2
  path.vec <- path.vec[(ceiling(length(path.vec)/2)):length(path.vec)]
  path.vec <- rev(path.vec)
  for (i in 1:length(path.vec)){
    str <- paste0('list.files(path = path.vec[i], pattern = "', file.name, '", recursive = T)')
    out <- eval(parse(text = str))
    if (length(out) > 0) { break }
  }
  if(length(out) > 0){out <- paste0(path.vec[i],out)}
  return(out)
}
