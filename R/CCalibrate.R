#' CCalibrate
#' @description Calibrates a CCseq DSBList.
#' @param DSBList A DSBList formatted input.
#' @param Nc A vector describing the number of control cells, or amount of control DNA  added to each sample. The absolute number is not very important.
#' @param Nx A vector describing the number of experimental cells, or amount of experimental DNA in each sample. The absolute number is not very important.
#' @param Mreads_X A vector of Mreads aligning to experimental genome.
#' @param Mreads_C A vector of Mreads aligning to control genome.
#' @author Will Gittens
#' @export
CCalibrate <- function(DSBList,Nc = rep(1,length(DSBList)),Nx = rep(1000, length(DSBList)), Mreads_X = Mreads, Mreads_C){
  # this is the maths based on the Hu paper ("biological chromodynamics...")
  Fx = Mreads_X /(Mreads_X + Mreads_C)
   OR <<- (Nc * Fx)/(Nx * (1- Fx))
  DSBList <- mapply(function(x,y) {
    x$Watson <- y*x$Watson
    x$Crick <- y*x$Crick
    return(x)
  }, DSBList, OR, SIMPLIFY = F)
  return(DSBList)
}
#' mitoCCalibrate
#' @description Special mitochondrial DNA-based calibration a CCseq DSBList.
#' @param DSBList A DSBList formatted input.
#' @param Nc A vector describing the number of control cells, or amount of control DNA  added to each sample. The absolute number is not very important.
#' @param Nx A vector describing the number of experimental cells, or amount of experimental DNA in each sample. The absolute number is not very important.
#' @author Will Gittens
#' @export
mitoCCalibrate <- function(DSBList,Nc = rep(1,length(DSBList)),Nx = rep(1, length(DSBList))){
  experiment <- lapply(DSBList, function(x){
    x <- x[x$Chr %in% c(1:16,18),]
    return(x)
  })

  control <- lapply(DSBList, function(x){
    x <- x[x$Chr == 17,]
    return(x)
  })

  Mreads_X <<- sapply(experiment, function(x){
    sum(x$Watson+x$Crick)/1000000
  }, simplify = T)

  Mreads_C <<- sapply(control, function(x){
    sum(x$Watson+x$Crick)/1000000
  }, simplify = T)

  # this is the maths based on the Hu paper ("biological chromodynamics...")
  Fx = Mreads_X /(Mreads_X + Mreads_C)
  OR <<- (Nc * Fx)/(Nx * (1- Fx))
  DSBList <- mapply(function(x,y) {
    x$Watson <- y*x$Watson
    x$Crick <- y*x$Crick
    return(x)
  }, DSBList, OR, SIMPLIFY = F)
  return(DSBList)
}
#' mitoCCalibrate2
#' @description Special mitochondrial DNA-based calibration a CCseq DSBList.
#' @param DSBList A DSBList formatted input.
#' @param Nc A vector describing the number of control cells, or amount of control DNA  added to each sample. The absolute number is not very important.
#' @param Nx A vector describing the number of experimental cells, or amount of experimental DNA in each sample. The absolute number is not very important.
#' @author George Brown
#' @export
mitoCCalibrate2 <- function(DSBList,Nc = rep(33,length(DSBList)),Nx = rep(4, length(DSBList)),returnOR=F){
   chrom_size=c(230218,813184,320870,1531933,576874,270161,1090940,562643,439888,745751,666816,1078177,924431,784333,1091291,948066,85779,6318)
   sf=chrom_size[17]/sum(chrom_size[-17])

   experiment <- lapply(DSBList, function(x){
    x <- x[x$Chr %in% c(1:16,18),]
    return(x)
  })

  control <- lapply(DSBList, function(x){
    x <- x[x$Chr == 17,]
    return(x)
  })

  Mreads_X <<- sapply(experiment, function(x){
    sum(x$Watson+x$Crick)/1000000
  }, simplify = T)

  Mreads_C <- sapply(control, function(x){
    sum(x$Watson+x$Crick)/1000000
  }, simplify = T)

  # this is the maths based on the Hu paper ("biological chromodynamics...")
  Fx <<- Mreads_X /(Mreads_X + Mreads_C)
  OR <<- (Nc * Fx)/(Nx * (1- Fx))
  DSBList <- mapply(function(x,y) {
    x$Watson <- y*x$Watson*sf
    x$Crick <- y*x$Crick*sf
    return(x)
  }, DSBList, OR, SIMPLIFY = F)
  if (returnOR ==T){return(OR)}else{return(DSBList)}
}
#' rDNACCalibrate
#' @description Special rDNA DNA-based calibration a CCseq DSBList.
#' @param DSBList A DSBList formatted input.
#' @param Nc A vector describing the number of control cells, or amount of control DNA  added to each sample. The absolute number is not very important.
#' @param Nx A vector describing the number of experimental cells, or amount of experimental DNA in each sample. The absolute number is not very important.
#' @author George Brown
#' @export
rDNACCalibrate <- function(DSBList,Nc = rep(75,length(DSBList)),Nx = rep(4, length(DSBList)),returnOR=F){
  chrom_size=c(230218,813184,320870,1531933,576874,270161,1090940,562643,439888,745751,666816,1078177,924431,784333,1091291,948066,85779,6318)
  rs=469000-451000
  sf=rs/(sum(chrom_size)-rs)
  # x<-x[with(x, (Chr == 12 & Pos >= 451000 & Pos <= 469000)),]
  experiment <- lapply(DSBList, function(x){
    x<-x[with(x, (Chr != 12 & Pos >= 451000 & Pos <= 469000)),]
    return(x)
  })
  
  control <- lapply(DSBList, function(x){
    x<-x[with(x, (Chr == 12 & Pos >= 451000 & Pos <= 469000)),]
    return(x)
  })
  
  Mreads_X <<- sapply(experiment, function(x){
    sum(x$Watson+x$Crick)/1000000
  }, simplify = T)
  
  Mreads_C <- sapply(control, function(x){
    sum(x$Watson+x$Crick)/1000000
  }, simplify = T)
  
  # this is the maths based on the Hu paper ("biological chromodynamics...")
  Fx = Mreads_X /(Mreads_X + Mreads_C)
  OR <<- (Nc * Fx)/(Nx * (1- Fx))
  DSBList <- mapply(function(x,y) {
    x$Watson <- y*x$Watson*sf
    x$Crick <- y*x$Crick*sf
    return(x)
  }, DSBList, OR, SIMPLIFY = F)
  if (returnOR ==T){return(OR)}else{return(DSBList)}
  
}
