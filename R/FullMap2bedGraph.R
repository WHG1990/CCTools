#' @title FullMap2BedGraph
#' @description Converts FullMaps downloaded from GEO (SPECIFICALLY), into bedGraphs which can be uploaded as custom tracks to UCSC.
#' @param path.name Character string defining the pathname to the fullmap file. Files must be in the format SAMPLENAME.FullMap.txt
#' @param scaling.factor a scaling factor used to boost the values. This is useful in calibration, for example. Defaults to 1.
#' @param out.dir A character string defining the pathname to the output directory. Defaults to current working directory.
#' @param chrom A vector defining the chromosomes to process. Defaults to all chromosomes, including chromosome 3, 17 and 18, which may cause compatibility issues with UCSC.
#' @param mode Can be "break" or "base", depending on whether you want to output the position of the DNA break (between bases), or the position of the base which is covalently linked to the protein.
#' @examples
#' @author Will Gittens
#' @export
FullMap2bedGraph <- function(path.name, scaling.factor = 1, out.dir = getwd(), chrom = 1:18, mode = "base") {
  library(data.table)
  setwd(out.dir)
  if(mode == "base"){
    data <- fread(path.name)
    data$Watson <- data$Watson*scaling.factor
    data$Crick <- data$Crick*scaling.factor
    file.name <- basename(path.name)
    file.name <- sub("\\..*", "", file.name)
    data <- data[data$Chr %in% chrom,]
    data2 <- data[data$Watson>0,]
    track_definition <- paste0("track type=bedGraph name=", file.name, "_WatsonBase visibility=full color=205,0,0")
    a =  data.frame(paste0("chr",as.roman(data2$Chr)), data2$Pos-1, data2$Pos, data2$Watson)
    cat(paste0(track_definition, "\n"),file=paste0(file.name,"_WatsonBase.bedGraph"))
    write.table(a, paste0(file.name,"_WatsonBase.bedGraph"),sep="\t",append=TRUE, col.names=F, row.names = F, quote = F)
    data2 <- data[data$Crick>0,]
    track_definition <- paste0("track type=bedGraph name=", file.name, "_CrickBase visibility=full color=28,134,238")
    b =  data.frame(paste0("chr",as.roman(data2$Chr)), data2$Pos-1, data2$Pos, data2$Crick)
    cat(paste0(track_definition, "\n"),file=paste0(file.name,"_CrickBase.bedGraph"))
    write.table(b, paste0(file.name,"_CrickBase.bedGraph"),sep="\t",append=TRUE, col.names=F, row.names = F, quote = F)
  }
  if(mode == "break") {
    data <- fread(path.name)
    data$Watson <- data$Watson*scaling.factor
    data$Crick <- data$Crick*scaling.factor
    file.name <- basename(path.name)
    file.name <- sub("\\..*", "", file.name)
    data <- data[data$Chr %in% chrom,]
    data2 <- data[data$Watson>0,]
    track_definition <- paste0("track type=bedGraph name=", file.name, "_WatsonBreak visibility=full color=205,0,0")
    a =  data.frame(paste0("chr",as.roman(data2$Chr)), data2$Pos-1, data2$Pos-1, data2$Watson)
    cat(paste0(track_definition, "\n"),file=paste0(file.name,"_WatsonBreak.bedGraph"))
    write.table(a, paste0(file.name,"_WatsonBreak.bedGraph"),sep="\t",append=TRUE, col.names=F, row.names = F, quote = F)
    data2 <- data[data$Crick>0,]
    track_definition <- paste0("track type=bedGraph name=", file.name, "_CrickBreak visibility=full color=28,134,238")
    b =  data.frame(paste0("chr",as.roman(data2$Chr)), data2$Pos, data2$Pos, data2$Crick)
    cat(paste0(track_definition, "\n"),file=paste0(file.name,"_CrickBreak.bedGraph"))
    write.table(b, paste0(file.name,"_CrickBreak.bedGraph"),sep="\t",append=TRUE, col.names=F, row.names = F, quote = F)
  }
}
