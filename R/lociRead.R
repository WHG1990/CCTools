#' @title lociRead
#' @description Module for reading interpreting loci argument. Used by various CC functions, including CCMap (implemented), CCBin and CCPileup
#' @param in.mode Reading in .RDS file or a folder of FullMaps? Choose "rds" or "fullmap".  Defaults to "rds".
#' @param path.name A character string defining the RDS file or folder containing FullMaps.
#' @param exp.name Manually specify a name for the experiment.
#' @param genome Which genome has the data been aligned to? "Cer3H4L2", "W303" or "hg19". Defaults to "Cer3H4L2".
#' @param bin.width What size bins, in bp? Defaults to 100.
#' @param os Do you want to offset the Crick strand relative to the Watson strand prior to binning? Choose 3 for Top2, or 1 for Spo11. Choose 0 for no ofsetting. Defaults to 3.
#' @param out.mode1 "sparse" or "full" format. Defaults to "full". Defaults to "full"
#' @param out.mode2 Do you want to sum the Watson and Crick strands ("total") or report them seperately ("strand")? Defaults to "strand".
#' @param plot.mode Do you want to plot the binned chromosome maps? (T/F). Defaults to TRUE.
#' @param ylims A single number specifying the maximum limit of the y axis, for both Watson and Crick. OR "auto", to automatically define for the entire set.
#' @param gene.track Do you want to plot the gene track? This is pretty pointless on whole chromosome plots. Awaiting a more sensible implementation. Defaults to FALSE.
#' @examples
#' @author Will Gittens
#' @import data.table
#' @export
lociRead <- function(loci = loci, features = features, chrom. = chrom., pos = pos){
if(is.na(loci[2])) {loci <- c(loci,"mid")}
if(loci[1] == "manual") {
  loci.table <- data.frame(Chr = chrom., Pos = pos)
}
else if(file.exists(loci[1])){
  loci.table <- fread(loci[1])
  if(identical(colnames(loci.table), c("Chr","Start","End","Strand","ExpressionValue","GeneLength","GeneName","HpM","RegionLength","HpMpBP","ntile.length","GeneLengthXExpression"))) {
    loci.table <- data.frame(Chr = loci.table$Chr, Pos = round((loci.table$End-loci.table$Start)/2+loci.table$Start), Name = loci.table$GeneName)
  }
    if(identical(colnames(loci.table), c("sysname", "chr", "type", "start", "stop", "orientation", "genename",
                                 "RNAseq_Brar_00_EXP", "RNAseq_Brar_01_preMeiosis", "RNAseq_Brar_02_DNArep",
                                 "RNAseq_Brar_03_LepZyg", "RNAseq_Brar_04_PachDip", "RNAseq_Brar_05_Met1",
                                 "RNAseq_Brar_06_Ana1", "RNAseq_Brar_07_Met2", "RNAseq_Brar_08_Ana2",
                                 "RNAseq_Brar_09_SporeEarly", "RNAseq_Brar_10_SporeLate"))) {
      loci.table <- data.frame(Chr = loci.table$chr, Pos = round((loci.table$stop-loci.table$start)/2+loci.table$start), Name = loci.table$genename)
  }
  loci.table <- data.frame(loci.table)
}
  else if(loci[1] == "ARS") {
    loci.table <- CCAnnotate(ARS_consensus)
  }
   else if (loci[1] %in% levels(features$type)){

    tempf= subset(features, type==(loci[1]))
    tstart=tempf$start
    tstop=tempf$stop
    tstart[tempf$orientation=="-"] <- tempf$stop[tempf$orientation=="-"]
    tstop[tempf$orientation=="-"] <- tempf$start[tempf$orientation=="-"]

    if (tolower(loci[2])== "m" || loci[2]=="middle"||loci[2]=="mid"){
      loci.table <- data.frame(Chr = tempf$chr, Pos = round((tstop-tstart)/2)+tstart, Name = tempf$genename)
      fplotpos="Mid"
    }
    else if (tolower(loci[2])== "s" ||loci[2]== "start"){
      loci.table <- data.frame(Chr = tempf$chr, Pos = tstart, Name = tempf$genename)
      fplotpos="Start"
    }
    else if (tolower(loci[2])== "stop" ||loci[2]== "end"||loci[2]=="e"){
      loci.table <- data.frame(Chr = tempf$chr, Pos = tstop, Name = tempf$genename)
      fplotpos="End"
    }
   }
    else if (loci[1] %in% features$sysname){

    tempf= subset(features, sysname==loci[1])
    tstart=tempf$start
    tstop=tempf$stop
    tstart[tempf$orientation=="-"] <- tempf$stop[tempf$orientation=="-"]
    tstop[tempf$orientation=="-"] <- tempf$start[tempf$orientation=="-"]

    if (tolower(loci[2])== "m" || loci[2]=="middle"||loci[2]=="mid"){
      loci.table <- data.frame(Chr = tempf$chr, Pos = round((tstop-tstart)/2)+tstart, Name = tempf$sysname)
      fplotpos="Mid"
    }
    else if (tolower(loci[2])== "s" ||loci[2]== "start"){
      loci.table <- data.frame(Chr = tempf$chr, Pos = tstart, Name = tempf$sysname)
      fplotpos="Start"
    }
    else if (tolower(loci[2])== "stop" ||loci[2]== "end"||loci[2]=="e"){
      loci.table <- data.frame(Chr = tempf$chr, Pos = tstop, Name = tempf$sysname)
      fplotpos="End"
    }
  }

  else if (loci[1]%in% features$genename){

    tempf= subset(features, genename== loci[1])
    tstart=tempf$start
    tstop=tempf$stop
    tstart[tempf$orientation=="-"] <- tempf$stop[tempf$orientation=="-"]
    tstop[tempf$orientation=="-"] <- tempf$start[tempf$orientation=="-"]

    if (tolower(loci[2])== "m" || loci[2]=="middle"||loci[2]=="mid"){
      loci.table <- data.frame(Chr = tempf$chr, Pos = round((tstop-tstart)/2)+tstart, Name = tempf$genename)
      fplotpos="Mid"
    }
    else if (tolower(loci[2])== "s" ||loci[2]== "start"){
      loci.table <- data.frame(Chr = tempf$chr, Pos = tstart, Name = tempf$genename)
      fplotpos="Start"
    }
    else if (tolower(loci[2])== "stop" ||loci[2]== "end"||loci[2]=="e"){
      loci.table <- data.frame(Chr = tempf$chr, Pos = tstop, Name = tempf$genename)
      fplotpos="End"
    }
  }

else if(loci[1] == "CTSS") {
  loci.table <- TSS.df
}
else if(loci[1] == "CEN") {
  CEN <- CCAnnotate(Cer3H4L2_centro)
  CEN$Pos <- CEN$chromStart + 0.5*(CEN$chromEnd - CEN$chromStart)
  loci.table <- data.frame(Chr = CEN$Chr, Pos = ceiling(CEN$Pos))

}
  loci.table$Pos <- ceiling(loci.table$Pos)
return(loci.table)
  }


