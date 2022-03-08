#' @title CCMap
#' @description Maps CCs in a specific chromosome region, alongside various annotations.
#' @param in.mode Reading in .RDS file or a folder of FullMaps? Choose "rds" or "fullmap".  Defaults to "rds".
#' @param path.name A character string defining the RDS file or folder containing FullMaps.
#' @param exp.name Manually specify a name for the experiment.
#' @param genome Which genome has the data been aligned to? "Cer3H4L2", "W303" or "hg19". Defaults to "Cer3H4L2".
#' @param phase "meiosis" or "vegetative"
#' @param chrom. Which chromosome to plot?
#' @param pos Which coordinate to centre plot on, in bp?
#' @param window.w How wide a region to plot in bp?
#' @param os Do you want to offset the Crick strand relative to the Watson strand prior to binning? Choose 3 for Top2, or 1 for Spo11. Choose 0 for no ofsetting. Defaults to 3.
#' @param samples Which samples to plot, as a vector of indexes in the DSBList.
#' @param out.mode "pdf" or "png" file format for output.
#' @param ylims A single number specifying the maximum limit of the y axis, for both Watson and Crick. OR "auto", to automatically define for each sample.
#' @param gene.track Do you want to plot the gene track? Defaults to TRUE.
#' @param region.mode Do you want to manually specify a region, or do you want to multiple maps centred on every TSS. CHoose "manual" (Default) or "TSS".
#' @param smooth.mode This option controls whether to apply VariableX ("VX") or Hann window ("HW") smoothing to the broad overlay. Deafults to no smoothing.
#' @param Fsize Fsize for VariableX (See the documentation for VX function). Only used if smooth.mode == "VX".
#' @param win Hanning window width. Only used if smooth.mode == "HW".
#' @param norm.factor This can be used to scale the height of the smoothed track (works for Hann  or VariableX smoothing). Defaults to 3.
#' @param plot.dimensions Specify the plot size. Can be any size from the "A sizing" scale (e.g "A1"). However this is not implemented properly yet. Defaults to "wide", which is fine.
#' @examples
#' @author Will Gittens
#' @import data.table
#' @export
CCTally <- function(in.mode = "rds", path.name = getwd(), exp.name, genome = "Cer3H4L2", phase = "meiosis", samples, os = 3, mcal = F, cal.path, chrom., pos, up.w, down.w, ylims = 50, region.mode = "manual", smooth.mode = "none", Fsize = 101, win = 50, norm.factor = 3, gene.track = T, MNase.track = T, out.mode = "png", plotdimensions = "wide", loci) {
{
  library("data.table")
  library("e1071")
  library("shape")
  # if((window.w %% 2)>0){window.w = window.w+1}

  if(genome == "W303") {print("WARNING: W303 Gene/TSS Annotations are not available. Cer3H4L2 annotations are used instead (for now).")}
  if(ylims == "auto") {y.mode <- "auto"} else {y.mode <- "fixed"}
  if(y.mode == "fixed") {fixed.y <- ylims}
    ################# data loading
    if(in.mode == "rds") {
      load(path.name)
      exp.name <- substr(basename(path.name), 0, nchar(basename(path.name))-4 )
      parent.directory <- dirname(path.name)
    }
    if(in.mode == "fullmap") {
      parent.directory <- dirname(path.name)
      CCList(path.name)
      exp.name <- exp.name
    }
  if(exists("DSBList.bf")){DSBList <- DSBList.bf}
  if(exists("DSBList.pooled")){DSBList <- DSBList.pooled}
  if(exists("pooled.DSBList")){DSBList <- pooled.DSBList}
DSBList <- lapply(DSBList, data.table)
DSBList <- lapply(DSBList, setkey, Chr, Pos)

if(missing(samples)){samples <- 1:length(DSBList)}

################# TSS loading
if(genome == "Cer3H4L2") {
  TSS.df <- CCAnnotate(Cer3_TSS_positions_forwhich_GSE36958_vegetative_available)
  features <- CCAnnotate(Cer3H4L2_AllElementsDUB)
  chrom.lengths <- CCAnnotate(Cer3H4L2_chromlengths)
  if(phase == "meiosis"){MNase <- CCAnnotate(MNase_meiosis_GSM1424408)}
  if(phase == "vegetative"){MNase <- CCAnnotate(MNase_Veg_GSM1849297_GSM1849298_R1R2)}
}

setkey(TSS.df, Chr, Pos)
################# gene loading=
if(genome == "W303") {
  print("Warning! This is a Cer3 Gene track!")
  TSS.df <- CCAnnotate(Cer3_TSS_positions_forwhich_GSE36958_vegetative_available)
  features <- CCAnnotate(Cer3H4L2_AllElementsDUB)
  chrom.lengths <- CCAnnotate(W303_chromlengths)
  if(phase == "vegetative"){MNase <- CCAnnotate(MNase_Veg_GSM1849297_GSM1849298_R1R2); print("Warning! This is a Cer3 MNase track!")}
}


  #########################################################################################################
  ############ START HERE ONCE DATAFRAMES ARE LOADED ######################
  #########################################################################################################

  if(region.mode == "manual") {
    loci.table <- data.frame(Chr = chrom., Pos = pos)
  }

  if(region.mode == "CTSS") {
    loci.table <- TSS.df
    motifsW <- loci.table[loci.table$Strand == "+",]
    motifsC <- loci.table[loci.table$Strand == "-",]
    motifsW$Start <- motifsW$Pos
    motifsW$End <- motifsW$Pos2
    motifsC$Start <- motifsC$Pos2
    motifsC$End <- motifsC$Pos
    loci.table <- rbind.data.frame(motifsW,motifsC)
    loci.table = data.frame(loci.table)

  }

  if(region.mode == "multi") {
    loci.table <- fread(loci)
    loci.table <<- data.frame(loci.table)
  }

  setwd(parent.directory)
  ifelse(!dir.exists(file.path(parent.directory, "Tally")), dir.create(file.path(parent.directory, "Tally")), FALSE); setwd(file.path(parent.directory, "Tally")); directory2 <- getwd()
  dir.create(file.path(directory2, paste0("Upwindow_", up.w,"Downwindow_", down.w,  "_fsize-", Fsize)), showWarnings = FALSE)
  setwd(file.path(directory2, paste0("Upwindow_", up.w,"Downwindow_", down.w,  "_fsize-", Fsize)))#Plot range minimum (bp); # Plot range width (bp); #Plot range maximum (bp)

loci.columns <- ncol(loci.table)
out.loci <- data.frame(matrix(data = 0, nrow = nrow(loci.table), ncol = length(samples)))
colnames(out.loci) <- DSBListNames[samples]

  ###############################################################################################################################################
  for (n in 1:nrow(loci.table)){
    if(region.mode == "CTSS") {

      chrom = loci.table[n,"Chr"];
      if(loci.table[n,"Strand"] == "+"){
      xl1 <- loci.table[n,"Start"]-up.w;
      xl2 <-loci.table[n,"End"]+down.w
      }

      if(loci.table[n,"Strand"] == "-"){
        xl1 <- loci.table[n,"Start"]-down.w;
        xl2 <-loci.table[n,"End"]+up.w
      }

      } #Plot range minimum (bp); # Plot range width (bp); #Plot range maximum (bp)

counter = 0
    ###############################################################################################################################################
    for (k in samples){
      counter = counter+1

      sae2.0=DSBList[[k]][J(chrom, xl1:xl2)]
      sae2.0$Total <- sae2.0$Watson+sae2.0$Crick
      #Decompression code here
      sae2.1 <- data.frame(Chr=chrom, Pos=(xl1:xl2)) # Creates expanded empty dataframe with Chr and Pos locations
      sae2.1 <- merge(sae2.1,sae2.0, all=TRUE) # Merge expanded empty dataframe with compressed sae2.1 dataframe
      if (os != 0){sae2.1$Crick <- c(sae2.1$Crick[-(seq(os))], rep(NA, os))}  # OFFSET CRICK IF YOU WANT
      sae2.1[is.na(sae2.1)] <- 0 # Convert all NA values to zero
      sae2.1$Total=sae2.1$Watson+sae2.1$Crick

      out.loci[n,counter]  <- sum(sae2.1$Watson+sae2.1$Crick)/Mreads[k]

    }
   }

loci.table <<- cbind.data.frame(loci.table, out.loci)
}
}
