#' @title CCMap
#' @description Maps CCs in a specific chromosome region, alongside various annotations.
#' @param in.mode Reading in .RDS file or a folder of FullMaps? Choose "rds" or "fullmap".  Defaults to "rds".
#' @param path.name A character string defining the RDS file or folder containing FullMaps.
#' @param exp.name Manually specify a name for the experiment.
#' @param genome Which genome has the data been aligned to? "Cer3H4L2", "W303" or "hg19". Defaults to "Cer3H4L2".
#' @param phase "meiosis" or "vegetative"
#' @param chrom. Which chromosome to plot?
#' @param pos Which coordinate to centre plot on, in bp?
#' @param featureplot (default=0) changes plot cords to match feature for example c("TOP2","start") plots the at the start position of TOP2
#' @param window.w How wide a region to plot in bp?
#' @param os Do you want to offset the Crick strand relative to the Watson strand prior to binning? Choose 3 for Top2, or 1 for Spo11. Choose 0 for no ofsetting. Defaults to 3.
#' @param samples Which samples to plot, as a vector of indexes in the DSBList.
#' @param out.mode "pdf" or "png" file format for output.
#' @param ylims A single number specifying the maximum limit of the y axis, for both Watson and Crick. OR "auto", to automatically define for each sample.
#' @param gene.track Do you want to plot the gene track? Defaults to TRUE.
#' @param loci Do you want to manually specify a region ("manual"), a specific single gene ("GENENAME"), or multiple maps centred on a types of features in the AllElementsDUB table (e.g "gene", "CEN", "ARS_consensus_sequence". For more options run CCFeatures() ). You can also provide a vector where the second element is "start", "mid", or "end", to centre the pileup on , for example, the end of a gene.
#' @param smooth.mode This option controls whether to apply VariableX ("VX") or Hann window ("HW") smoothing to the broad overlay. Defults to no smoothing.
#' @param Fsize Fsize for VariableX (See the documentation for VX function). Only used if smooth.mode == "VX".
#' @param win Hanning window width. Only used if smooth.mode == "HW".
#' @param norm.factor This can be used to scale the height of the smoothed track (works for Hann  or VariableX smoothing). Defaults to 3.
#' @param plot.dimensions Specify the plot size. Can be any size from the "A sizing" scale (e.g "A1"). However this is not implemented properly yet. Defaults to "wide", which is fine.
#' @param origin.line Do you want to plot a dotted line in the centre? Useful for demonstrating offset at finescale resolution.
#' @param plot.scale Do you want to keep the overall plot dimensions constant? If so set to "fixed". Defaults to unfixed, which means that the overal height of the plot increases with more samples. This was George's solution to the overlapping label issue.
#' @param loci.arrow Do you want to plot an arrow highlighting the centre of the plot? E.g the locus
#' @param gene.names Do you want the gene track to be labelled with "genename" (e.g TEL1), or "sysname" (e.g YBL088C)?
#' @examples
#' @author Will Gittens,George Brown
#' @import data.table
#' @export
CCMap <- function(in.mode = "fullmap", path.name = getwd(), exp.name, genome = "Cer3H4L2", phase = "meiosis", samples, os = 0, mcal = F, cal.path, chrom., pos, window.w, ylims = 50, loci = "manual", smooth.mode = "none", Fsize = 101, win = 50, norm.factor = 3, gene.track = T, ARS.track = T, Ty.track = T, MNase.track = T, scc1.track = T, disparity = F, out.mode = "png", plotdimensions = "wide", featureplot=0, names, origin.line = F, portrait = F, basic.mode = F, plot.scale = F, resolution=400, loci.arrow = F, gene.names = "genename") {

  library("data.table")
  library("e1071")
  library("shape")
  if((window.w %% 2)>0){window.w = window.w+1}

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
  if(samples[1] == "all"){samples <- 1:length(DSBList)}
  if(!missing(names)){DSBListNames <- names}
  ################# TSS loading
  if(genome == "Cer3H4L2") {
    TSS.df <- CCAnnotate(Cer3H4L2_CTSS_positions_forwhich_GSE36958_vegetative_available)
    features <- CCAnnotate(Cer3H4L2_AllElementsDUB)
    chrom.lengths <- CCAnnotate(Cer3H4L2_chromlengths)
    setkey(TSS.df, Chr, Pos)
  }

  if(genome == "pombase220208"){
  features <- CCAnnotate(pombase220208_AllElements)
  chrom.lengths <- CCAnnotate(pombase220208_chromlengths)
  }

  ################# MNase loading
  if(genome == "Cer3H4L2") {
    if(phase == "meiosis"){MNase <- CCAnnotate(MNase_meiosis_GSM1424408); Nucleosomes = "Meiotic\nMNAse"}
    if(phase == "vegetative"){MNase <- CCAnnotate(MNase_Veg_GSM1849297_GSM1849298_R1R2); Nucleosomes = "Vegetative\nMNAse"}
  }

  ################# Scc1 loading
  if(genome == "Cer3H4L2") {
    if(phase == "vegetative"){scc1 <- CCAnnotate(GSM1712311_Fig4_G1_releasing_60min_IP);Cohesin="Scc1"}
    if(phase == "meiosis"){scc1 <- CCAnnotate(GSM739669_Rec8_4h_Cer3H4L2);Cohesin="Rec8"}
  }

  ################# gene loading=
  if(genome == "W303") {
    print("Warning! This is a Cer3 Gene track!")
    TSS.df <- CCAnnotate(Cer3H4L2_CTSS_positions_forwhich_GSE36958_vegetative_available)
    features <- CCAnnotate(Cer3H4L2_AllElementsDUB)
    chrom.lengths <- CCAnnotate(W303_chromlengths)
    if(phase == "vegetative"){MNase <- CCAnnotate(MNase_Veg_GSM1849297_GSM1849298_R1R2); print("Warning! This is a Cer3 MNase track!")}
    if(phase == "vegetative"){scc1 <- CCAnnotate(GSM1712311_Fig4_G1_releasing_60min_IP); print("Warning! This is a Cer3 Scc1 track!")}
  }


  # ################# CTCF loading (human only) NOT IMPLEMENTED
  # if(genome == "human"){
  #   setwd("/Users/wg45/Dropbox/Work/Top2Seq/WD/")
  #   load("Properly_averaged_RPE1_CTCF_Motifs_signalval.rds")
  #   setkey(motifs, Chr, Pos)
  # }

  # c("sysname", "genename", "description", "type", "start", "stop", "Chr", "Strand")
  if ("Strand" %in% colnames(features)) {colnames(features)[colnames(features) == "Strand"] <- "orientation"}
  if ("Chr" %in% colnames(features)) {colnames(features)[colnames(features) == "Chr"] <- "chr"}

    #########################################################################################################
    ############ START HERE ONCE DATAFRAMES ARE LOADED ######################
    #########################################################################################################
    # if(loci == "manual" && length(pos) == 2) {
    #   xl1 <- pos[[1]]
    #   xl2 <- pos[[2]]
    #   window.w = xl2-xl1
    # }

    if(is.na(loci[2])) {loci <- c(loci,"mid")}
    loci.table <- lociRead(loci = loci, features = features, chrom. = chrom., pos = pos)
    if (tolower(loci[2])== "m" || loci[2]=="middle"||loci[2]=="mid"){
      fplotpos="Mid"
    }
    else if (tolower(loci[2])== "s" ||loci[2]== "start"){
      fplotpos="Start"
    }
    else if (tolower(loci[2])== "stop" ||loci[2]== "end"||loci[2]=="e"){
      fplotpos="End"
    }

    vec <- c(66.2,46.8,33.1,23.4,16.5,11.7,8.3,5.8,4.1,2.9,2.0,1.5,1, 7.5, 4.6, 3.3)
    vec2 <- c(93.6,66.2,46.8,33.1,23.4,16.5,11.7,8.3,5.8,4.1,2.9,2.0,1.5, 13.3, 8.3, 5.8)
    vec3 <- c("4A0", "2A0", "A0", "A1", 'A2', 'A3', 'A4', 'A5', 'A6', "A7", 'A8', 'A9', 'A10', "wide", "wide2", "wide3")
    papersize <- data.frame(vec3,vec2,vec)
    if(portrait){papersize <- data.frame(vec3,vec,vec2)}
    plotdimension <- as.numeric(papersize[match(plotdimensions, vec3),2:3])

    setwd(parent.directory)
    ifelse(!dir.exists(file.path(parent.directory, "Mapping")), dir.create(file.path(parent.directory, "Mapping")), FALSE); setwd(file.path(parent.directory, "Mapping")); directory2 <- getwd()
    ifelse(!dir.exists(file.path(directory2, exp.name)), dir.create(file.path(directory2, exp.name)), FALSE); setwd(file.path(directory2, exp.name)); directory3 <- getwd()
    dir.create(file.path(directory3, paste0("window_", window.w, "_fsize-", Fsize)), showWarnings = FALSE)
    setwd(file.path(directory3, paste0("window_", window.w, "_fsize-", Fsize)))#Plot range minimum (bp); # Plot range width (bp); #Plot range maximum (bp)

    # xl1<-loci.table[1,"Pos"]-1/2*window.w; xl2<-xl1+window.w
    #
    # if(ylims == "auto.set") {
    #   sae2.0.list <- lapply(DSBList, function(x) {
    #     x <- x[J(chrom, xl1:xl2), nomatch=0L]
    #     return(x)
    #   })
    #
    #   ylims2 = c(-max(max(x$Watson/Mreads[k]), max(x$Crick/Mreads[k])),max(max(x$Watson/Mreads[k]), max(x$Crick/Mreads[k])))}
    #
    # }

    ###############################################################################################################################################
    for (n in 1:nrow(loci.table)){
      chrom = loci.table[n,"Chr"]; xl1 <- loci.table[n,"Pos"]-1/2*window.w; xl2 <- xl1+window.w  #Plot range minimum (bp); # Plot range width (bp); #Plot range maximum (bp)

      if(loci[1] == "manual") {plot.name <- paste0(exp.name,"_Chr", chrom., "_Pos", xl1, "-", xl2, "_OS", os)} else {plot.name <- paste0(exp.name,"_",n, "_", loci.table$Name[n],"_", fplotpos, "_", window.w, "bp_OS", os)}
      if(mcal) {plot.name <- paste0(plot.name, "_mitocal")
      } else {if(!(missing("cal.path"))) {plot.name <- paste0(plot.name, "_spikecal")}}
      if(smooth.mode == "VX"){plot.name <- paste0(plot.name, "_Fsize", Fsize)}
      if(smooth.mode == "HW"){plot.name <- paste0(plot.name, "_Hann", win)}
      if(y.mode == "auto") {plot.name <- paste0(plot.name, "_autoY")}
      if(y.mode == "fixed") {plot.name <- paste0(plot.name, "_ylim", fixed.y)}

      sf=5;if (length(samples) > 5){sf=length(samples)}


      if (plot.scale){
      if (out.mode == 'pdf') {pdf(file = paste0(plot.name, ".pdf"), width=plotdimension[1], height=plotdimension[2])
      }
      if (out.mode == 'png') {
        png(file = paste0(plot.name, ".png"), width=plotdimension[1], height=plotdimension[2],  units = "in", res = 400)
      }
      } else {
      if (out.mode == 'pdf') {pdf(file = paste0(plot.name, ".pdf"), width=plotdimension[1], height=(plotdimension[2]/5)*sf)
      }
      if (out.mode == 'png') {
        png(file = paste0(plot.name, ".png"), width=plotdimension[1], height=(plotdimension[2]/5)*sf,  units = "in", res = resolution)
      }
      }


      plotnumber=length(samples) # Number from 1 to 5

      if(basic.mode) {layout(matrix(rep(1:plotnumber,each = 3)))} else {layout(matrix(c(rep(1:plotnumber,each = 3), (plotnumber+1):(plotnumber+3))))}

      par(mar=c(0,7,0,0),oma = c(3, 1, 3, 1),las=1, mgp = c(2,1,0), bty = "n") # Sets margins per graph and outside margins per grouped set (order is bottom, left, top, right)

      if(!basic.mode) {
        # TSS.1 = TSS.df[J(chrom, xl1:xl2), nomatch = 0L]
      features.1 <- features[features$chr == chrom & features$start >(xl1-10000) & features$stop<(xl2+10000),]}

      ###############################################################################################################################################
      for (k in samples){
        #Subset for region of interest
        # sae3.0=DSBList[[k]][J(chrom, (xl1-(100*Fsize)):(xl2+(100*Fsize))), nomatch=0L]
        # sae3.0$Total <- sae3.0$Watson+sae3.0$Crick

        sae2.0=DSBList[[k]][J(chrom, xl1:xl2), nomatch=0L]
        sae2.0$Total <- sae2.0$Watson+sae2.0$Crick
        #Decompression code here
        sae2.1 <- data.frame(Chr=chrom, Pos=(xl1:xl2)) # Creates expanded empty dataframe with Chr and Pos locations
        sae2.1 <- merge(sae2.1,sae2.0, all=TRUE) # Merge expanded empty dataframe with compressed sae2.1 dataframe
        if (os != 0){sae2.1$Crick <- c(sae2.1$Crick[-(seq(os))], rep(NA, os))}  # OFFSET CRICK IF YOU WANT
        sae2.1[is.na(sae2.1)] <- 0 # Convert all NA values to zero
        sae2.1$Total=sae2.1$Watson+sae2.1$Crick

        # ########################################################################################################################################################
        # VatiableX Smoothing function #### New version creates two smoothed plots for each profile for overlaying
        #Window width (this is the number of values that will be averaged, NOT the  distance in bp).
        if (smooth.mode == "VX"){
          sae3.0=DSBList[[k]][J(chrom, (xl1-(100*Fsize)):(xl2+(100*Fsize))), nomatch=0L]
          sae3.0$Total <- sae3.0$Watson+sae3.0$Crick
          VX.out <- VX(sae3.0, Fsize = Fsize, norm.factor = norm.factor)
          b <- VX.out[[1]]
          d <- VX.out[[2]]
          e <- VX.out[[3]]
        }

        ##########################################################################################################################################################
        # Han Smoothing function #### temp and smooth are just two temporary vectors. New version creates two smoothed plots for each profile for overlaying
        # adjust this if needed when adjusting hann window smoothing
        if (smooth.mode == "HW") {
          hw=hanning.window(win) #create hanning window (require package e1071 to be loaded)
          scalar=1/sum(hw) #scalar [1]

          temp=NULL
          for (j in 3:5){
            temp=c(rep(0,win),sae2.1[1:nrow(sae2.1),j], rep(0,win)) # Create vector length of chromosome and extend by the length of the slidign window with zeros at both ends
            smooth=norm.factor*scalar*(stats::filter(temp,hw)) # smooth the temp vector using the hann window
            smooth=smooth[(win+1):(length(smooth)-win)] # trim smooth to correct lengthsmooth2=smooth2[(win2+1):(length(smooth2)-win2)] # trim smooth to correct length
            sae2.1[j+3]=smooth
          }

          colnames(sae2.1)=c("Chr", "Pos", "Watson", "Crick", "Total", "watson.s", "crick.s", "total.s")
        }


        if (disparity) {
          sae2.2 <- sae2.1
          hw=rep(1,win)
          # hw = hanning.window(win) #create hanning window (require package e1071 to be loaded)
          scalar=1/sum(hw) #scalar [1]

          temp=NULL
          for (j in 3:5){
            temp=c(rep(0,win),sae2.2[1:nrow(sae2.2),j], rep(0,win)) # Create vector length of chromosome and extend by the length of the slidign window with zeros at both ends
            smooth=norm.factor*scalar*(stats::filter(temp,hw)) # smooth the temp vector using the hann window
            smooth=smooth[(win+1):(length(smooth)-win)] # trim smooth to correct lengthsmooth2=smooth2[(win2+1):(length(smooth2)-win2)] # trim smooth to correct length
            sae2.2[j+3]=smooth
          }

          colnames(sae2.2)=c("Chr", "Pos", "Watson", "Crick", "Total", "watson.s", "crick.s", "total.s")
        }

        ##########################################################################################################################################################
        # Plot boundaries:
        # par(pty = "s")
        if(y.mode == "fixed") {ylims2 = c(-ylims,ylims)} # ylims for plotting
        if(y.mode == "auto"){ylims2 = c(-max(max(sae2.1$Watson/Mreads[k]), max(sae2.1$Crick/Mreads[k])),max(max(sae2.1$Watson/Mreads[k]), max(sae2.1$Crick/Mreads[k])))}
        sae2.1$Pos <- as.numeric(sae2.1$Pos)
        plot(sae2.1$Pos,sae2.1$Total/Mreads[k], type="n", xlim=c(xl1,xl2), ylim=ylims2, ylab=wrap.it(paste0(DSBListNames[k]," (HpM)"), 20),  axes=FALSE, frame.plot=TRUE, xaxs = "i", yaxs = "i",cex.lab = 0.75) #plot the start histogram
        Axis(side=2)
        # if(k == length(samples)) {Axis(side=1, at = pretty(c(xl1, xl2), 10), xaxs = "i", yaxs = "i")}

        # Broad Overlays of VariableX-smoothed data:
        if (smooth.mode == "VX") {
          xx <- c(b$PosAve, rev(b$PosAve))
          yy <- c(rep(0, nrow(b)), rev(b$WatsonAveNorm))
          polygon(xx, yy, col=rgb(240,128,128,alpha = 120, maxColorValue = 255), border=NA)

          xx2 <- c(d$PosAve, rev(d$PosAve))
          yy2 <- c(rep(0, nrow(d)), rev(d$CrickAveNorm))
          polygon(xx2, -1*yy2, col=rgb(173,216,230,alpha=120, maxColorValue = 255), border=NA)

          # xx3 <- c(e$PosAve, rev(e$PosAve))
          # yy3 <- c(rep(0, nrow(e)), rev(e$TotalAveNorm))
          # polygon(xx3, yy3+ylims2[1], col=rgb(100,100,100,alpha=100, maxColorValue = 255), border=NA)
        }
# abline(h=0)
if(origin.line) {points(xl1:xl2, rep(0,length(xl1:xl2)), pch = ".", col = "lightgrey")}
        if(loci.arrow) {arrows(loci.table$Pos[n], (-(ylims2[2]/3)), loci.table$Pos[n], (-(ylims2[2]/10)))}

        # Broad Overlays of han-smoothed
        if (smooth.mode == "HW"){
          lines(sae2.1$Pos,sae2.1$watson.s/Mreads[k], type="h", xlim=c(xl1,xl2), col=rgb(240,128,128,alpha = 120, maxColorValue = 255), lend = "butt") #plot the start histogram
          lines(sae2.1$Pos,-sae2.1$crick.s/Mreads[k], type="h", xlim=c(xl1,xl2), col=rgb(173,216,230,alpha=120, maxColorValue = 255), lend = "butt") #plot the start histogram
          # lines(sae2.1$Pos,(0.5*sae2.1$total.s/Mreads[k])+ylims2[1], type="l", xlim=c(xl1,xl2), col="grey") #plot the start histogram
        }

        # Overlay the raw data
        if (xl1 != 0) {
          lines(sae2.1$Pos,sae2.1$Watson/(Mreads[k]), type="h", xlim=c(xl1,xl2), col="red3", lwd=1, lend = "butt") #plot the start histogram
          lines(sae2.1$Pos,-sae2.1$Crick/(Mreads[k]), type="h", xlim=c(xl1,xl2), col="dodgerblue2", lwd=1, lend = "butt") #plot the start histogram
          # lines(sae2.1$Pos,sae2.1$Total/(Mreads[k])+ylim[1], type="l", xlim=c(xl1,xl2), col="grey27", lwd=0.5) #plot the start histogram
        }

        if (disparity){
          lines(sae2.2$Pos,(0.4*ylims2[2]*(log2(sae2.2$watson.s/sae2.2$crick.s))), type="l", xlim=c(xl1,xl2), col="black", lend = "butt") #plot the start histogram
          # lines(sae2.1$Pos,(0.5*sae2.1$total.s/Mreads[k])+ylims2[1], type="l", xlim=c(xl1,xl2), col="grey") #plot the start histogram
        }

        # plot the orientated TSS present in the region
        # if (nrow(TSS.1) > 0) { col.index <- cbind(c("#FF000080", NA), c("+", "-"))
        # col.vec <- col.index[match(TSS.1$Strand, col.index[,2]),1]
        # y.index <- cbind(c(ylims2[1], ylims2[1]), c("+", "-"))
        # y.vec <- as.numeric(y.index[match(TSS.1$Strand, y.index[,2]),1]) -ylims2[1]/10
        # Arrows(TSS.1$Pos, y.vec, TSS.1$Pos+1, y.vec, arr.length = 2/plotdimension[1], arr.type = "triangle", arr.col = col.vec, lcol = NA)
        # text(TSS.1$Pos, y.vec, labels = TSS.1$SYMBOL, cex = 0.5, col = "black", pos = 3)
        # col.index <- cbind(c(NA,"#0000FF80"), c("+", "-"))
        # col.vec <- col.index[match(TSS.1$Strand, col.index[,2]),1]
        # Arrows(TSS.1$Pos, y.vec, TSS.1$Pos-1, y.vec, arr.length = 2/plotdimension[1], arr.type = "triangle", arr.col = col.vec, lcol = NA)
        # }
      }
      title(main = paste(c(as.character(chrom), " : ", format(xl1,big.mark=",",scientific=FALSE), " - ", format(xl2,big.mark=",",scientific=FALSE), "  Fsize = ", Fsize), sep = " " , collapse = ""), outer=T)

      ##########################################################################################################################################################
      #Now plot the gene datatrack
      #First subset the relevant data
      # features <- subset(features,chr==chrom & start>(xl1-10000) & stop<(xl2+10000)) #Make a sub-table of ALLElements where chr = 1 and has limits just beyond plot range

      if(!basic.mode){
      genes <- features.1[features.1$type %in% c("gene", "ncRNA_gene"),] #Make a sub-table of ALLElements
      rRNA <- features.1[features.1$type == "rRNA_gene",] #Make a sub-table of ALLElements
      rRNA <- rRNA[rRNA$sysname != "RDN37-2" & rRNA$sysname != "RDN37-1",] #Make a sub-table of ALLElements

      tRNA <- features.1[features.1$type == "tRNA_gene",] #Make a sub-table of ALLElements
      snRNA <- features.1[features.1$type == "snRNA_gene",] #Make a sub-table of ALLElements
      MAS <- features.1[features.1$type == "matrix_attachment_site",] #Make a sub-table of ALLElements

      LTR_retro= features.1[features.1$type == "LTR_retrotransposon",] #Make a sub-table of ALLElements
      LTR= features.1[features.1$type == "long_terminal_repeat",] #Make a sub-table of ALLElements
      CEN <- features.1[features.1$type == "centromere",]
      ARS <- features.1[features.1$type == "ARS_consensus_sequence",]
      #Now perform the plot
      plot(sae2.1$Pos,sae2.1$Total, xaxt="n",yaxt="n",type="n", ylab=paste("Genes"),cex.lab=0.75,font=2, xlim=c(xl1,xl2),  xaxs = "i", yaxs = "i", ylim=c(-60,155),axes=F) #set up empty plot
      ##########################################################################################################################################################
      ########### STOP HERE IF YOU ARE PLOTTING WHOLE CHROMOSOMES!!! #############
      # ##########################################################################################################################################################
      # Following module draws arrows for each element
      # segments(centro[chrom.,"start"], 0, centro[chrom.,"stop"], lwd = 3)
      # text((centro[chrom.,"start"]+centro[chrom.,"stop"])/2,-20, font=3, paste0("CEN", chrom.), cex=0.9)
      if(gene.names == "genename"){
      if(gene.track){
        genesW=subset(genes,genename=="Dubious_ORF" & orientation =="+") #Make a sub-table of ALLElements
        featureplotter(genesW, col = "wheat", lty = 2, strand = "+", av = 75, pos1 = xl1, pos2 = xl2)
        genesW=subset(genes,genename !="Dubious_ORF" & orientation =="+") #Make a sub-table of ALLElements
        featureplotter(genesW, col = "wheat", lty = 1, strand = "+", av = 75, pos1 = xl1, pos2 = xl2)
        genesC=subset(genes,genename=="Dubious_ORF" & orientation =="-") #Make a sub-table of ALLElements
        featureplotter(genesC, col = "thistle", lty = 2, strand = "-", av = 25, pos1 = xl1, pos2 = xl2)
        genesC=subset(genes,genename !="Dubious_ORF" & orientation =="-") #Make a sub-table of ALLElements
        featureplotter(genesC, col = "thistle", lty = 1, strand = "-", av = 25, pos1 = xl1, pos2 = xl2)

        rRNAW=subset(rRNA,orientation =="+") #Make a sub-table of ALLElements
        featureplotter(rRNAW, col = "wheat", lty = 2, strand = "+", av = 75, tv = 100, pos1 = xl1, pos2 = xl2)
        rRNAC=subset(rRNA,orientation =="-") #Make a sub-table of ALLElements
        featureplotter(rRNAC, col = "thistle", lty = 2, strand = "-", av = 25, tv = 0, pos1 = xl1, pos2 = xl2)

        tRNAW=subset(tRNA,orientation =="+") #Make a sub-table of ALLElements
        featureplotter(tRNAW, col = "darkolivegreen1", lty = 2, strand = "+", av = 75, tv = 125, pos1 = xl1, pos2 = xl2)
        tRNAC=subset(tRNA,orientation =="-") #Make a sub-table of ALLElements
        featureplotter(tRNAC, col = "darkolivegreen1", lty = 2, strand = "-", av = 25, tv = -25, pos1 = xl1, pos2 = xl2)

        snRNAW=subset(snRNA,orientation =="+") #Make a sub-table of ALLElements
        featureplotter(snRNAW, col = "gold", lty = 2, strand = "+", av = 75, tv = 125, pos1 = xl1, pos2 = xl2)
        snRNAC=subset(snRNA,orientation =="-") #Make a sub-table of ALLElements
        featureplotter(snRNAC, col = "gold", lty = 2, strand = "-", av = 25, tv = -25, pos1 = xl1, pos2 = xl2)

        MASW=subset(MAS,orientation =="+") #Make a sub-table of ALLElements
        featureplotter(MASW, col = "slateblue1", lty = 1, strand = "*", av = 0, tv = -40, pos1 = xl1, pos2 = xl2)
      }
      }

      if(gene.names == "sysname"){
        if(gene.track){
          genesW=subset(genes,genename=="Dubious_ORF" & orientation =="+") #Make a sub-table of ALLElements
          featureplotter(genesW, col = "wheat", coln = "sysname", lty = 2, strand = "+", av = 75, pos1 = xl1, pos2 = xl2)
          genesW=subset(genes,genename !="Dubious_ORF" & orientation =="+") #Make a sub-table of ALLElements
          featureplotter(genesW, col = "wheat", coln = "sysname", lty = 1, strand = "+", av = 75, pos1 = xl1, pos2 = xl2)
          genesC=subset(genes,genename=="Dubious_ORF" & orientation =="-") #Make a sub-table of ALLElements
          featureplotter(genesC, col = "thistle", coln = "sysname", lty = 2, strand = "-", av = 25, pos1 = xl1, pos2 = xl2)
          genesC=subset(genes,genename !="Dubious_ORF" & orientation =="-") #Make a sub-table of ALLElements
          featureplotter(genesC, col = "thistle", coln = "sysname", lty = 1, strand = "-", av = 25, pos1 = xl1, pos2 = xl2)

          rRNAW=subset(rRNA,orientation =="+") #Make a sub-table of ALLElements
          featureplotter(rRNAW, col = "wheat", coln = "sysname", lty = 2, strand = "+", av = 75, tv = 100, pos1 = xl1, pos2 = xl2)
          rRNAC=subset(rRNA,orientation =="-") #Make a sub-table of ALLElements
          featureplotter(rRNAC, col = "thistle", coln = "sysname", lty = 2, strand = "-", av = 25, tv = 0, pos1 = xl1, pos2 = xl2)

          tRNAW=subset(tRNA,orientation =="+") #Make a sub-table of ALLElements
          featureplotter(tRNAW, col = "darkolivegreen1", coln = "sysname", lty = 2, strand = "+", av = 75, tv = 125, pos1 = xl1, pos2 = xl2)
          tRNAC=subset(tRNA,orientation =="-") #Make a sub-table of ALLElements
          featureplotter(tRNAC, col = "darkolivegreen1", coln = "sysname", lty = 2, strand = "-", av = 25, tv = -25, pos1 = xl1, pos2 = xl2)

          snRNAW=subset(snRNA,orientation =="+") #Make a sub-table of ALLElements
          featureplotter(snRNAW, col = "gold", coln = "sysname", lty = 2, strand = "+", av = 75, tv = 125, pos1 = xl1, pos2 = xl2)
          snRNAC=subset(snRNA,orientation =="-") #Make a sub-table of ALLElements
          featureplotter(snRNAC, col = "gold", coln = "sysname",lty = 2, strand = "-", av = 25, tv = -25, pos1 = xl1, pos2 = xl2)

          MASW=subset(MAS,orientation =="+") #Make a sub-table of ALLElements
          featureplotter(MASW, col = "slateblue1", coln = "sysname",  lty = 1, strand = "*", av = 0, tv = -40, pos1 = xl1, pos2 = xl2)
        }
      }



      if(ARS.track){
        ARSW=subset(ARS,orientation =="+") #Make a sub-table of ALLElements
        featureplotter(ARSW, col = "tomato3", coln = "sysname", lty = 1, strand = "*", av = 0, tv = -40, pos1 = xl1, pos2 = xl2)
        ARSC=subset(ARS,orientation =="-") #Make a sub-table of ALLElements
        featureplotter(ARSC, col = "royalblue3", coln = "sysname", lty = 1, strand = "*", av = 0, tv = -40,pos1 = xl1, pos2 = xl2)
      }

      if(Ty.track){
        LTRW=subset(LTR, orientation =="+") #Make a sub-table of ALLElements
        featureplotter(LTRW, col = "palegreen1", coln = "sysname", lty = 1, strand = "+", av = 125, pos1 = xl1, pos2 = xl2)
        LTR_retroW=subset(LTR_retro, orientation =="+") #Make a sub-table of ALLElements
        featureplotter(LTR_retroW, col = "palegreen1", coln = "sysname", lty = 1, strand = "+", av=75, pos1 = xl1, pos2 = xl2)
        LTR_retroC=subset(LTR_retro,orientation =="-") #Make a sub-table of ALLElements
        featureplotter(LTR_retroC, col = "seagreen1", coln = "sysname", lty = 1, strand = "-", av = 25, pos1 = xl1, pos2 = xl2)
        LTRC=subset(LTR,orientation =="-") #Make a sub-table of ALLElements
        featureplotter(LTRC, col = "seagreen1", coln = "sysname", lty = 1, strand = "-", av = -25, pos1 = xl1, pos2 = xl2)
      }

      featureplotter(CEN, coln="sysname",col = "snow4", lty = 1, strand = "*", av = 0, tv= -30, pos1 = xl1, pos2 = xl2)

      # Now perform the plot
      if(MNase.track == T){
        MNase.1 <- MNase[J(chrom, xl1:xl2), nomatch=0L]
        if(nrow(MNase.1) != 0){
          plot(MNase.1$Pos,MNase.1$MNaseSignal, col = "lightsteelblue3", xaxt="n",yaxt="n",type="h", ylab=Nucleosomes,cex.lab=0.75,font=2, xlim=c(xl1,xl2),  xaxs = "i", yaxs = "i", axes= F) #set up empty plot
          Axis(side=2, xaxs = "i", yaxs = "i")
        }
      }
      ##########################################################################################################################################################


      if(scc1.track == T){
        scc1.1 <- scc1[J(chrom, xl1:xl2), nomatch=0L]
        plot(scc1.1$Pos,scc1.1$Value, col = "tomato2", xaxt="n",yaxt="n",type="h", ylab=paste(Cohesin),cex.lab=0.75,font=2, xlim=c(xl1,xl2),  xaxs = "i", yaxs = "i", axes= F) #set up empty plot
        Axis(side=1, at = pretty(c(xl1, xl2), 10), xaxs = "i", yaxs = "i")
        Axis(side=2, xaxs = "i", yaxs = "i")
      }



      }

    Axis(side=1, at = pretty(c(xl1, xl2), 10), xaxs = "i", yaxs = "i")
    Axis(side=2, xaxs = "i", yaxs = "i")
    dev.off()
      }
  setwd(path.name)
}

