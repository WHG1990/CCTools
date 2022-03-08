#' @title CCBin
#' @description Bins yeast or human CC-seq data quickly using data.table. The outpout can be strand separated, or totalled, and can be in sparse or full format.
#' @param in.mode Reading in .RDS file or a folder of FullMaps? Choose "rds" or "fullmap".  Defaults to "rds".
#' @param path.name A character string defining the RDS file or folder containing FullMaps.
#' @param exp.name Manually specify a name for the experiment.
#' @param genome Which genome has the data been aligned to? "Cer3H4L2", "W303" or "hg19". Defaults to "Cer3H4L2".
#' @param bin.width What size bins, in bp? Defaults to 100.
#' @param os Do you want to offset the Crick strand relative to the Watson strand prior to binning? Choose 3 for Top2, or 1 for Spo11. Choose 0 for no ofsetting. Defaults to 3.
#' @param out.mode1 "sparse" or "full" format. Defaults to "full". Defaults to "full"
#' @param out.mode2 Do you want to sum the Watson and Crick strands ("total") or report them seperately ("strand")? Defaults to "strand".
#' @param plot.mode Do you want to plot the binned chromosome maps as "pdf" of "png"? Choose "none" for no plots
#' @param ylims A single number specifying the maximum limit of the y axis, for both Watson and Crick. OR "auto", to automatically define for the entire set.
#' @param gene.track Do you want to plot the gene track? This is pretty pointless on whole chromosome plots. Awaiting a more sensible implementation. Defaults to FALSE.
#' @param bg Choose the background for the plots. Can be "none" (transparent), or "white".
#' @param mcal normalise based of the mitochondria
#' @param rDNAcal normalise based of the rDNA
#' @examples
#' @author Will Gittens
#' @import data.table
#' @export
CCBin <- function(in.mode = "fullmap", path.name = getwd(), exp.name, genome = "Cer3H4L2", bin.width = 100, os = 0, mcal = F,rDNAcal=F, cal.path, out.mode1 = "full", out.mode2 = "strand", plot.mode = "png", violin = T, ylims ="auto", ARS.track = T, cen.track = T,names, plotdimensions = "wide", bg = "white", samples, lin.mod = F) {

  if(!(out.mode1 %in% c("sparse", "full"))) { stop("Error: out.mode1 must be one of 'sparse' or 'full'") }
  if(!(out.mode2 %in% c("strand", "total"))) { stop("Error: out.mode2 must be one of 'strand' or 'total'") }

  library("data.table")
  library("magick")
  library("ggplot2")
  if(!missing(names)){DSBListNames <- names}
  if(ylims == "auto") {y.mode <- "auto"} else {y.mode <- "fixed"}
  if(y.mode == "fixed") {fixed.y <- ylims}

  sett.name <- paste0((bin.width/1000),"Kbp_binned_OS", os, "_", out.mode1, "_", out.mode2)
  if(mcal) {sett.name <- paste0(sett.name, "_mitocal")
  } else {if(!(missing("cal.path"))) {sett.name <- paste0(sett.name, "_spikecal")}}
  if(rDNAcal) {sett.name <- paste0(sett.name, "_rDNAcal")
  } else {if(!(missing("cal.path"))) {sett.name <- paste0(sett.name, "_spikecal")}}

  parent.directory <- dirname(path.name)

  if(in.mode == "rds") {
    load(path.name)
    exp.name <- substr(basename(path.name), 0, nchar(basename(path.name))-4 )
  }
  if(in.mode == "fullmap") {
    parent.directory <- dirname(path.name)
    CCList(path.name)
    exp.name <- exp.name
  }
  DSBList <- lapply(DSBList, data.table)

  if(missing(samples)){samples <- 1:length(DSBList)}
  if(samples[1] == "all"){samples <- 1:length(DSBList)}
  DSBList <- DSBList[samples]
  DSBListNames <- DSBListNames[samples]
  Mreads <- Mreads[samples]

  ####

  if(mcal) {DSBList <- mitoCCalibrate2(DSBList)
  } else {if(!(missing("cal.path"))) {DSBList <- CCalibrate(DSBList, Mreads_X = Mreads, Mreads_C = MreadR(cal.path))}}
  if(rDNAcal) {DSBList <- rDNACCalibrate(DSBList)
  } else {if(!(missing("cal.path"))) {DSBList <- CCalibrate(DSBList, Mreads_X = Mreads, Mreads_C = MreadR(cal.path))}}
  ############## Import Annotations
  if(genome == "Cer3H4L2") {chroms <- 1:18; chromlengths <- CCAnnotate(Cer3H4L2_chromlengths); centro <- CCAnnotate(Cer3H4L2_centro); features <- CCAnnotate(Cer3H4L2_AllElementsDUB)}
  if(genome == "W303") {chroms <- 1:18; chromlengths <- CCAnnotate(W303_chromlengths); centro <- CCAnnotate(W303_centro); features <- CCAnnotate(Cer3H4L2_AllElementsDUB)}
  if(genome == "human") {chroms <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22", "chrX", "chrM"); chromlengths <- CCAnnotate(hg19_chromlengths); centro <- CCAnnotate(hg19_centro)}
  if(genome == "pombase220208") {chroms <- 1:4; chromlengths <- CCAnnotate(pombase220208_chromlengths); centro <- CCAnnotate(pombase220208_centro); features <- CCAnnotate(pombase220208_AllElements)}
  if(genome =="WS220"){chroms <- 1:7; chromlengths <- CCAnnotate(WS220_chromlengths)}

  #offset if you like
  if(os != 0) {
    DSBList <- lapply(DSBList, function(y){
      temp.w <- y[y$Watson > 0,]
      temp.c <- y[y$Crick > 0,]
      temp.w$Crick <- NULL
      temp.c$Watson <- NULL
      temp.c$Pos <- temp.c$Pos - os
      z <- merge(temp.w, temp.c, by = c("Chr", "Pos"), all = T)
      z[is.na(z)] <- 0
      z$Total = z$Watson + z$Crick
      return(z)
    })
  }

    DSBList <- lapply(DSBList, function(y){
      y$Total = y$Watson + y$Crick
      return(y)
    })

  binned.DSBList <- lapply(DSBList, function(x){
    x$Pos <- roundup(x$Pos, bin.width)
    return(x)
  })

  rm(DSBList)
  gc()
  binned.DSBList <- lapply(binned.DSBList, function(x){
    if(out.mode2 == "total"){x <- x[,.(Total = sum(Total)),by=.(Chr,Pos)]}
    if(out.mode2 == "strand"){x <- x[,.(Watson = sum(Watson), Crick = sum(Crick)),by=.(Chr,Pos)]}
    return(x)
  })
  # If you want full (not sparse) output, this creates an empty binned genome with which you can merge the binned data.
  if(out.mode1 == "full"){
    chromlengths$RoundDown<-rounddown(chromlengths$Length, bin.width)
    binlist<-list()
    for (chrom in chroms) {
      tempchromlength<- chromlengths[which(chromlengths$Chr == chrom), ]
      rmax<-tempchromlength$RoundDown
      bintemp<-(seq(0,(rmax), by=bin.width))
      binlist[[length(binlist)+1]] <- bintemp
    }
    names(binlist)<- chroms
    binlist <- lapply(binlist, function(x){
      x <- x[-1]
      return(x)
    })
    binlist <-  lapply(binlist, function(x){
      x <- data.frame(x)
    })
    binlist <- mapply(function(x,y) {
      x <- cbind(x, rep(y, nrow(x)))
      return(x)
    }
    , binlist, chroms, SIMPLIFY = F)
    bintable <- do.call(rbind, binlist)
    colnames(bintable) <- c("Pos", "Chr")
    bintable <- bintable[,c(2,1)]
    bintable <- data.table(bintable); setkey(bintable, Chr, Pos)
    binned.DSBList <- lapply(binned.DSBList, function(x){
      setkey(x, Chr, Pos)
      x <- merge(x, bintable, all.x =F, all.y = T)
      x[is.na(x)] <- 0L
      return(x)
    })
  }
  setwd(parent.directory)
  ifelse(!dir.exists(file.path(parent.directory, "Binned")), dir.create(file.path(parent.directory, "Binned")), FALSE); setwd(file.path(parent.directory, "Binned")); directory2 <- getwd()
  ifelse(!dir.exists(file.path(directory2, "/",exp.name)), dir.create(file.path(directory2, "/",exp.name)), FALSE); setwd(file.path(directory2, "/", exp.name)); directory2 <- getwd()
  ifelse(!dir.exists(file.path(directory2, "/",sett.name)), dir.create(file.path(directory2, "/",sett.name)), FALSE); setwd(file.path(directory2, "/", sett.name)); directory2 <- getwd()
  ifelse(!dir.exists(file.path(directory2, "Binned_FullMaps")), dir.create(file.path(directory2, "Binned_FullMaps")), FALSE); setwd(file.path(directory2, "Binned_FullMaps")); directory3 <- getwd()
  save(binned.DSBList, DSBListNames, nfiles, chroms, bin.width, Mreads, file = paste0(exp.name,"_", sett.name, ".rds"))
  CCWrite(binned.DSBList, filenames = paste0(DSBListNames,"_", sett.name))

  # Pairwise correlations. This will write a table containing all pairwise pearson correlation values between samples in the binned.DSBList
  if(out.mode1 == "full" && out.mode2 == "total" && length(DSBListNames) > 2) {
    ifelse(!dir.exists(file.path(directory2, "/Correlation")), dir.create(file.path(directory2, "/Correlation")), FALSE); setwd(file.path(directory2, "/Correlation")); directory3 <- getwd()
    ifelse(!dir.exists(file.path(directory3, "/Total")), dir.create(file.path(directory3, "/Total")), FALSE); setwd(file.path(directory3, "/Total"))

    names(binned.DSBList) <- DSBListNames
    write.table(outer(binned.DSBList, binned.DSBList, FUN = Vectorize(function(x,y) {cor(x$Total, y$Total, method = "pearson")})), file = paste0(exp.name, "_Total_pearson_matrix.txt"), sep = "\t", col.names = DSBListNames, row.names = DSBListNames)

    if(lin.mod) {write.table(outer(binned.DSBList, binned.DSBList, FUN = Vectorize(function(x,y) {
      linmod <- lm((x$Watson+x$Crick) ~ (y$Watson+y$Crick))
      return(summary(linmod)$r.squared)
            }
      )), file = paste0(exp.name, "_Total_lm_rsquared_matrix.txt"), sep = "\t", col.names = DSBListNames, row.names = DSBListNames)
    }

    test <- binned.DSBList
    test <- lapply(test, function(x){
      x$Chr <- NULL
      x$Pos <- NULL
      return(x)
    })
    #plot dendogram
    test <- do.call("cbind", test)
    hc <- hclust(as.dist(1-cor(test, method="pearson")), method="complete")
    png(file = paste0(exp.name,"_Total_pearson_dendogram.png"), width=10, height=20,  units = "in", res = 400)
    plot(hc)
    dev.off()
    library(magick)
    image1 <- image_read(paste0(exp.name,"_Total_pearson_dendogram.png"))
    image_rotate(image1, 90) %>% image_write(paste0(exp.name,"_Total_pearson_dendogram.png"))
  }



  #strand mode plots all correlations (including Total). This should be tidied up at some point...
  if(out.mode1 == "full" && out.mode2 == "strand" && length(DSBListNames) > 2) {
    names(binned.DSBList) <- DSBListNames

    ifelse(!dir.exists(file.path(directory2, "/Correlation")), dir.create(file.path(directory2, "/Correlation")), FALSE); setwd(file.path(directory2, "/Correlation")); directory3 <- getwd()
    ifelse(!dir.exists(file.path(directory3, "/Total")), dir.create(file.path(directory3, "/Total")), FALSE); setwd(file.path(directory3, "/Total"))
    #TOTAL
    write.table(outer(binned.DSBList, binned.DSBList, FUN = Vectorize(function(x,y) {cor((x$Watson+x$Crick), (y$Watson+y$Crick), method = "pearson")})), file = paste0(exp.name, "_Total_pearson_matrix.txt"), sep = "\t", col.names = DSBListNames, row.names = DSBListNames)
    if(lin.mod) {write.table(outer(binned.DSBList, binned.DSBList, FUN = Vectorize(function(x,y) {
      linmod <- lm((x$Watson+x$Crick) ~ (y$Watson+y$Crick))
      return(summary(linmod)$r.squared)
    }
    )), file = paste0(exp.name, "_Total_lm_rsquared_matrix.txt"), sep = "\t", col.names = DSBListNames, row.names = DSBListNames)
    }
    test <- binned.DSBList
    test <- lapply(test, function(x){
      x$Total <- x$Watson+x$Crick
      x$Chr <- NULL
      x$Pos <- NULL
      x$Crick <- NULL
      x$Watson <- NULL
      return(x)
    })
    #plot dendogram
    test <- do.call("cbind", test)
    hc <- hclust(as.dist(1-cor(test, method="pearson")), method="complete")
    png(file = paste0(exp.name,"_Total_pearson_dendogram.png"), width=10, height=20,  units = "in", res = 400)
    plot(hc)
    dev.off()
    library(magick)
    image1 <- image_read(paste0(exp.name,"_Total_pearson_dendogram.png"))
    image_rotate(image1, 90) %>% image_write(paste0(exp.name,"_Total_pearson_dendogram.png"))

    #WATSON
    ifelse(!dir.exists(file.path(directory3, "/Watson")), dir.create(file.path(directory3, "/Watson")), FALSE); setwd(file.path(directory3, "/Watson"));
    write.table(outer(binned.DSBList, binned.DSBList, FUN = Vectorize(function(x,y) {cor(x$Watson, y$Watson, method = "pearson")})), file = paste0(exp.name, "_Watson_pearson_matrix.txt"), sep = "\t", col.names = DSBListNames, row.names = DSBListNames)
    if(lin.mod) {write.table(outer(binned.DSBList, binned.DSBList, FUN = Vectorize(function(x,y) {
      linmod <- lm(x$Watson ~ y$Watson)
      return(summary(linmod)$r.squared)
    }
    )), file = paste0(exp.name, "_Watson_lm_rsquared_matrix.txt"), sep = "\t", col.names = DSBListNames, row.names = DSBListNames)
    }

    test <- binned.DSBList
    test <- lapply(test, function(x){
      x$Chr <- NULL
      x$Pos <- NULL
      x$Crick <- NULL
      return(x)
    })
    #plot dendogram
    test <- do.call("cbind", test)
    hc <- hclust(as.dist(1-cor(test, method="pearson")), method="complete")
    png(file = paste0(exp.name,"_Watson_pearson_dendogram.png"), width=10, height=20,  units = "in", res = 400)
    plot(hc)
    dev.off()
    library(magick)
    image1 <- image_read(paste0(exp.name,"_Watson_pearson_dendogram.png"))
    image_rotate(image1, 90) %>% image_write(paste0(exp.name,"_Watson_pearson_dendogram.png"))


  # CRICK
    ifelse(!dir.exists(file.path(directory3, "/Crick")), dir.create(file.path(directory3, "/Crick")), FALSE); setwd(file.path(directory3, "/Crick"));
    write.table(outer(binned.DSBList, binned.DSBList, FUN = Vectorize(function(x,y) {cor(x$Crick, y$Crick, method = "pearson")})), file = paste0(exp.name, "_Crick_pearson_matrix.txt"), sep = "\t", col.names = DSBListNames, row.names = DSBListNames)
    if(lin.mod) {write.table(outer(binned.DSBList, binned.DSBList, FUN = Vectorize(function(x,y) {
      linmod <- lm(x$Crick ~ +y$Crick)
      return(summary(linmod)$r.squared)
    }
    )), file = paste0(exp.name, "_Crick_lm_rsquared_matrix.txt"), sep = "\t", col.names = DSBListNames, row.names = DSBListNames)
    }
    test <- binned.DSBList
    test <- lapply(test, function(x){
      x$Chr <- NULL
      x$Pos <- NULL
      x$Watson <- NULL
      return(x)
    })
    #plot dendogram
    test <- do.call("cbind", test)
    hc <- hclust(as.dist(1-cor(test, method="pearson")), method="complete")
    png(file = paste0(exp.name,"_Crick_pearson_dendogram.png"), width=10, height=20,  units = "in", res = 400)
    plot(hc)
    dev.off()
    library(magick)
    image1 <- image_read(paste0(exp.name,"_Crick_pearson_dendogram.png"))
    image_rotate(image1, 90) %>% image_write(paste0(exp.name,"_Crick_pearson_dendogram.png"))

}


  #plot violin plots of distribution
  if (violin){
    ifelse(!dir.exists(file.path(directory2, "/Violins")), dir.create(file.path(directory2, "/Violins")), FALSE); setwd(file.path(directory2, "/Violins"))
  plotdata <- lapply(binned.DSBList, function(x){
    x$Total <- x$Watson + x$Crick
    return(x)
  })
  plotdata <- mapply(function(x,y) {
    x$Condition <- y
    return(x)
  }, plotdata, DSBListNames, SIMPLIFY = F)

  plotdata <- mapply(function(x,y) {
    x$HpMpbp <-( x$Total/y)/bin.width
    return(x)
  }, plotdata, Mreads, SIMPLIFY = F)
  plotdata <- do.call("rbind", plotdata)
  plotdata <- plotdata[plotdata$Chr < 17,]
  plotdata <- plotdata[!(plotdata$Chr == 12 & plotdata$Pos >= 451000 & plotdata$Pos <= 469000), ]
  plotdata <- plotdata[!(plotdata$Chr == 8 & plotdata$Pos >= 211000 & plotdata$Pos <= 216000), ]
  plotdata[plotdata$HpM == max(plotdata$HpM),]
  # my.cols <- c("#66C2A5", "#F06B42")
    p <- ggplot(plotdata, aes(x=Condition, y=HpMpbp, fill = Condition)) +
    geom_violin() +
    # scale_fill_manual(values=my.cols) +
    scale_y_log10() +
    geom_pointrange(aes(fill = Condition),
                    stat = "summary",
                    fun.ymin = function(z) {quantile(z,0.25)},
                    fun.ymax = function(z) {quantile(z,0.75)},
                    fun.y = median,
                    position = position_dodge(width = 0.9),
                    size = 0.2,
                    color = "black") +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  p
  ggsave(paste0(exp.name, "_Total_Violin.pdf"), width = 30, height = 30, units = "cm")
}

  # browser()
    if(plot.mode != "none") {
    # optional plotting module below. This will plot all chroms.
    vec <- c(66.2,46.8,33.1,23.4,16.5,11.7,8.3,5.8,4.1,2.9,2.0,1.5,1, 7.5, 4.6, 3.3, 3.75)
    vec2 <- c(93.6,66.2,46.8,33.1,23.4,16.5,11.7,8.3,5.8,4.1,2.9,2.0,1.5, 13.3, 8.3, 5.8, 6.65)
    vec3 <- c("4A0", "2A0", "A0", "A1", 'A2', 'A3', 'A4', 'A5', 'A6', "A7", 'A8', 'A9', 'A10', "wide", "wide2", "wide3", "wide4")
    papersize <- data.frame(vec3,vec2,vec)
    plotdimension <- as.numeric(papersize[match(plotdimensions, vec3),2:3])

    setwd(directory2)
    sf=5;if (length(DSBListNames) > 5){sf=length(DSBListNames)}
    ifelse(!dir.exists(file.path(directory2, "ChromosomePlots")), dir.create(file.path(directory2, "ChromosomePlots")), FALSE); setwd(file.path(directory2, "ChromosomePlots")); directory3 <- getwd()
    for(i in 1:length(chroms)){

      plot.name <- paste0("Chr", i, "_", exp.name, "_", sett.name)
      if(y.mode == "auto") {plot.name <- paste0(plot.name, "_autoY")}
      if(y.mode == "fixed") {plot.name <- paste0(plot.name, "_ylim", fixed.y)}

      # png(file = paste0(plot.name, ".png"), width=plotdimension[1], height=(plotdimension[2]/5)*sf,  units = "in", res = 400)

      if(plot.mode == "pdf"){pdf(file = paste0(plot.name, ".pdf"), width=plotdimension[1], height=(plotdimension[2]/5)*sf)}
      if(plot.mode == "png"){png(file = paste0(plot.name, ".png"), width=plotdimension[1], height=(plotdimension[2]/5)*sf,  units = "in", res = 400)}

      # png(file = paste0(plot.name, ".png"), width=plotdimension[1], height=(plotdimension[2]/5)*sf,  units = "in", res = 600)
      plotnumber=length(binned.DSBList) # Number from 1 to 5
      # if(ARS.track) {layout(matrix(c(rep(1:plotnumber,each = 3), (plotnumber+1):(plotnumber+3))))}

      layout(matrix(c(rep(1:plotnumber,each = 3),(plotnumber+1) )))

      if(bg == "none"){par(mar=c(0,7,0,0),oma = c(1, 1, 1, 1),mgp = c(3,1,0),  bty = "n", bg = NA)} # Sets margins per graph and outside margins per grouped set (order is bottom, left, top, right)
      if(bg == "white"){par(mar=c(0,7,0,0),oma = c(1, 1, 1, 1),mgp = c(3,1,0),  bty = "n")} # Sets margins per graph and outside margins per grouped set (order is bottom, left, top, right)

      data.chrom <- lapply(binned.DSBList, function(x){
        x <- x[x$Chr == chroms[i],]
        return(x)
      })

      if("Total" %in% colnames(data.chrom[[1]])) {ylims <- max(mapply(function(x,y) {max(x$Total/ (y*bin.width))}, data.chrom, Mreads))}
      else {ylims <- max(mapply(function(x,y) {max(x$Watson/ (y*bin.width))}, data.chrom, Mreads))}
      if(y.mode == "fixed") { ylims <- fixed.y}
      xbreaks <- pretty(data.chrom[[1]]$Pos, n = 10)
      for(y in 1:length(data.chrom)){
        if("Total" %in% colnames(data.chrom[[1]])) {plot(data.chrom[[y]]$Pos, data.chrom[[y]]$Total/ (Mreads[y]*bin.width),  ylim = c(0,ylims), ylab=wrap.it(paste0(DSBListNames[y]," (HpM)"), 20), xlab ="",  axes=FALSE, frame.plot=TRUE, xaxs = "i", yaxs = "i",cex.lab = 0.6, type = "h", col = "dodgerblue3")}
        else{
          plot(data.chrom[[y]]$Pos, data.chrom[[y]]$Watson/(Mreads[y]*bin.width),  ylim = c(-ylims,ylims), ylab=wrap.it(paste0(DSBListNames[y]," (HpM)"), 20), lend=1,axes=FALSE, frame.plot=TRUE, xaxs = "i", yaxs = "i",cex.lab = 0.6, type = "h", col = "red3")
          lines(data.chrom[[y]]$Pos, -data.chrom[[y]]$Crick/(Mreads[y]*bin.width),  ylim = c(-ylims,ylims), ylab=wrap.it(paste0(DSBListNames[y]," (HpM)"), 20),  lend=1,axes=FALSE, frame.plot=TRUE, xaxs = "i", yaxs = "i",cex.lab = 0.6, type = "h", col = "dodgerblue2")
          # polygon2(data.chrom[[y]]$Pos, data.chrom[[y]]$Watson/(Mreads[y]*bin.width),  col = "red3")
          # polygon2(data.chrom[[y]]$Pos, -data.chrom[[y]]$Crick/(Mreads[y]*bin.width),  col = "dodgerblue2")

        }
        Axis(side = 2)
        if(y == length(data.chrom)) {
          Axis(side=1, at = xbreaks, xaxs = "i", yaxs = "i")
        }
      }
      mtext(text=paste0("Position on chromosome ", chroms[i], " / bp"),side=1,line=0,outer=TRUE,cex=0.8)

      # genetrack
      ################# gene loading
      plot(data.chrom[[y]]$Pos, data.chrom[[y]]$Watson, xaxt="n",yaxt="n",type="n", ylab=paste(" "),cex.lab=1.5,font=2,  xaxs = "i", yaxs = "i", ylim=c(-100,120),axes=F) #set up empty plot
      # Axis(side=1, at = pretty(c(xl1, xl2), 10), xaxs = "i", yaxs = "i")
      if(cen.track){
        if(i <= nrow(centro)){segments(as.numeric(centro[centro$Chr == chroms[i],"chromStart"]), 0, as.numeric(centro[centro$Chr == chroms[i],"chromEnd"]), 0, lwd = 3)}
      }
      if(ARS.track){
        features.1 <- features[features$chr == i,]
        ARS <- features.1[features.1$type == "ARS_consensus_sequence",]
        if(nrow(ARS) > 0){
          ARS$orientation = "*"
          featureplotter(ARS, col = "tomato3", coln = "sysname", lty = 1, strand = "*", av = 50, tv = 75, pos1 = 1, pos2 = chromlengths[i,"Length"], lab = F)
        }
      }

      dev.off()
    }
  }
  setwd(path.name)
  return(binned.DSBList)
}

