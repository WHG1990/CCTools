#' @title CCFold
#' @description Bins, loess smooths and calculates folds of yeast or human CC-seq data quickly using data.table.
#' @param in.mode Reading in .RDS file or a folder of FullMaps? Choose "rds" or "fullmap".  Defaults to "rds".
#' @param s.dir A character string defining the RDS file or folder containing the SAMPLE FullMaps.
#' @param c.dir A character string defining the RDS file or folder containing the CONTROL FullMaps.
#' @param exp.name Manually specify a name for the experiment.
#' @param genome Which genome has the data been aligned to? "Cer3H4L2", "W303" or "hg19". Defaults to "Cer3H4L2".
#' @param bin.width What size bins, in bp? Defaults to 100.
#' @param sc Parameter which changes the "span" of the smooth. Higher numbers will give a broader smooth.
#' @param si Parameter which alters the final number of points of the smooth curve to be plotted. Defaults to 100 bp.
#' @param os Do you want to offset the Crick strand relative to the Watson strand prior to binning? Choose 3 for Top2, or 1 for Spo11. Choose 0 for no ofsetting. Defaults to 3.
#' @param out.mode1 "sparse" or "full" format. Defaults to "full". Defaults to "full"
#' @param out.mode2 Do you want to sum the Watson and Crick strands ("total") or report them seperately ("strand")? Defaults to "strand".
#' @param plot.mode Do you want to plot the binned chromosome maps? (T/F). Defaults to TRUE.
#' @param ylims A single number specifying the maximum limit of the y axis, for both Watson and Crick. OR "auto", to automatically define for the entire set.
#' @param gene.track Do you want to plot the gene track? This is pretty pointless on whole chromosome plots. Awaiting a more sensible implementation. Defaults to FALSE.
#' @param names Can specify your own vector of alias names for the libraries.
#' @examples
#' @author Will Gittens
#' @import data.table
#' @export
CCFold <- function(in.mode = "fullmap", s.dir = getwd(), c.dir,  exp.name, genome = "Cer3H4L2", bin.width = 100, sc = 80, si = 100, os = 0, mcal = F,rDNAcal=F, cal.path, out.mode1 = "full", out.mode2 = "strand", plot.mode = T, ylims ="auto", gene.track = F, ARS.track = T,names) {
library("data.table")
library("magick")
library("RColorBrewer")
if(ylims == "auto") {y.mode <- "auto"} else {y.mode <- "fixed"}
if(y.mode == "fixed") {fixed.y <- ylims}

sett.name <- paste0((bin.width/1000),"Kbp_binned_OS", os, "_SC", sc, "_SI", si, out.mode1, "_", out.mode2)
if(mcal) {sett.name <- paste0(sett.name, "_mitocal")
} else {if(!(missing("cal.path"))) {sett.name <- paste0(sett.name, "_spikecal")}}
if(rDNAcal) {sett.name <- paste0(sett.name, "_rDNAcal")
} else {if(!(missing("cal.path"))) {sett.name <- paste0(sett.name, "_spikecal")}}

parent.directory <- dirname(s.dir)

if(in.mode == "rds") {
  load(s.dir)
  exp.name <- substr(basename(s.dir), 0, nchar(basename(s.dir))-4 )
}
  if(in.mode == "fullmap") {
    parent.directory <- dirname(s.dir)
    CCList(s.dir)
    exp.name <- exp.name
  }
DSBList <- lapply(DSBList, data.table)

####

if(mcal) {DSBList <- mitoCCalibrate2(DSBList)
} else {if(!(missing("cal.path"))) {DSBList <- CCalibrate(DSBList, Mreads_X = Mreads, Mreads_C = MreadR(cal.path))}}
if(rDNAcal) {DSBList <- rDNACCalibrate(DSBList)
} else {if(!(missing("cal.path"))) {DSBList <- CCalibrate(DSBList, Mreads_X = Mreads, Mreads_C = MreadR(cal.path))}}
############## Import Annotations
if(genome == "Cer3H4L2") {chroms <- 1:18; chromlengths <- CCAnnotate(Cer3H4L2_chromlengths); centro <- CCAnnotate(Cer3H4L2_centro); features <- CCAnnotate(Cer3H4L2_AllElementsDUB)}
if(genome == "W303") {chroms <- 1:18; chromlengths <- CCAnnotate(W303_chromlengths); centro <- CCAnnotate(W303_centro); features <- CCAnnotate(Cer3H4L2_AllElementsDUB)}
if(genome == "human") {c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22", "chrX", "chrY", "chrM"); chroms <- chromlengths <- CCAnnotate((hg19_chromlengths))}

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
return(z)
})
}

binned.DSBList <- lapply(DSBList, function(x){
  x$Total = x$Watson + x$Crick
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

samp.DSBList <- binned.DSBList
samp.DSBListNames <- DSBListNames

#### now for controls

if(in.mode == "rds") {
  load(c.dir)
  exp.name <- substr(basename(s.dir), 0, nchar(basename(s.dir))-4 )
}
if(in.mode == "fullmap") {
  CCList(c.dir)
  exp.name <- exp.name
}
DSBList <- lapply(DSBList, data.table)


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
    return(z)
  })
}

binned.DSBList <- lapply(DSBList, function(x){
  x$Total = x$Watson + x$Crick
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

ctrl.DSBList <- binned.DSBList
ctrl.DSBListNames <- DSBListNames


# folds
samp.Mreads <- MreadR(samp.DSBList)
ctrl.Mreads <- MreadR(ctrl.DSBList)

samp.DSBList <- mapply(function(x,y){
  x$samp.HpM <- x$Total/y
  return(x)
}, samp.DSBList, samp.Mreads, SIMPLIFY = F)

ctrl.DSBList <- mapply(function(x,y){
  x$ctrl.HpM <- x$Total/y
  return(x)
}, ctrl.DSBList, ctrl.Mreads, SIMPLIFY = F)

out.DSBList <- rep(list(samp.DSBList), length(ctrl.DSBList))

# out.DSBList2 <- out.DSBList
# out.DSBList <- out.DSBList2
DSBListNames <- c()
setwd(parent.directory)
ifelse(!dir.exists(file.path(parent.directory, "Fold")), dir.create(file.path(parent.directory, "Fold")), FALSE); setwd(file.path(parent.directory, "Fold")); directory2 <- getwd()
for (c in 1:length(ctrl.DSBList)){
  for (i in 1:length(samp.DSBList)){
      out.DSBList[[c]][[i]]$log2fold = log2(samp.DSBList[[i]]$samp.HpM/ctrl.DSBList[[c]]$ctrl.HpM)
      out.DSBList[[c]][[i]] = out.DSBList[[c]][[i]][is.finite(out.DSBList[[c]][[i]]$log2fold),]
      write.table(out.DSBList[[c]][[i]],file=paste0(samp.DSBListNames[[i]],"_OVER_", ctrl.DSBListNames[[c]],"_",sett.name, ".txt"),quote = FALSE, sep="\t",row.names = FALSE,)
  }
}


DSBListNames <- paste(rep(samp.DSBListNames, length(ctrl.DSBListNames)), rep(ctrl.DSBListNames, each = length(samp.DSBListNames)), sep = "_OVER_")



out.DSBList <- do.call("c", out.DSBList)
names(out.DSBList) <- DSBListNames
out.DSBList <- out.DSBList


if(plot.mode) {
  # optional plotting module below. This will plot all chroms.
  if(!missing(names)){DSBListNames <- names}
  setwd(directory2)

  ifelse(!dir.exists(file.path(directory2, "ChromosomePlots")), dir.create(file.path(directory2, "ChromosomePlots")), FALSE); setwd(file.path(directory2, "ChromosomePlots")); directory3 <- getwd()
  for(i in 1:length(chroms)){
    features.1 <- features[features$chr == i,]
    plot.name <- paste0("Chr", i, "_", exp.name, "_", sett.name)
    if(y.mode == "auto") {plot.name <- paste0(plot.name, "_autoY")}
    if(y.mode == "fixed") {plot.mode <- paste0(plot.name, "_ylim", fixed.y)}

    # png(file = paste0(plot.name,".png"), width=9.75, height=6.1,  units = "in", res = 400)
    pdf(file = paste0(plot.name,".pdf"), width=9.75, height=6.1)


    plotnumber=length(out.DSBList) # Number from 1 to 5
    layout(matrix(c(rep(1:plotnumber,each = 3), (plotnumber+1):(plotnumber+3))))

    par(mar=c(0,5,0,0),oma = c(1, 1, 1, 1),mgp = c(1,1,0), las=1, bty = "n") # Sets margins per graph and outside margins per grouped set (order is bottom, left, top, right)


      data.chrom <- lapply(out.DSBList, function(x){
      x <- x[x$Chr == chroms[i],]
      return(x)
    })


    data.chrom <- lapply(data.chrom, function(x){
      x$log2fold[x$log2fold == -Inf] <- 0
      x$log2fold[x$log2fold == Inf] <- 0
      loessMod10 <- loess(control=loess.control(surface="direct"),formula= as.numeric(log2fold) ~ as.numeric(Pos), data=x, span=sc/nrow(x))
      xo <- data.frame(Pos = seq(1, max(x$Pos),si)) #change to 100bp
      temp=transform(xo, log2fold = predict(loessMod10, xo))# 10% smoothing span
      temp$Chr=i
      return(temp)
    })


    mycols <- brewer.pal(11, "PRGn")

    if("log2fold" %in% colnames(data.chrom[[1]])) {ylims <- max(mapply(function(x,y) {max(x$log2fold/ (y*bin.width))}, data.chrom, Mreads))}
    else {ylims <- max(mapply(function(x,y) {max(x$Watson/ (y*bin.width))}, data.chrom, Mreads))}
    if(y.mode == "fixed") { ylims <- fixed.y}
    xbreaks <- pretty(data.chrom[[1]]$Pos, n = 10)
    for(y in 1:length(data.chrom)){
     plot(data.chrom[[y]]$Pos, data.chrom[[y]]$log2fold/ (Mreads[y]*bin.width), ylim = c(-ylims,ylims), ylab=wrap.it(paste0(DSBListNames[y]," (log2 fold)"), 40), xlab ="",  axes=FALSE, frame.plot=TRUE, xaxs = "i", yaxs = "i",cex.lab = 0.6, type = "n", col = "dodgerblue3")

        pos.data <- data.chrom[[y]][data.chrom[[y]]$log2fold >0,]
        neg.data <- data.chrom[[y]][data.chrom[[y]]$log2fold <=0,]
        polygon2(pos.data$Pos, pos.data$log2fold, col = mycols[9], border = NA)
        polygon2(neg.data$Pos, neg.data$log2fold, col = mycols[3], border = NA)

      Axis(side = 2)
      if(y == length(data.chrom)) {
        Axis(side=1, at = xbreaks, xaxs = "i", yaxs = "i")
      }
    }
    mtext(text=paste0("Position on chromosome ", chroms[i], " / bp"),side=1,line=0,outer=TRUE,cex=0.8)

    # genetrack
    ################# gene loading
    plot(data.chrom[[y]]$Pos, data.chrom[[y]]$Watson, xaxt="n",yaxt="n",type="n", ylab=paste("Genes"),cex.lab=1.5,font=2,  xaxs = "i", yaxs = "i", ylim=c(-100,120),axes=F) #set up empty plot
    # Axis(side=1, at = pretty(c(xl1, xl2), 10), xaxs = "i", yaxs = "i")
    if(i <= nrow(centro)){segments(centro[i,"chromStart"], 0, centro[i,"chromEnd"], 0, lwd = 3)}
    if(gene.track){
      genes <- features.1[features.1$type == "gene",]
      genesW=subset(genes,genename=="Dubious_ORF" & orientation =="+") #Make a sub-table of ALLElements
      featureplotter(genesW, col = "wheat", lty = 2, strand = "+", av = 75)
      genesW=subset(genes,genename !="Dubious_ORF" & orientation =="+") #Make a sub-table of ALLElements
      featureplotter(genesW, col = "wheat", lty = 1, strand = "+", av = 75)
      genesC=subset(genes,genename=="Dubious_ORF" & orientation =="-") #Make a sub-table of ALLElements
      featureplotter(genesC, col = "thistle", lty = 2, strand = "-", av = 25)
      genesC=subset(genes,genename !="Dubious_ORF" & orientation =="-") #Make a sub-table of ALLElements
      featureplotter(genesC, col = "thistle", lty = 1, strand = "-", av = 25)
    }

    # plot.new()
    # plot(c(1,chromlengths[1,2]), c(-100,100))
    #
    if(ARS.track){
      ARS <- features.1[features.1$type == "ARS_consensus_sequence",]
      if(nrow(ARS) > 0){
        ARS$orientation = "*"
        featureplotter(ARS, col = "tomato3", coln = "sysname", lty = 1, strand = "*", av = 0, tv = 0, pos1 = 1, pos2 = chromlengths[i,"Length"], lab = F)
      }
    }

    dev.off()
  }
}

}


