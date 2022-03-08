#' CCHSplotter
#' @description Plots Hotspot data from a hotspot table
#' @param path.name location of the hotspot tables
#' @param names names of the tracks
#' @param plot.mode "png" or "pdf"
#' @param bg "white" or "none"
#' @param coln column to use in the hotspot table
#' @param smooth whether to add a loess smooth
#' @param sc smooth constant
#' @author George Brown
#' @export
CCHSplotter <-function(path.name=getwd(),names,exp.name="NA",plot.mode = "png",bg="white",coln="NormHpMChr",smooth=F,sc=10)
{
  
  features <- CCAnnotate(Cer3H4L2_AllElementsDUB)
  
  chroms <- 1:16
  
  chromlengths <- CCAnnotate(Cer3H4L2_chromlengths); centro <- CCAnnotate(Cer3H4L2_centro)
  
  
                library("data.table")
                library("magick")
                if(!missing(names)){DSBListNames <- names}
  
                sett.name <- paste0("Hotspot")
  
                DSBList<-NULL
                parent.directory <- dirname(path.name)
                setwd(path.name)
                HSfiles=list.files(pattern="Hotspot.Table")
                DSBListNames<- substr(HSfiles, 15, nchar(HSfiles)-4)
                for (i in 1:length(HSfiles)){
                  DSBList[[i]]=fread(HSfiles[i])
                }
  
  
                binned.DSBList <- lapply(DSBList, data.table)
                ####
  
                #offset if you like
  
  
                # If you want full (not sparse) output, this creates an empty binned genome with which you can merge the binned data.
                setwd(parent.directory)
                dir.create("Hotspot")
  
                # Pairwise correlations. This will write a table containing all pairwise pearson correlation values between samples in the binned.DSBList
  
                  names(binned.DSBList) <- DSBListNames
                  write.table(outer(binned.DSBList, binned.DSBList, FUN = Vectorize(function(x,y) {cor(x[[coln]], y[[coln]], method = "pearson")})), file = paste0(exp.name, "_Total_pearson_matrix.txt"), sep = "\t", col.names = DSBListNames, row.names = DSBListNames)
                  test <- binned.DSBList
                  test <- lapply(test, function(x){
                    x <- x[[coln]]
                    return(x)
                  })
                  #plot dendogram
                  setwd("Hotspot")
                  test <- do.call("cbind", test)
                  hc <- hclust(as.dist(1-cor(test, method="pearson")), method="complete")
                  png(file = paste0(exp.name,"_HS_",coln,"_pearson_dendogram.png"), width=10, height=20,  units = "in", res = 400)
                  plot(hc)
                  dev.off()
                  library(magick)
                  image1 <- image_read(paste0(exp.name,"_HS_",coln,"_pearson_dendogram.png"))
                  image_rotate(image1, 90) %>% image_write(paste0(exp.name,"_HS_",coln,"_pearson_dendogram.png"))
                  plotdimensions="wide"
  
                if(plot.mode != "none") {
                  # optional plotting module below. This will plot all chroms.
                  vec <- c(66.2,46.8,33.1,23.4,16.5,11.7,8.3,5.8,4.1,2.9,2.0,1.5,1, 7.5, 4.6, 3.3, 3.75)
                  vec2 <- c(93.6,66.2,46.8,33.1,23.4,16.5,11.7,8.3,5.8,4.1,2.9,2.0,1.5, 13.3, 8.3, 5.8, 6.65)
                  vec3 <- c("4A0", "2A0", "A0", "A1", 'A2', 'A3', 'A4', 'A5', 'A6', "A7", 'A8', 'A9', 'A10', "wide", "wide2", "wide3", "wide4")
                  papersize <- data.frame(vec3,vec2,vec)
                  plotdimension <- as.numeric(papersize[match(plotdimensions, vec3),2:3])
                  
                  
                  sf=5;if (length(DSBListNames) > 5){sf=length(DSBListNames)}
                  dir.create("ChromosomePlots")
                  setwd("ChromosomePlots")
                  DSBListNames<-names
                  for(i in 1:length(chroms)){
                    features.1 <- features[features$chr == i,]
                    
                    plot.name <- paste0("Chr", i, "_",exp.name,"_", coln, "_", sett.name)
                    
                    
                    # png(file = paste0(plot.name, ".png"), width=plotdimension[1], height=(plotdimension[2]/5)*sf,  units = "in", res = 400)
                    
                    if(plot.mode == "pdf"){pdf(file = paste0(plot.name, ".pdf"), width=plotdimension[1], height=(plotdimension[2]/5)*sf)}
                    if(plot.mode == "png"){png(file = paste0(plot.name, ".png"), width=plotdimension[1], height=(plotdimension[2]/5)*sf,  units = "in", res = 400)}
                    
                    # png(file = paste0(plot.name, ".png"), width=plotdimension[1], height=(plotdimension[2]/5)*sf,  units = "in", res = 600)
                    plotnumber=length(binned.DSBList) # Number from 1 to 5
                    # if(ARS.track) {layout(matrix(c(rep(1:plotnumber,each = 3), (plotnumber+1):(plotnumber+3))))}
                    
                    layout(matrix(c(rep(1:plotnumber,each = 3),(plotnumber+1) )))
                    
                    if(bg == "none"){par(mar=c(0,7,0,0),oma = c(1, 1, 1, 1),mgp = c(3,1,0),  bty = "n", bg = NA)} # Sets margins per graph and outside margins per grouped set (order is bottom, left, top, right)
                    if(bg == "white"){par(mar=c(0,7,0,0),oma = c(1, 1, 1, 1),mgp = c(3,1,0),  bty = "n")} # Sets margins per graph and outside margins per grouped set (order is bottom, left, top, right)
                    
                    data.chrom <<- lapply(binned.DSBList, function(x){
                    x <- x[x$Chr == chroms[i],]
                    return(x)
                    })
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    si=1000
                    xbreaks <- pretty(data.chrom[[1]]$Midpoint, n = 10)
                    chrpl=NULL
                    ylims=1
                    for(y in 1:length(data.chrom)){
                    if (ylims<max(data.chrom[[y]][[coln]])){
                    ylims=max(data.chrom[[y]][[coln]])}
                    print(ylims)}
                    ylims=ylims*1.2
                    for(y in 1:length(data.chrom)){
                    plot(data.chrom[[y]]$Midpoint, data.chrom[[y]][[coln]],  ylim = c(0,ylims),xlim = c(0,chromlengths[i,2]), ylab=wrap.it(paste0(DSBListNames[y]," (HpM)"), nchar(DSBListNames[y])), xlab ="",  axes=FALSE, frame.plot=TRUE, xaxs = "i", yaxs = "i",cex.lab = 0.6, type = "h", col = "dodgerblue3")
                    # Axis(side = 2)
                    if(smooth==T){
                      chrp=na.omit(data.chrom[[y]])
                      chrp[[coln]][chrp[[coln]] == -Inf] <- 0
                      chrp[[coln]][chrp[[coln]] == Inf] <- 0
                      chrp$tosmooth=chrp[[coln]]
                      loessMod10 <- loess(control=loess.control(surface="direct"),formula= as.numeric(tosmooth) ~ as.numeric(Midpoint), data=chrp, span=sc/nrow(chrp))
                      chrpo <- data.frame(Midpoint = seq(1,chromlengths[i,2],si)) #change to 100bp
                      chrpl[[y]]=transform(chrpo, smoothed = predict(loessMod10, chrpo))# 10% smoothing span
                      chrpl[[y]]=subset(chrpl[[y]],Midpoint>=min(chrp$Midpoint)) 
                      chrpl[[y]]=subset(chrpl[[y]],Midpoint<=max(chrp$Midpoint)) 
                      chrpl[[y]]$Chr=i
                      
                    lines(chrpl[[y]]$Midpoint,chrpl[[y]]$smoothed)
                    }
                    # Axis(side = 2)
                    if(y == length(data.chrom)) {
                      Axis(side=1, at = xbreaks, xaxs = "i", yaxs = "i")
                    }
                    }
                    mtext(text=paste0("Position on chromosome ", chroms[i], " / bp"),side=1,line=0,outer=TRUE,cex=0.8)
                    
                    # genetrack
                    ################# gene loading
                    plot(data.chrom[[y]]$Midpoint, data.chrom[[y]][[coln]], xaxt="n",yaxt="n",type="n", ylab=paste(" "),cex.lab=1.5,font=2,  xaxs = "i", yaxs = "i", ylim=c(-100,120),axes=F) #set up empty plot
                    # Axis(side=1, at = pretty(c(xl1, xl2), 10), xaxs = "i", yaxs = "i")
                    if(i <= nrow(centro)){segments(centro[i,"chromStart"], 0, centro[i,"chromEnd"], 0, lwd = 3)}
                    
                    
                    # if(ARS.track){
                    #   ARS <- features.1[features.1$type == "ARS_consensus_sequence",]
                    #   if(nrow(ARS) > 0){
                    #     ARS$orientation = "*"
                    #     featureplotter(ARS, col = "tomato3", coln = "sysname", lty = 1, strand = "*", av = 50, tv = 75, pos1 = 1, pos2 = chromlengths[i,"Length"], lab = F)
                    #   }
                    # }
                    
                    dev.off()
                  }
                }
                setwd(path.name)
                return(chrpl)
}
