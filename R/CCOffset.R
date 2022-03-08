#' @title CCOffset
#' @description Generates plots of a pearson offset analysis plots from CCseq fullmaps
#' @param exp.name name of experiment
#' @param path.name path to fullmaps
#' @param os.range range to plot the offset over
#' @param mode specify whether to analyse offsets of Crick relative to Watson ("WC"), or Watson relative to Watson ("WW")
#' @examples CCOffset()
#' @author Will Gittens, George Brown
#' @export
CCOffset<-function(path.name=getwd(),exp.name="NA",os.range = seq(-10, 10, by = 1),mode = "WC") {


  library("doParallel")
  library("data.table")

  mrg_lst <- function(x, y) {
    out <- list()
    for (n in 1:length(x)) {
      out[[n]] <- merge(x[[n]], y[[n]], nomatch = 0L)
    }
    return(out)
  }
  sum_lst <- function(x, y) {
    out <- c()
    # x=temp
    # y="Watson"
    for (n in 1:length(x)) {
      df=x[[n]]
      col=df[, ..y]
      out[n] <- as.numeric(sum(col))
    }
    return(out)
  }
  nrow_lst <- function(x, y) {
    out <- list()
    for (n in 1:length(x)) {
      out[[n]] <- nrow(x[[n]])
    }
    return(out)
  }

  #############WILL'S REAL FUNCTIONS :-)

  # nrow_lst <- function (x)
  # {
  #   out <- c()
  #   for (n in 1:length(x)) {
  #     out <- c(out, nrow(x[[n]]))
  #   }
  #   return(out)
  # }
  #
  # mrg_lst <- function (x, y)
  # {
  #   out <- list()
  #   for (n in 1:length(x)) {
  #     out[[n]] <- merge(x[[n]], y[[n]])
  #   }
  #   return(out)
  # }
  #
  # sum_lst <- function (x, y)
  # {
  #   out <- c()
  #   for (n in 1:length(x)) {
  #     out <- c(out, sum(x[[n]][y]))
  #   }
  #   return(out)
  # }
  #




  setwd(path.name)
  CCList(getwd())
  ###############################################################################################################################################
  # Import histogram FullMap files for each strain in working directory and tally up the total number of Million mapped reads
  parent.directory <- dirname(getwd())

  # DSBList <- DSBList.1HpM
  # nfiles <- length(DSBList)

  # NOTE THE FOLLOWING MAKES EVERYTHING RELATIVE TO MILLION MAPPED READS (MMREADS - AFTER BLACKLIST FILTERING). THIS IS AN IMPORTANT DISTINCTION WHEN CONSIDERING THRESHOLDED DATA, FOR EXAMPLE
  for (i in 1:length(DSBList)){
    DSBList[[i]]=subset(DSBList[[i]],Chr != 17)
    DSBList[[i]]=subset(DSBList[[i]],Chr != 18)
    DSBList[[i]]$Watson <- DSBList[[i]]$Watson/Mreads[i]
    DSBList[[i]]$Crick <- DSBList[[i]]$Crick/Mreads[i]
  }

  ##################
  Offsets <- data.frame(matrix(0, nrow = length(os.range), ncol = 3 * length(DSBList) + 1))
  Offsets[1:length(os.range), 1] <- os.range #create an empty data frame
  names(Offsets) <- c("Pos", paste("Sites_", DSBListNames), paste("HpM_", DSBListNames), paste("Pearson_", DSBListNames)) # give it some column names
  ##################
  if (mode == "WC") {
    osDSBList <- DSBList2 <- DSBList
    osDSBList <- lapply(osDSBList, function(x) {subset(x, x$Crick > 0) }) # drop all rows from osDSBList where Crick is 0
    DSBList2 <-  lapply(DSBList2, function(x) {subset(x, x$Watson > 0) }) # drop all rows from DSBList2 where Watson is 0
    DSBList2 <- lapply(DSBList2, data.table)
    DSBList2 <- lapply(DSBList2, setkey, Chr, Pos)
    osDSBList <-  lapply(osDSBList, function(x) { x$Watson <- NULL; return(x)})
    DSBList2 <-   lapply(DSBList2, function(x) { x$Crick <- NULL; return(x)}) # drop Crick column from DSBList2 and Watson column from osDSBList
    time1 <- Sys.time()
    osDSBList <- lapply(osDSBList, data.table)
    q = 0L
    no_cores <- detectCores() -1
    cl <- makeCluster(no_cores)
    registerDoParallel(cl)
    Offsetsout=Offsets
    Offsetsout = foreach (q = 1:length(os.range),.packages=c("data.table"),.combine = rbind) %dopar% {
      i= os.range[q]
      #for (i in os.range) {
      # loop through the offset range, subtracting this value from the Pos column of osDSBList, to create osDSBList2.
      # q <- q
      # print(q)
      osDSBList2 <- lapply(osDSBList, function(x) {
        x$Pos <- x$Pos - i
        return(x)
      })

      osDSBList2 <- lapply(osDSBList2, setkey, Chr, Pos)
      temp <- mrg_lst(DSBList2, osDSBList2)# merge DSBlist2 and osDSBList2, sample by sample, keeping only rows which are present in both
      # temp <- lapply(temp, function(x) {
      #   x[is.na(x)] <- 0
      #   return(x)
      # })

      Offsets[q,(2:(1+length(DSBList)))]<-nrow_lst(temp)
      Offsets[q,((2+length(DSBList)):(1+2*length(DSBList)))]<-sum_lst(temp,'Watson')+sum_lst(temp,'Crick')
      Offsets[q,((2+2*length(DSBList)):(1+3*length(DSBList)))] <-
        sapply(temp, function(x){
          # library(tidyverse)
          # x$Watson %>% lead()

          cor(x$Watson, x$Crick, method ="pearson")
        })

      print(Sys.time()-time1)
      invisible(Offsets)
    }
    time2 <- Sys.time()
    run.time <- time2 - time1
    print(paste0("Run Time = ", run.time))

  }

  stopCluster(cl)
  Offsets=Offsetsout
  ##################
  if (mode == "WW") {
    DSBList2 <- lapply(DSBList, function(x) { subset(x, x$Watson > 0) })       # drop all rows from DSBList where Watson is 0
    DSBList2 <- lapply(DSBList2, data.table);   DSBList2 <- lapply(DSBList2, setkey, Chr, Pos)
    DSBList2 <-   lapply(DSBList2, function(x) { x$Crick <- NULL; return(x)}) # drop Crick column from DSBList2 and Watson column from osDSBList
    osDSBList <- DSBList2 # copy the data
    osDSBList <- lapply(osDSBList, function(x) { colnames(x)[3] <- "Watson2"; return(x)})
    osDSBList <- lapply(osDSBList, data.table)
    osDSBList <- lapply(osDSBList, setkey, Chr, Pos)
    time1 <- Sys.time()
    q = 0L
    for (i in os.range) {
      # loop through the offset range, subtracting this value from the Pos column of osDSBList, to create osDSBList2.
      q <- q + 1
      if (i == 0)
        next
      osDSBList2 <- lapply(osDSBList, function(x) {
        x$Pos <- x$Pos - i
        return(x)
      })

      osDSBList2 <- lapply(osDSBList2, data.table)
      osDSBList2 <- lapply(osDSBList2, setkey, Chr, Pos)

      temp <-
        mrg_lst(DSBList2, osDSBList2) # merge DSBlist2 and osDSBList2, sample by sample, keeping only rows which are present in both
      temp <- lapply(temp, function(x) {
        x[is.na(x)] <- 0
        return(x)
      })

      Offsets[q,(2:(1+length(DSBList)))]<-nrow_lst(temp)
      Offsets[q,((2+length(DSBList)):(1+2*length(DSBList)))]<-sum_lst(temp,'Watson')+sum_lst(temp,'Watson2')
      Offsets[q,((2+2*length(DSBList)):(1+3*length(DSBList)))] <-
        sapply(temp, function(x){
          cor(x$Watson, x$Watson2, method ="pearson")
        })
    }
    time2 <- Sys.time()
    run.time <- time2 - time1
    print(paste0("Run Time = ", run.time))
  }

  ifelse(!dir.exists(file.path(parent.directory, "Offsets")), dir.create(file.path(parent.directory, "Offsets")), FALSE); setwd(file.path(parent.directory, "Offsets"))
  save(Offsets, file = paste0("OffsetAnalysis_", mode, "_", exp.name, "_", Sys.Date(), ".rds")) #save Offsets as an R object
  write.table(
    Offsets,
    file = paste0("OffsetAnalysis_", mode, "_", exp.name, "_", Sys.Date(), ".txt"),
    sep = "\t",
    row.names = FALSE
  ) #save Offsets as a text file

  #
  #  plot(
  #    Offsets$Pos,
  #    Offsets[,4],
  #    type = "h"
  #  )
  #
  # lines(
  #    Offsets$Pos,
  #    Offsets[,5],
  #    type = "h",
  #    col = "red"
  #  )
  #
  #
  #  plot(
  #    Offsets$Pos,
  #    Offsets$`Reads Cer3H4L2_MJ315_WT_2A_6h`,
  #    type = "h",
  #    ylim = c(0, 1000000)
  #  )
  #
  #  barplot(
  #    Offsets$`Pearson_ top2yftop2mnavg_NOTHSonly`,
  #    main = "Watson-Crick Offset Correlation",
  #    names.arg = Offsets$Pos,
  #    col = "lightblue",
  #    xlab = "Watson-Crick Offset",
  #    ylab = "Pearson Correlation",
  #    ylim = c(0, 1)
  #  )
  #
  #  plot(
  #    temp[[3]]$Watson,
  #    temp[[3]]$Crick,
  #    ylim = c(0, 8),
  #    xlim = c(0, 8),
  #    pch = ".",
  #    col = rgb(1, 0, 0, 0.01)
  #  )
  #  title(main = "")
  #

  library("ggplot2")
  library("ggpubr")
#  Offsets=fread("/Users/georgebrown/Desktop/TC2/results/averages/Fullmaps/Offsets/OffsetAnalysis_WC_NA_2021-06-04.txt")
  Offsets=as.data.frame(Offsets)
  samplesno=ncol(Offsets)
  DSBListNames=rev(DSBListNames)
  for (i in 1:length(DSBListNames)){
    cono=samplesno-(i-1)
    toplot=data.frame(Pos=Offsets$Pos,Pearson=Offsets[,cono])

    p=ggplot(toplot,aes(Pos,Pearson))+geom_col(position = position_dodge(width = 0.9))+theme_bw()+labs(title=paste0(DSBListNames[i]," Watson-Crick Offset Correlation"),x="Watson-Crick Offset",y="Pearson Correlation")+ scale_x_continuous(breaks = os.range,, expand = c(0, 0))+
      theme(panel.grid.minor = element_blank(),panel.grid.major.x=element_blank(),plot.title = element_text(size=10)) +scale_y_continuous(limits = c(0,1), expand = c(0, 0))
    # +xlab("Length of the gene x transcription rate")

    ggsave(plot=p,filename = paste0(DSBListNames[i]," Watson-Crick Offset Correlation.png"), width = 20, height = 20, units = "cm" )
  }


  #  library("ggplot2")
  #  library("ggpubr")
  #  Offsets=fread("/Users/georgebrown/Desktop/TC2/results/averages/Fullmaps/Offsets/OffsetAnalysis_WC_NA_2021-06-03.txt")
  #  Offsets=as.data.frame(Offsets)
  #  samplesno=ncol(Offsets)
  #  for (i in 1:length(DSBListNames)){
  #    cono=samplesno-(i-1)
  #    toplot=data.frame(Pos=Offsets$Pos,Pearson=Offsets[,cono])
  #  toplotrev = toplot
  #  toplot$Strand="Watson"
  #  toplotrev$Strand="Crick"
  #  toplotrev$Pearson=rev(toplot$Pearson)
  #  new=rbind(toplot,toplotrev)
  #  group.colors=c(Watson="red",Crick="steelblue")
  #  p=ggplot(new,aes(Pos,Pearson,fill=Strand))+geom_bar(stat="identity",alpha = 0.5,position =position_identity())+theme_bw()+ylim(0,0.8)+labs(title=paste0(DSBListNames[i]," Watson-Crick Offset Correlation"),x="Offset",y="Pearson Correlation")+
  #  scale_fill_manual(values=group.colors)
  # +xlab("Length of the gene x transcription rate")

  #  ggexport(p,filename = paste0(DSBListNames[i]," Watson-Crick Offset Correlation.pdf"))
  #  }



  setwd(parent.directory)

}
