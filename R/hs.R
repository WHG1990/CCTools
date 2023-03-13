#' @title HS.totals

#' @description Generates a hotspot table from a full map file or rds
#' set the working directory to the location of the full map files
#' @param word.d Working directory where the fullmaps are.
#' @param geneTable location of the AllElementsDUB file
#' @param outdir output directory
#' @param use.BG whether to use BGreads or not
#' @param exp.name experiment name
#' @param readfilename name of the read pattern
#' @examples HS.totals(work.d=getwd(),outdir="/Users/georgebrown/Desktop/foldhsdev"
#' ,Pan="/Users/georgebrown/Dropbox/Annotations/Cer3H4L2/Pan.Hotspots.IGR.SacCer3_H4L2_2016.08.10a.txt")
#' @export
#' @import stringr
#' @import doParallel
#' @import data.table
#' @import e1071
HS.totals <- function(work.d=getwd(),outdir=getwd(),Pan,use.BG=T,readfilename="FullMap.",mode="txt",exp.name="NA",use.mreads=T,
                      geneTable="/Users/georgebrown/Dropbox/Annotations/AllElementsDUB_Cer3H4L2_Brar_2016.08.16.txt",neeman=F,WCequal=F,extend=300){
  setwd(work.d)
  ###############################################################################################################################################
  packages.check("e1071") # This pacakge permits smoothing functions (used later)
  packages.check("stringr")
  options(scipen=999) #Suppresses scientific notation appearing in plots/graphs etc
  packages.check("doParallel")
  packages.check("plyr")
  packages.check("data.table")
  ###############################################################################################################################################
  # Import histogram FullMap files for each strain in working directory and tally up the total number of Million mapped reads
  Mreads=NULL; DSBList=list();dflistNames=NULL
  # if (mode=="txt"){
  #   #Read in all tables with string "Full.Map."
  #   files = list.files(pattern=readfilename) # import files names with "FullMap." string into variable "files"
  #   dflistNames = substr(files, 9, nchar(files)-6) # Shorten filename by 8 characters from beginning and 6 characters form end (i.e. remove "FullMap." and "_c.txt")
  #   nfiles = length(files)
  # }else{
  #   files <- list.files(pattern = readfilename)
  #   file.to.load <- 1
  #   load(files[file.to.load])
  #   dflistNames=DSBListNames
  #
  #   nfiles = length(DSBList)
  # }
  # Count number of files
  # Load required files:

  Pan = fread(Pan)

  if(neeman==F){names(Pan) <- c("Chr","Start","End","Length","PanHits","Feature_name","Name","Midpoint", "Type","Direction","IGR","IGR.start","IGR.end")
  # Pan <- Pan[c("Chr","Start","End","Length","Midpoint","PanHits","Feature_name","Name","Type","Direction","IGR","IGR.start","IGR.end")] # reorder
  Pan <- Pan[,c(1,2,3,4,8,5,6,7,9,10,11,12,13)]
  }
  files = list.files(pattern=readfilename)
  if (use.BG==TRUE){
    bgreads(DUB=geneTable,filenames=readfilename)
    BG = fread("bgreads.txt")
    BGmean=BG$MeanCore
  }else{BGmean=rep(0,length(files))}
  hs=list() # list of tables containing info on the hotspots for each gene for each strain
  no_cores <- detectCores() -1
  cl <- makeCluster(no_cores)
  registerDoParallel(cl)
  k=1
  #writeLines(c(""), "log.txt")
  hs=NULL
  r=NULL
  DSBList=list()
  #step through sequentially each dataframe/strain
  Mreads=NULL
  # DSB map files now loaded in within each loop/instance
  # if (include.rDNA==T){
  #   startrDNA=451000
  #   endrDNA=469000
  #   lengthrDNA=endrDNA-startrDNA
  #   new_row<-c(12,startrDNA,endrDNA,lengthrDNA,0,"rDNA","rDNA",startrDNA+lengthrDNA)
  #   Pan <- rbind(Pan, new_row)
  # }
  if (mode=="txt"){
    for (k in 1:length(files)){
      DSBList[[k]] = fread(files[k], sep = "\t", header=TRUE)
      dflistNames = substr(files, 9, nchar(files)-4)}
  }else{load(files[file.to.load])
    dflistNames=DSBListNames}
  hs=list()
  # for (k in 1:length(DSBList)){
  foreach (k = 1:length(files),.packages=c("data.table")) %dopar% {
      extendwatsonname=paste0("WatsonHpM",extend)
      extendcrickname=paste0("CrickHpM",extend)
      extendtotalname=paste0("TotalHpM",extend)
      extendbgname=paste0("BGHpM",extend)
      extendtotalbgname=paste0("TotalBGHpM",extend)
      extendnormname=paste0("NormHpM",extend)
      extendtotalnormname=paste0("TotalNormHpM",extend)
      extendnormchrname=paste0("NormHpMChr",extend)
      extendtotalnormchrname=paste0("TotalNormHpMChr",extend)



    input=DSBList[[k]]
    if (use.mreads == T){Mreads=sum(input$Watson+input$Crick)/1000000}else{Mreads=1} # Calculate Million reads per sample for conveting to HpM
    setkey(input, Chr)




    hs1l=list()
    for (j in 1:16){
      hs1=NULL
      hs1=data.frame(subset(Pan, Chr==j)) # subsetting the datatables first by chromosome massively speeds the script up!

      # temp1=subset(strainDSBList, Chr==j) # subsetting the datatables first by chromosome massively speeds the script up! I am sure this is most important for the large DSBList tables!
      temp1 <- input[J(j),]
      setkey(temp1, Pos)

      # hs1[,c("WatsonHpM","CrickHpM","TotalHpM","BGHpM","Total-BGHpM","WatsonHpM300","CrickHpM300","TotalHpM300","BGHpM300","Total-BGHpM300","NormHpM","NormHpM300")]=1 #Fill in missing columns before rbind call


      # hs[[k]]=rbind(hs[[k]],hs1)
      # cat("\r", "Job", k, dflistNames[k], "Chromosome", j, "Hotspot", row,"of", nrow(hotspots)); flush.console() # Keep track of progress. cat "\r" overprints to same line of console

      for (i in 1:nrow(hs1)) {
        # temp=temp1[Pos>=(hs1[i,"Start"])& Pos<=(hs1[i,"End"])]

        temp <- temp1[J((as.numeric(hs1[i,"Start"])):(as.numeric(hs1[i,"End"]))), nomatch=0L]
        # temp=subset(temp1, Pos>=(hs1[i,"Start"]) & Pos<=(hs1[i,"End"]))
        # temp=subset(temp1, Pos>=(hs1[i,"Start"]) & Pos<=(hs1[i,"End"])) # Create temp vector with DSB hits across each hotspot in the table

        hs1[i,"WatsonHpM"]=sum(temp$Watson)/Mreads # Calculate sum of hits within this region/Mreads[k]

        hs1[i,"CrickHpM"]=sum(temp$Crick)/Mreads # Calculate sum of hits within this region/Mreads[k]

        if (WCequal == T){hs1[i,"WatsonHpM"] <- sum(pmax(temp$Watson, temp$Crick))/Mreads;hs1[i,"CrickHpM"] <- sum(pmax(temp$Watson, temp$Crick))/Mreads}


        hs1[i,"log2WatsonCrick"]=abs(log2((sum(temp$Crick)+1)/(sum(temp$Watson)+1)))
        temp <- temp1[J((as.numeric(hs1[i,"Start"])):(as.numeric(hs1[i,"Midpoint"]))), nomatch=0L]
        left =sum(temp$Watson,na.rm = T)+0.001
        print(left)
        temp <- temp1[J((as.numeric(hs1[i,"Midpoint"])):(as.numeric(hs1[i,"End"]))), nomatch=0L]
        right =sum(temp$Crick,na.rm = T)+0.001
        hs1[i,"log2leftright"] = abs(log2((right+1)/(left+1)))
        print(j)
        temp <- temp1[J((as.numeric(hs1[i,"Start"])-extend):(as.numeric(hs1[i,"End"])+extend)), nomatch=0L]

        # temp=subset(temp1, Pos>=((hs1[i,"Start"])-extend) & Pos<=((hs1[i,"End"])+extend)) # Create temp vector with DSB hits across hotspot in table +/-300bp
        hs1[i,extendwatsonname]=sum(temp$Watson)/Mreads # Calculate sum of hits within this region/Mreads[k]
        hs1[i,extendcrickname]=sum(temp$Crick)/Mreads# Calculate sum of hits within this region/Mreads[k]
        #
        # row=row+1 # Increment counter

      }



      # hs1$WatsonHpM=hs1$WatsonHpM/Mreads # Calculate sum of hits within this region/Mreads[k]
      # hs1$CrickHpM=hs1$CrickHpM/Mreads # Calculate sum of hits within this region/Mreads[k]
      # hs1$WatsonHpM300=hs1$WatsonHpM300/Mreads # Calculate sum of hits within this region/Mreads[k]
      # hs1$CrickHpM300=hs1$CrickHpM300/Mreads # Calculate sum of hits within this region/Mreads[k]

      hs1$TotalHpM=hs1$WatsonHpM+hs1$CrickHpM

      hs1[,extendtotalname]=hs1[,extendwatsonname]+hs1[,extendcrickname]

      hs1$BGHpM=BGmean[k]*hs1$Length
      hs1[,extendbgname]=BGmean[k]*(hs1$Length+(extend*2))

      hs1$TotalBGHpM=hs1$TotalHpM-hs1$BGHpM
      hs1[,extendtotalbgname]=hs1[,extendtotalname]-hs1[,extendbgname]

      hs1$NormHpM=hs1$TotalHpM-hs1$BGHpM
      hs1$NormHpM[hs1$NormHpM < 0] = 0 # Convert all -ve values to zero


      FinalSum=sum(hs1$TotalHpM)
      hs1$TotalNormHpMChr=hs1$TotalHpM/FinalSum*1000000

      FinalSum=sum(hs1$NormHpM)
      hs1$NormHpMChr=hs1$NormHpM/FinalSum*1000000

      hs1[,extendnormname]=hs1[,extendtotalname]-hs1[,extendbgname]
      hs1[,extendnormname][hs1[,extendnormname] < 0] = 0 # Convert all -ve values to zero

      FinalSum=sum(hs1[,extendtotalname])
      hs1[,extendtotalnormchrname]=hs1[,extendtotalname]/FinalSum*1000000

      FinalSum=sum(hs1[,extendnormname])
      hs1[,extendnormchrname]=hs1[,extendnormname]/FinalSum*1000000
      hs1l[[j]]=hs1
    }
    hs[[k]]=do.call(rbind, hs1l)
    # View(hs[[2]])
    FinalSum=sum(hs[[k]]$NormHpM)
    hs[[k]]$NormHpM=hs[[k]]$NormHpM/FinalSum*1000000
    FinalSum=sum(hs[[k]]$TotalHpM)
    hs[[k]]$TotalNormHpM=hs[[k]]$TotalHpM/FinalSum*1000000
    FinalSum=sum(hs[[k]][,extendnormname])
    hs[[k]][,extendnormname]=hs[[k]][,extendnormname]/FinalSum*1000000
    FinalSum=sum(hs[[k]][,extendtotalname])
    hs[[k]][,extendtotalnormname]=hs[[k]][,extendtotalname]/FinalSum*1000000


    extendnormmitoname=paste0("NormHpMmitocal",extend)
    extendnormmitochrname=paste0("NormHpMChrmitocal",extend)
    extendnormrdnaname=paste0("NormHpMrDNAcal",extend)
    extendnormrdnachrname=paste0("NormHpMChrrDNAcal",extend)


    rORo=rDNACCalibrate(list(input),returnOR=T)
    mORo=mitoCCalibrate2(list(input),returnOR=T)
    hs[[k]]$NormHpMrDNAcal =hs[[k]]$NormHpM*rORo

    hs[[k]]$NormHpMmitocal =hs[[k]]$NormHpM*mORo


    hs[[k]]$NormHpMChrrDNAcal =hs[[k]]$NormHpM*rORo

    hs[[k]]$NormHpMChrmitocal =hs[[k]]$NormHpM*mORo


    hs[[k]][,extendnormrdnaname] =hs[[k]][,extendnormname]*rORo

    hs[[k]][,extendnormmitoname] =hs[[k]][,extendnormname]*mORo


    hs[[k]][,extendnormrdnachrname] =hs[[k]][,extendnormchrname]*rORo

    hs[[k]][,extendnormmitochrname] =hs[[k]][,extendnormchrname]*mORo

    extendtotalmitoname=paste0("TotalHpMmitocal",extend)
    extendtotalrdnaname=paste0("TotalHpMrDNAcal",extend)

    hs[[k]]$TotalHpMrDNAcal =hs[[k]]$TotalHpM*rORo

    hs[[k]]$TotalHpMmitocal =hs[[k]]$TotalHpM*mORo

    hs[[k]][,extendtotalrdnaname] =hs[[k]][,extendtotalname]*rORo

    hs[[k]][,extendtotalmitoname] =hs[[k]][,extendtotalname]*mORo

    # hs[[k]][,extendtotalnormname]
    hs[[k]][,paste0(extendtotalnormname,"rDNAcal")]=hs[[k]][,extendtotalnormname]*rORo
    hs[[k]][,paste0(extendtotalnormname,"mitocal")]=hs[[k]][,extendtotalnormname]*mORo

    hs[[k]]$mOr=mORo

    hs[[k]]$rOr=rORo



    #hs[[k]][1:nrow(hs[[k]]),"Total-BG"]=hs[[k]][1:nrow(hs[[k]]),"Watson"]+hs[[k]][1:nrow(hs[[k]]),"Crick"]-hs[[k]][1:nrow(hs[[k]]),"BGhits"]
    #hs[[k]][1:nrow(hs[[k]]),"Total300-BG"]=hs[[k]][1:nrow(hs[[k]]),"Watson300"]+hs[[k]][1:nrow(hs[[k]]),"Crick300"]-hs[[k]][1:nrow(hs[[k]]),"BGhits300"]

    #hs[[k]][1:nrow(hs[[k]]),"Hit.Increase"]=hs[[k]][1:nrow(hs[[k]]),"Total300-BG"]-hs[[k]][1:nrow(hs[[k]]),"Total-BG"]
    #hs[[k]][1:nrow(hs[[k]]),"Fold.Increase"]=hs[[k]][1:nrow(hs[[k]]),"Total300-BG"]/hs[[k]][1:nrow(hs[[k]]),"Total-BG"]

    #hs[[k]][1:nrow(hs[[k]]),"WC.ratio"]=hs[[k]][1:nrow(hs[[k]]),"Watson"]/hs[[k]][1:nrow(hs[[k]]),"Crick"]
    #hs[[k]][1:nrow(hs[[k]]),"WC.ratio300"]=hs[[k]][1:nrow(hs[[k]]),"Watson300"]/hs[[k]][1:nrow(hs[[k]]),"Crick300"]

    # cat("\r", "Job", k, "COMPLETED", dflistNames[k], "Chromosome", j, "Hotspot", row-1, "of", nrow(hotspots)); flush.console()
    # print(ncol(hs[[k]]))
    # hs[[k]][1:ncol(hs[[k]])]=hs[[k]][1:ncol(hs[[k]])] # For unknown reasons this code is ESSENTIAL to get the script to populate hs[[k]] with anything. Otherwise it returns "NULL"
    # hs[[k]][12:ncol(hs[[k]])]=round(hs[[k]][12:ncol(hs[[k]])], digits=2)

    ######### Write tables to text file ---- Now part of loop ##########

    out = paste(outdir,"/","Hotspot.Table.",dflistNames[k],".txt", sep="")
    write.table(hs[[k]], out, col.names =TRUE, row.names = FALSE, quote = FALSE, sep="\t", append=F)

  }

  # save(hs,file=paste0(outdir,"/Hotspots_",exp.name,".rds"))

  stopCluster(cl)
  # return(hs[[1]])
  sink()
}
#' @title CCfoldHS.calc

#' @description Performs log2(x/y) for a column in a hotspot table
#' @param p.dir Parent directory
#' @param s.dir sample directory containing hotspot tables to be x in log2(x/y)
#' @param c.dir sample directory containing hotspot tables to be y in log2(x/y)
#' @param coln name of the column to be used in log2(x/y)
#' @param readfilename name of the read pattern
#' @examples CCfoldHS.calc(c.dir="/Users/georgebrown/Desktop/foldhsdev/oligosae2d",
#' s.dir="/Users/georgebrown/Desktop/foldhsdev/s1",
#' p.dir=getwd(),coln="NormHpM")
#' @export
CCfoldHS.calc<-function(p.dir=getwd(),s.dir,c.dir,coln="NormHpMChr",readfilename="Hotspot.Table",constant=0){


  setwd(c.dir);controls <- list.files(pattern = readfilename);print(controls)

  setwd(s.dir);files <- list.files(pattern = readfilename);print(files)
  setwd(p.dir)



  hs.renamer<-function(input){
    output=strsplit(input, "Hotspot.Table.")[[1]][2]
    output = strsplit(output, "[, _ -]+")[[1]][1]
    output = strsplit(output, ".txt")[[1]][1]
    print(output)
    return(output)
  }
  packages.check("data.table")
  for (c in 1:length(controls)){
    for (i in 1:length(files)){
      setwd(s.dir);sv=fread(files[i])
      setwd(c.dir);cv=fread(controls[c])
      if (files[i]!=controls[c]){
        sv$log2fold=log2((sv[[coln]]+constant)/(cv[[coln]]+constant))
        fg=data.frame(Chr=sv$Chr,Midpoint=sv$Midpoint,log2fold=sv$log2fold,minname=hs.renamer(controls[c]),plusname=hs.renamer(files[i]))
        print(paste0(hs.renamer(files[i]),"/",hs.renamer(controls[c])))
        fg[sapply(fg, is.infinite)] <- NA
        setwd(p.dir);write.table(fg,file=paste0(hs.renamer(files[i]),"-",hs.renamer(controls[c]),"_",coln,".txt"),quote = FALSE, sep="\t",row.names = FALSE)
      }
    }
  }
}
#' @title CCfoldHS.smoother

#' @description Smooths a log2(x/y) values from a hotspot column using loess
#' @param input filename for a hotspot log2(x/y) table
#' @param sc smoothing constant
#' @param si smooth interval default 100
#' @examples packages.check("doParallel")
#' cl <- makeCluster(8)
#' registerDoParallel(cl)
#' foreach (i = 1:length(files),.packages=c("nealeLabFunctions")) %dopar% {
#' CCfoldHS.smoother(files[i])
#'}stopCluster(cl)
#' @export
CCfoldHS.smoother<-function(input,sc=80,si=100,txtout=F){
  chrSize=c(230218,813184,320870,1531933,576874,270161,1090940,562643,439888,745751,666816,1078177,924431,784333,1091291,948066)
  outputdf=fread(input)
 outputdf<-na.omit(outputdf)
  chrpl=list()
  for (c in 1:16){
    print(c)
    chrp=subset(outputdf, Chr ==c)

    chrp$log2fold[chrp$log2fold == -Inf] <- 0
    chrp$log2fold[chrp$log2fold == Inf] <- 0
    loessMod10 <- loess(control=loess.control(surface="direct"),formula= as.numeric(log2fold) ~ as.numeric(Midpoint), data=chrp, span=sc/nrow(chrp))
    chrpo <- data.frame(Midpoint = seq(1,chrSize[c],si)) #change to 100bp
    chrpl[[c]]=transform(chrpo, log2fold = predict(loessMod10, chrpo))# 10% smoothing span
    chrpl[[c]]$Chr=c

  }
  out.smooth<- do.call("rbind",  chrpl)
  out.smooth$log2fold=round(out.smooth$log2fold,3)
  dir.create("smoothed")

  # write.table(jk,file=paste0("smoothed/smooth_",input), sep = "\t", row.names = F, quote = F)
  save(out.smooth,file=paste0("smoothed/smooth_",substr(input, 1, nchar(input)-4) ,"_sc",sc,".Rbin"), compress=T)
  if(txtout==T){write.table(out.smooth, file = paste0("smoothed/smooth_",substr(input, 1, nchar(input)-4) ,"_sc",sc,".txt"),  sep = "\t", row.names = F, quote = F)}
}
#' @title CCfoldHS.binning

#' @description Bins smoothed hotspot log2(x/y) data
#' @param files a vector of filenames for a smoothed hotspot log2(x/y) table
#' @param i position in the files vector
#' @param binwidth.v A used of the binning binwidths in Kb
#' @param work.d working directory
#' @param si smooth interval default 100
#' @param usesmooth Whether to use smoothed HS data or not. When TRUE plotmode is disabled
#' @param plotmode Whether to generate a plot
#' @param ARSline whether to add ARS positions to the plot
#' @param control.col What colour do you want the control
#' @param mutant.col What colour do you want the mutant
#' @examples files <- list.files(pattern = "m")
#' packages.check("doParallel")
#' cl <- makeCluster(7)
#' registerDoParallel(cl)
#' packages.check("GenomicRanges",BiocM=T)
#' packages.check("ggplot2")
#' packages.check("ggpubr")
#' packages.check("GenomicFeatures",BiocM=T)
#' foreach (i = 1:length(files),.packages=c("nealeLabFunctions","GenomicRanges","ggplot2","ggpubr")) %dopar% {
#'  CCfoldHS.binning(files,i)
#'  }stopCluster(cl)
#' @export
CCfoldHS.binning <-function(files,i=1,binwidth.v=c(10,25,50,100),work.d=getwd(),si=100,plotmode=T,ARSline=F,ylims=2,control.col="lightblue",mutant.col="red",usesmooth=T) {
  if (usesmooth==T){HS <- loadRData(files[i])}else{HS=fread(files[i]);plotmode=F}


      sfiles = strsplit(files[i], ".txt")[[1]][1]
      if (usesmooth==T){sfiles = strsplit(sfiles, "[, _ ]+")[[1]][2]}
  s1=strsplit(sfiles, "[, _ -]+")[[1]][1]
  s2=strsplit(sfiles, "[, _ -]+")[[1]][2]
  for (v in 1:length(binwidth.v)){
    binwidth=binwidth.v[v]
    # Mreads[k]<-sum(temp$Watson + temp$Crick)/ 100000

    convert.dsblist <- function(input){
      output <- data.frame("Chr" = seqnames(input), "start" = start(input), "end" = end(input), "log2fold" = input$log2fold)    #convert back to DSBList format do num of reads in each stacked bar chart.
      return (output)
    }

    # genome <- BSgenome.Scerevisiae.UCSC.sacCer3
    # Chr=c("chrI","chrII","chrIII","chrIV","chrV","chrVI","chrVII","chrVIII","chrIX","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI")
    # Chr=as.character(c(1:16))
    chrSize=c("chrI"=230138,"chrII"=813197,"chrIII"=316643,"chrIV"=1531994,"chrV"=439575,"chrVI"=270161,"chrVII"=1090940,"chrVIII"=562643,"chrIX"=439888,"chrX"=745751,"chrXI"=666816,"chrXII"=1078177,"chrXIII"=924431,"chrXIV"=784333,"chrXV"=1091291,"chrXVI"=948066)
    chrSizeChr=c("1"=230218,"2"=813184,"3"=320870,"4"=1531933,"5"=576874,"6"=270161,"7"=1090940,"8"=562643,"9"=439888,"10"=745751,"11"=666816,"12"=1078177,"13"=924431,"14"=784333,"15"=1091291,"16"=948066)

    chrSize=c(230218,813184,320870,1531933,576874,270161,1090940,562643,439888,745751,666816,1078177,924431,784333,1091291,948066)
    hj=list()
    binwidthbp=binwidth*1000
    for (c in 1:16){
      hj[[c]]=subset(HS,Chr==c)

      # hj[[c]]=tail(hj[[c]],-round(chrSize[c]%%binwidthbp/(si*2)))

      # hj[[c]]=head(hj[[c]],-round(chrSize[c]%%binwidthbp/(si*2)))
       chrSizeChr[c]=chrSize[c] %/% binwidthbp*binwidthbp
    }
    print(chrSizeChr[c])
    cd=NULL
    cd=do.call(rbind, hj)
    cd$Chr=as.character(cd$Chr)
    if (usesmooth==T){cd$log2fold=cd$log2fold*si}else{cd$log2fold=cd$log2fold*(sum(chrSize)/nrow(cd))}
    cd$log2fold[cd$log2fold == -Inf] <- 0
    cd$log2fold[cd$log2fold == Inf] <- 0
    cd[is.na(cd)] <- 0
    CCseq_gr <- makeGRangesFromDataFrame(cd,                           # generate a GRange for your CCseq data
                                         keep.extra.columns=T,
                                         ignore.strand=T,
                                         seqinfo=NULL,
                                         seqnames.field=c("Chr"),
                                         start.field="Midpoint",
                                         end.field=c("Midpoint"),
                                         starts.in.df.are.0based=FALSE)

    cd=NULL
    gr.data.cov <- GenomicRanges::coverage(CCseq_gr,weight = "log2fold")
    gr.windows=tileGenome(chrSizeChr,tilewidth=binwidth*1000,cut.last.tile.in.chrom=T)
    # ranges(gr.windows)=IRanges(start=230218%%100000/2:230218,width=99999)
    # gr <- GRanges(seqnames = c(1:2), ranges = IRanges(start = 1:200000, width = 1000))
    # startx=seq(1,100001,100000)
    # x <- Seqinfo(seqnames=as.character(c(1:16)),
    # seqlengths=chrSize,
    # genome="toy")
    # gr1 <- GRanges("3:15-25", seqinfo=x)
    # gr1 <- GRanges("3:15-25", seqinfo=x)
    # gr <- GRanges(seqnames = 1, ranges = IRanges(start = 1:20000, Endwidth = 2000))
    # gr <- GRanges(seqnames = c("chr1", "chrU345"),
    # ranges = IRanges(start = startx, end = startx+100000))
    # gr.windows=tileGenome(seqinfo(genome),tilewidth=binwidth*1000,cut.last.tile.in.chrom=T)
    # gr.windows=unlist(gr.windows)
    seqlevels(gr.windows, pruning.mode="coarse") <- names(gr.data.cov)

    gr.data.binnedAvg <- binnedAverage(gr.windows, gr.data.cov, "log2fold")
    df=convert.dsblist(gr.data.binnedAvg)
    gr.data.cov=NULL
    CCseq_gr=NULL
    gc()
    print(df[1,4])

    df$mid=(df$end-df$start)/2+df$start
    temp=c()
    temp= ifelse(df$log2fold >0,paste0(s1),paste0(s2))
    temp2= ifelse(df$log2fold ==0,NA,df$log2fold)
    temp3= ifelse(is.na(temp2),temp2,temp)
    df$Library=temp3
    # df$Chr=as.numeric(df$Chr)
    dir.create("binData");dir.create(paste0("binData/",binwidth,"Kb"))
    dir.create("binPlots");dir.create(paste0("binPlots/",binwidth,"Kb"))
    write.table(df,file=paste0("binData/",binwidth,"Kb/",binwidth,"Kb_",substr(files[i], 1, nchar(files[i])-5),".txt"), sep = "\t", row.names = F, quote = F)
    if(plotmode==T){
    centro=CCAnnotate(Cer3H4L2_centro)
    CEN=round((centro$chromEnd+centro$chromStart)/2)
    # CEN=c(151523, 238265, 114443, 449766,152045, 148568, 496979,105644, 355687, 436366,440187, 150887, 268090, 628816, 326643, 556015) # CENtromere positions
    c=10
    options(scipen=10000)
    binplot=df
    print(colnames(binplot))
    # dir.create(paste0("binPlots/",binwidth))
    multi.page=list()
    ARS=CCAnnotate(ARS_consensus)

    test <- CCAnnotate(Cer3H4L2_AllElementsDUB)
    levels(test$type)
    ARS <- test[test$type == "ARS_consensus_sequence",]
    parent.directory <- dirname(work.d)
    if (usesmooth==T){setwd(parent.directory)}
    unsmoothed=list.files(pattern=sfiles)
    library(data.table)
    df=fread(unsmoothed[1])
    setwd(work.d)
    # test <- test[test$chr == 11,]
    for (c in 1:16){
      binchr=subset(binplot,Chr==c)
      print(colnames(binchr))
      toplot=subset(df,Chr==c)
      smooth=subset(HS,Chr==c)
      ARSchr=subset(ARS,chr==c)
      arsmidpoint=as.vector((ARSchr$stop-ARSchr$start)+ARSchr$start)
      ARSv=as.vector(ARSchr$Pos)
      ARSv=arsmidpoint

      # toplot=head(toplot,-1)
      # toplot$mid=toplot$mid/1000
      # toplot$start=toplot$start/1000
      # toplot$end=toplot$end/1000
      if (ARSline==F){ARSv=-1000}
      toplot$Midpoint=toplot$Midpoint/1000
      smooth$Midpoint=smooth$Midpoint/1000
      binchr$mid=binchr$mid/1000
      multi.page[[c]]=ggplot(smooth,aes(Midpoint,log2fold))+
        geom_line()+
        geom_bar(data=binchr,aes(mid,log2fold,fill=Library),stat="identity",width=binwidth,alpha=0.3)+
         geom_vline(xintercept = ARSv/1000)+

        annotate(geom = "label",x=CEN[c]/1000,y=-3, label = "Centromere", show.legend = FALSE,alpha=0.2)+
        # annotate(geom = "label",x=ARSv/1000,y=2, label = "ARS", show.legend = FALSE,alpha=0.2)+

         geom_histogram(stat="identity",data=toplot, aes(Midpoint,log2fold),width=0.5)+

        # geom_line(data=smooth, aes(Midpoint,log2fold))+
        # geom_vline(xintercept = ARSv)+ geom_line(HS[[v]])
        # geom_vline(xintercept = CEN[c])+
        theme_bw(base_size = 18) +
        # +theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
        ylim(-ylims,ylims)+xlim(0,max(smooth$Midpoint))+

        labs(x = paste0("Chromosome ", c, " (Kb)"),y=paste0("Fold log2(",strainTimeSplit(s1),"/",strainTimeSplit(s2),")")) + scale_fill_manual(values = c(control.col, mutant.col))

      # +scale_x_continuous(breaks=seq(0,max(toplot$end),100))

    }
    ggexport(multi.page, filename = (paste0("binPlots/",binwidth,"Kb/",binwidth,"Kb_", substr(files[i], 1, nchar(files[i])-5) ,".pdf")),width = 14,height = 7)
    # chrSize=c("1"=230218,"2"=813184,"3"=320870,"4"=1531933,"5"=576874,"6"=270161,"7"=1090940,"8"=562643,"9"=439888,"10"=745751,"11"=666816,"12"=1078177,"13"=924431,"14"=784333,"15"=1091291,"16"=948066)
    }
  }
}
#' @title CCfoldHS.binpileup

#' @description Pilesup binned smoother log2(x/y) data generated from hotspot tables
#' @param p.dir parent directory containing the binData folder
#' @param binwidth.v A used of the binning binwidths in Kb
#' @param log2FoldMax max and min log2 values to plot
#' @param centromere should the output be centered on the centromere
#' @param EXT should the output be .png or .pdf
#' @param collow Choose negative colour
#' @param colhigh Choose positive colour
#' @param invert.col should the colours be inverted
#' @param fixmax should log2fold values higher than the max be fixed to the max. If set to false such values are displayed in gray.
#' @param log2foldoffset value to offset log2fold values
#' @examples setwd(paste0(getwd(),"/smoothed"));CCfoldHS.binpileup()
#' @import lemon
#' @export
CCfoldHS.binpileup <-function(p.dir=getwd(),binwidth.v=c(100,50,25,10),log2FoldMax=2,centromere=T,EXT=".png",invert.col=F,fixmax=F,log2foldoffset=0,colhigh="red",collow="navy",usesmooth=T){
  packages.check("ggplot2");packages.check("ggpubr");packages.check("lemon")
  for (i in 1:length(binwidth.v)){
    setwd(paste0(p.dir,"/binData/",binwidth.v[i],"Kb"));d.dir=getwd()
    readfilename=paste0(binwidth.v[i],"Kb")
    files = list.files(pattern=readfilename)

    for (f in 1:length(files)){
      sfiles = strsplit(files[f], ".txt")[[1]][1]
      if (usesmooth==F){ti=1}else{ti=0}
      s1=strsplit(sfiles, "[, _ -]+")[[1]][3-ti]
      s2=strsplit(sfiles, "[, _ -]+")[[1]][4-ti]
      coln=strsplit(sfiles, "[, _ -]+")[[1]][5-ti]
      HS <- fread(paste0(files[f]))
      HS$log2fold <- HS$log2fold - log2foldoffset
      if (fixmax==T){
        HS$log2fold=ifelse(HS$log2fold>=log2FoldMax,log2FoldMax,HS$log2fold)
        HS$log2fold=ifelse(HS$log2fold<=-log2FoldMax,-log2FoldMax,HS$log2fold)

      }

      if (invert.col==T){
        temptg=s1
        s1=s2
        s2=temptg
        HS$log2fold=-HS$log2fold
      }
      binwidth=binwidth.v[i]

      HS$length=HS$end-HS$start
      hsl=list()

      CEN=c(151523, 238265, 114443, 449766,152045, 148568, 496979,105644, 355687, 436366,440187, 150887, 268090, 628816, 326643, 556015) # CENtromere positions
      if (centromere==F){CEN=rep(1,16)}
      for (c in 1:16){
        HSt=subset(HS, Chr ==c)
        e=subset(HSt,end<CEN[c]+binwidth*1000 &start>CEN[c]-binwidth*1000)
        cenr=round(e$mid)
        HSt$cen=round((HSt$mid-cenr)/1000)
        hsl[[c]]=HSt
      }
      HSc=do.call("rbind", hsl)

      max=max(HSc$end)
      HSc$cons=1


      chrSize=c("1"=230218,"2"=813184,"3"=320870,"4"=1531933,"5"=576874,"6"=270161,"7"=1090940,"8"=562643,"9"=439888,"10"=745751,"11"=666816,"12"=1078177,"13"=924431,"14"=784333,"15"=1091291,"16"=948066)
      chrsorder =c(4, 15,7,12,16,13,2,14,10,11,5,8,9,3,6,1)
      options(scipen=10000)

      xlims.max=1200
      xlims.min=-xlims.max
      axis.labal.max=1100
      axis.labal.min=-600
      x.label="Position relative to centromere (kb)"
      titlefront="CEN"
      if (centromere==F){
        xlims.max=1800
        xlims.min=0
        axis.labal.max=1500
        axis.labal.min=0
        x.label="Position relative to the left telomere (kb)"
        titlefront="TEL"

      }

      HSc=na.omit(HSc)
      HSc$Chr <- factor(HSc$Chr,      # Reordering group factor levels
                        levels = rev(chrsorder))
      g=ggplot(HSc,aes(x=cen,y=cons,fill=log2fold))+geom_col(width=binwidth)+facet_grid(margins=F,rows=vars(Chr),switch="y")+xlim(xlims.min,xlims.max)+scale_fill_gradient2(
        low = collow, high = colhigh,limits=c(-log2FoldMax,log2FoldMax),breaks=c(-log2FoldMax,0, log2FoldMax),labels=c(paste0(-log2FoldMax,"  ",strainTimeSplit(s2)),0,paste0(log2FoldMax,"  ",strainTimeSplit(s1))))+
        theme_classic2(base_size = 15)+theme(panel.background=element_rect("lightgrey"),strip.text.y.left = element_text(angle = 0),axis.ticks.y=element_blank(),axis.text.y=element_blank(),axis.line.y=element_blank(),axis.title.y = element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(),plot.title = element_text(size=12))+
        xlab(x.label)+labs(title=paste0(coln," Log2 fold change (",strainTimeSplit(s1),"/",strainTimeSplit(s2),") binned at ",binwidth,"kb"))+
        coord_capped_cart(xlim=c(axis.labal.min,axis.labal.max),bottom='both')



      setwd(p.dir)
      dir.create(paste0("binPlots/",binwidth.v[i],"Kb/",titlefront,"pileup"))
      ggsave(g, filename = (paste0("binPlots/",binwidth.v[i],"Kb/",titlefront,"pileup/",titlefront,"pileup_",binwidth,"Kb_",s1,"-",s2,"_",coln,EXT)),width = 8,height = 6)
      # ggexport(g, filename = (paste0("binPlots/",binwidth.v[i],"Kb/",titlefront,"pileup_",binwidth,"Kb_",s1,"-",s2,"_",coln,".pdf")),width = 8,height = 6)

      setwd(d.dir)
    }
  }
  setwd(p.dir)
}
#' @title strainTimeSplit

#' @description Adds a character in a string between the time and strain
#' @param input a string
#' @param toADD character to add
#' @export
strainTimeSplit <-function(input,timeString="hr",toADD="_"){
  te1 = strsplit(input, timeString)[[1]][1]
  if (nchar(te1)<2){
    print("Please rename input files to strainHr rather than HrStrain")
  }
  if (input!=te1){
    te1=substr(te1,1,nchar(te1)-1)
    output = substr(input, nchar(te1)+1, nchar(input))
    output=paste0(te1,toADD,output)
  }else{output=input}

  return(output)
}
#' @title smooth_plotter

#' @description Plots all the log2fold smoothed datasets on the same plot
#' @param work.d path to smoothed.Rbin files
#' @export
smooth_plotter<-function(work.d=getwd()){
  packages.check("ggpubr")
  packages.check("ggplot2")
  packages.check("stringr")
  save.d=getwd()
  setwd(work.d)
  CEN=c(151523, 238265, 114443, 449766,152045, 148568, 496979,105644, 355687, 436366,440187, 150887, 268090, 628816, 326643, 556015) # CENtromere positions

  options(scipen=10000)

  files=list.files(pattern="Rbin")
  nam=c()
  filemam=nam
  for (i in 1:length(files)){
  HS=loadRData(files[i])
  nam[i]=strsplit(files[i], "[, _ ]+")[[1]][2]
  filemam[i]=nam[i]
  nam[i]=str_replace(nam[i], "-", "/")

  assign(nam[i], HS,inherits=T)
  }
  multi.page=NULL
  fg=AppendMe(nam)

  for (c in 1:16){
  toplot=subset(fg,Chr==c)
  multi.page[[c]]=ggplot(toplot,aes(Midpoint,log2fold,col=Library))+geom_line()+
    ylim(-2,2)+xlim(0,max(toplot$Midpoint))+
    annotate(geom = "label",x=CEN[c],y=-2, label = "Centromere", show.legend = FALSE,alpha=0.2)+
    theme_bw(base_size = 18) +
  labs(x = paste0("Chromosome ", c, " (bp)"))
  }
  samples=paste(filemam, collapse='' )
  ggexport(multi.page, filename = (paste0("log2fold",samples,".pdf")),width = 14,height = 7)
  setwd(save.d)
}
