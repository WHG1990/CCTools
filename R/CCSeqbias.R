#' CCSeqbias
#' @description generates a sequence bias plot.
#' @param path.name point to a folder containing fullmaps
#' @param minThresholdv A vector containing the min threshold
#' @param maxThreshold A value containing the max
#' @param w Window width in bp
#' @param rDNA include rDNA
#' @param chroms A Vector of chromosomes to process.
#' @param ymin Min Y axis value
#' @param ymax Max Y axis value
#' @param dyad Dyad to plot center lines from (0.5 for Spo11, 1.5 for Top2)
#' @param genome location to genome.fa
#' @param loci this can be set to a specific genename. The sequence bias will be calcuated only for sites within the genebody of this gene. HpM thresholding occurs before this, but you will probably still need to modify min max thresholds.
#' @author George Brown, Matt Neale
#' @export
CCSeqbias <-function(path.name=getwd(),minThresholdv=c(0,1,5),maxThreshold="MAX",w=20,rDNA=F,chroms=c(1:16),ymin=0.05,ymax=0.6,dyad=0.5,genome="Cer3H4L2", loci){
  
  library(stringr)
  library(readr)
  library(data.table)
  library(stringi)
  library(doParallel)
  
  # optional loci subset
  if(!missing(loci)) {
    features <- CCAnnotate(Cer3H4L2_AllElementsDUB)
    genes <- features[features$type == "gene",]
    generegion <- genes[genes$genename == loci,]
    pos1 <- generegion$start
    pos2 <- generegion$stop
    # if(is.na(loci[2])) {loci <- c(loci,"mid")}
    # loci.table <- lociRead(loci = loci, features = features, chrom. = chrom., pos = pos)
    # if (tolower(loci[2])== "m" || loci[2]=="middle"||loci[2]=="mid"){
    #   fplotpos="Mid"
    # }
    # else if (tolower(loci[2])== "s" ||loci[2]== "start"){
    #   fplotpos="Start"
    # }
    # else if (tolower(loci[2])== "stop" ||loci[2]== "end"||loci[2]=="e"){
    #   fplotpos="End"
    # }
    
  }else{loci=0}
  
  
  script="SeqBiasR_v02a"
  parent.directory=dirname(path.name)
  setwd(parent.directory);dir.create("seqBias")
  setwd(path.name)
  CCList(getwd())
  DSBList=DSBList
  DSBListNames=DSBListNames
  # for (i in 1:length(DSBListNames)){
  #
  #   dir.create(DSBListNames[i])
  # }
  
  
  
  #Convert to HpM
  no_cores <- detectCores() -1
  cl <- makeCluster(no_cores)
  registerDoParallel(cl)
  
  
  if (genome == "Cer3H4L2"){
    ref=CCAnnotate(Cer3H4L2_genome)
  }else if(genome == "W303"){
    ref=CCAnnotate(W303_genome)
  }else if (genome == "Cer3H4L2MATa"){
    ref=CCAnnotate(Cer3H4L2MATa_genome)
  }else if (genome == "pombase220208"){
    ref  <- readChar("/Users/wg45/Dropbox/Work/Top2Seq/Genomes_and_Annotations/Pombe/temppombe2/pombase220208.fa", file.info("/Users/wg45/Dropbox/Work/Top2Seq/Genomes_and_Annotations/Pombe/temppombe2/pombase220208.fa")$size)
  }else{
    ref=read_file(genome)
  }
  
  
  
  
  
  
  
  chromsC=vec2str(chroms)
  
  
  foreach (k = 1:length(DSBList),.packages=c("data.table","stringr","stringi","readr")) %dopar% {
    data=DSBList[[k]]
    
    Mreads=sum(data$Watson+data$Crick)/1000000
    data$Watson=data$Watson/Mreads
    data$Crick=data$Crick/Mreads
    ########################################################################################################
    
    #Read in Reference and tidy up
    #ref=read_file("Cer3H4L2_Chr1.txt")
    
    #ref=read_file("W303_MJN.fa")
    #ref=read_file("Cer3H4L2_1-3only.txt")
    
    #Tidy up fasta
    ref=str_replace_all(ref, "[\r\n]" , "") #Strip out newlines
    refA=str_split(ref,">") #Split into chromosomes by ">" character
    refB=NULL
    refC=NULL
    for (n in 1:18){ #Work through each chromosome
      refB[n]=refA[[1]][n+1] #Vector of each chromosome as a string
    }
    
    for (n in 1:9){ #Work through each chromosome
      refC[n]=str_sub(refB[n],2) #Remove leading chromosome number in each string
    }
    
    for (n in 10:18){ #Work through each chromosome
      refC[n]=str_sub(refB[n],3) #Remove leading chromosome number in each string
    }
    
    #########################################################################################################
    #########################################################################################################
    ######   START HERE ONCE DATA READ and REFERNCE MADE  #######################
    #########################################################################################################
    #########################################################################################################
    
    #Set loop of different min thresholds
    
    #Loop
    t2a=maxThreshold
    if (maxThreshold=="MAX"|| 1000000){maxThreshold=1000000;t2a="MAX"}
    for (j in minThresholdv){
      # j=0
      #Analysis Options
      minThreshold=j # Minimum value for position to be sampled NOTE: Use at least 1
      
      # chroms=c(1:16) #All
      #chroms=c(12) #Vector of chromosomes to analyse
      #chroms=c(17) #Vector of chromosomes to analyse
      
      #Code
      # mapS=substr(map, 9, nchar(map)-4) #Shorten output filename
      # mapS=strsplit(map, "FullMap.")[[1]][3]
      # mapS=strsplit(mapS, ".txt")[[1]][1]
      mapS=DSBListNames[k]
      width=2*w+1 #Full width of analysed region
      
      
      #Empty base results vectors
      #Franklin
      Af=NULL
      Cf=NULL
      Gf=NULL
      Tf=NULL
      c1f=NULL
      lociF=NULL
      #Rosalind
      Ar=NULL
      Cr=NULL
      Gr=NULL
      Tr=NULL
      c1r=NULL
      lociR=NULL
      
      # if( rDNA ==F) {data<-data[with(data, (Chr == 12 & Pos >= 451000 & Pos <= 469000)),]}
      
      #Main loop
      #Franklin
      # x=data
      # x=subset(x,Chr==12)
      # x<-x[with(x, Pos <= 451000 | Pos >= 469000),]
      # x<-x[with(x, (Chr == 12 & Pos >= 451000 & Pos <= 469000)),]
      maxpossiblesitesFv=c()
      maxpossiblesitesRv=c()
      for (i in chroms){ #Step through each chromosome
        
        chrX=subset(data,Chr==i) #Subset data by each chromosome
        
        if( rDNA ==F & i == 12) {chrX<-chrX[with(chrX, Pos <= 451000 | Pos >= 469000),]}
        #Franklin=subset(chrX,Watson>0) #First Remove all positions with zero hits
        Franklin=subset(chrX,Watson>0)
        if(loci !=0)  {Franklin <- Franklin[Franklin$Pos > pos1 & Franklin$Pos < pos2,] }
        maxpossiblesitesFv[i]=nrow(Franklin)
        Franklin=subset(Franklin,Watson>minThreshold&Watson<=maxThreshold) #Subset by Rosalind>=minThreshold&<=maxThreshold
        Rosalind=subset(chrX,Crick>0)
        if(loci !=0)  {Rosalind <- Rosalind[Rosalind$Pos > pos1 & Rosalind$Pos < pos2,] }
        maxpossiblesitesRv[i]=nrow(Rosalind)
        Rosalind=subset(Rosalind,Crick>minThreshold&Crick<=maxThreshold) #Subset by Rosalind>=minThreshold&<=maxThreshold
        
        
        xF=Franklin$Pos #Resulting vector of Franklin positions
        # xR=Rosalind$Pos #Resulting vector of Rosalind positions
        a=refC[i] #Current chromosome fasta string
        lociF=c(lociF,str_sub(a,xF-w,xF+w)) #This concatenates a vector (loci) of substrings centred on position (x) ± width (w) for current chromosome (i)
        # lociR=c(lociR,str_sub(a,xR-w,xR+w)) #This concatenates a vector (loci) of substrings centred on position (x) ± width (w) for current chromosome (i)
        x=Rosalind$Pos #Resulting vector of Rosalind positions
        a=refC[i] #Current chromosome fasta string
        lociR=c(lociR,str_sub(a,x-w,x+w)) #This concatenates a vector (loci) of substrings centred on position (x) ± width (w) for current chromosome (i)
      }
      SampledSitesF=length(lociF)
      maxpossiblesitesF=sum(maxpossiblesitesFv)
      
      # for (i in chroms){ #Step through each chromosome
      #   chrX=subset(data,Chr==i) #Subset data by each chromosome
      #   if( rDNA ==F & i == 12) {chrX<-chrX[with(chrX, Pos <= 451000 | Pos >= 469000),]}
      #
      #   #Rosalind=subset(chrX,Crick>0) #Remove positions with zero hits
      #   Rosalind=subset(chrX,Crick>0)
      #   maxpossiblesitesRv[i]=nrow(Rosalind)
      #   Rosalind=subset(Rosalind,Crick>minThreshold&Crick<=maxThreshold) #Subset by Rosalind>=minThreshold&<=maxThreshold
      #   x=Rosalind$Pos #Resulting vector of Rosalind positions
      #   a=refC[i] #Current chromosome fasta string
      #   lociR=c(lociR,str_sub(a,x-w,x+w)) #This concatenates a vector (loci) of substrings centred on position (x) ± width (w) for current chromosome (i)
      # }
      SampledSitesR=length(lociR)
      
      maxpossiblesitesR=sum(maxpossiblesitesRv)
      
      SMF=round(SampledSitesF/maxpossiblesitesF,2)
      SMR=round(SampledSitesR/maxpossiblesitesR,2)
      SMT=round((SampledSitesR+SampledSitesF)/(maxpossiblesitesR+maxpossiblesitesF),2)
      
      #Calculate %base at each Franklin position
      start.total=Sys.time()
      
      for (n in 1:width){ #Step through each position in the sampled width region
        start.time=Sys.time()
        c1f=str_sub(lociF,n,n) #c1f=vector through lociF of characters at position n
        # c1r=stri_sub(lociR,n,n)
        # c1f=lociF
        end.time=Sys.time()
        print(1)
        print(end.time-start.time)
        start.time=Sys.time()
        Af[n]=length(stri_subset(c1f,fixed="A"))/SampledSitesF #Count number of A in c1f; divide by total sites (length of lociF)
        # Ar[width-n+1]=length(stri_subset(c1f,fixed="A"))/SampledSitesR
        end.time=Sys.time()
        print(2)
        print(end.time-start.time)
        start.time=Sys.time()
        Cf[n]=length(stri_subset(c1f,fixed="C"))/SampledSitesF
        # Cr[width-n+1]=length(stri_subset(c1f,fixed="C"))/SampledSitesR
        end.time=Sys.time()
        print(3)
        print(end.time-start.time)
        start.time=Sys.time()
        Gf[n]=length(stri_subset(c1f,fixed="G"))/SampledSitesF
        # Gr[width-n+1]=length(stri_subset(c1f,fixed="G"))/SampledSitesR
        end.time=Sys.time()
        print(4)
        print(end.time-start.time)
        start.time=Sys.time()
        Tf[n]=length(stri_subset(c1f,fixed="T"))/SampledSitesF
        # Tr[width-n+1]=length(stri_subset(c1f,fixed="T"))/SampledSitesR
        end.time=Sys.time()
        print(5)
        print(end.time-start.time)
        start.time=Sys.time()
        cat("\r", "Pos", n); flush.console() # Keep track of progress. cat "\r" overprints to same line of console
      }
      
      #Calculate %base at each Rosalind position. Note these are now reversed, bottom-strand bases
      # for (n in 1:width){
      #   c1r=stri_sub(lociR,n,n)
      #   Ar[width-n+1]=length(stri_subset(c1r,fixed="A"))/SampledSitesR
      #   Cr[width-n+1]=length(stri_subset(c1r,fixed="C"))/SampledSitesR
      #   Gr[width-n+1]=length(stri_subset(c1r,fixed="G"))/SampledSitesR
      #   Tr[width-n+1]=length(stri_subset(c1r,fixed="T"))/SampledSitesR
      #   cat("\r", "Pos_R", n); flush.console() # Keep track of progress. cat "\r" overprints to same line of console
      # }
      
      for (n in 1:width){
        c1r=str_sub(lociR,n,n)
        Ar[width-n+1]=length(stri_subset(c1r, fixed="T"))/SampledSitesR
        Cr[width-n+1]=length(stri_subset(c1r, fixed="G"))/SampledSitesR
        Gr[width-n+1]=length(stri_subset(c1r, fixed="C"))/SampledSitesR
        Tr[width-n+1]=length(stri_subset(c1r, fixed="A"))/SampledSitesR
        cat("\r", "Pos_R", n); flush.console() # Keep track of progress. cat "\r" overprints to same line of console
      }
      
      
      end.total=Sys.time()
      end.total-start.total
      
      
      ##########################################################################################################################################################
      ##########################################################################################################################################################
      #Plotting - START here to replot after modifying ylims etc.
      
      #Plotting OPTIONS:
      
      vert=seq(from =-round(w,-1), to=round(w,-1), by=5)+dyad #Add vertical grid lines centred on dyad
      plotT="o" #PlotType: "l", "p", "o"
      
      ##########################################################################################################################################################
      setwd(parent.directory)
      dir.create(paste0("seqBias/",mapS))
      wd = getwd();
      if(loci ==0) {out = paste(wd,"/","seqBias","/",mapS,"/",script,"_",mapS,"_rDNA=",rDNA,"_±",w,"_T>",minThreshold, "to",t2a,"_HpM_dyad=",dyad,"_C",chromsC,".pdf",sep="")}
      if(loci !=0)  {out = paste(wd,"/","seqBias","/",mapS,"/",script,"_",mapS,"_rDNA=",rDNA,"_±",w,"_T>",minThreshold, "to",t2a,"_HpM_dyad=",dyad,"_",loci, "only.pdf",sep="")}
      ; pdf(file=out, width=9,height=18) # Plot PDF
      layout(matrix(c(1,1,1,2,2,2,3,3,3),9, 1, byrow = T))
      
      #Plotting_F
      if(loci ==0)  {plot(0,0,xlim=c(-w,w), ylim=c(ymin,ymax), type="n", xlab="Position (Zero is the 5'end)", ylab="Base Frequency", main=paste0(mapS," Mreads = ",round(Mreads,2),"\nFranklin_Sites = ",SampledSitesF,"(",SMF*100,"%)"," / Thresh > ",minThreshold, " HpM to ",t2a," HpM / Chroms = ",paste0(chromsC, collapse=" ")))}
      if(loci !=0)  {plot(0,0,xlim=c(-w,w), ylim=c(ymin,ymax), type="n", xlab="Position (Zero is the 5'end)", ylab="Base Frequency", main=paste0(mapS," Mreads = ",round(Mreads,2),"\nFranklin_Sites = ",SampledSitesF,"(",SMF*100,"%)"," / Thresh > ",minThreshold, " HpM to ",t2a," HpM / Gene = ", loci))}
      polygon(c(-w,w,w,-w),c(ymin, ymin, ymax, ymax),col="grey", border=NA)
      abline(v=dyad, col="white", lwd=4)
      abline(v=vert, col="white")
      lines((-w:w),Af[1:width], type=plotT, pch=19, col="yellow",cex=0.01)
      lines((-w:w),Cf[1:width], type=plotT, pch=19, col="green",cex=0.01)
      lines((-w:w),Gf[1:width], type=plotT, pch=19, col="red",cex=0.01)
      lines((-w:w),Tf[1:width], type=plotT, pch=19, col="blue",cex=0.01)
      legend(w*0.85, ymax*0.95, legend=c("A", "T", "C", "G"), col=c("yellow", "blue", "green", "red"), lty=1, cex=1, bg="white")
      
      #Plotting_R
      if(loci ==0)  {plot(0,0,xlim=c(-w,w), ylim=c(ymin,ymax), type="n", xlab="Position (Zero is the 5'end)", ylab="Base Frequency", main=paste0(mapS," Mreads = ",round(Mreads,2),"\nRosalind_Sites = ",SampledSitesF,"(",SMF*100,"%)"," / Thresh > ",minThreshold, " HpM to ",t2a," HpM / Chroms = ",paste0(chromsC, collapse=" ")))}
      if(loci !=0)  {plot(0,0,xlim=c(-w,w), ylim=c(ymin,ymax), type="n", xlab="Position (Zero is the 5'end)", ylab="Base Frequency", main=paste0(mapS," Mreads = ",round(Mreads,2),"\nRosalind_Sites = ",SampledSitesF,"(",SMF*100,"%)"," / Thresh > ",minThreshold, " HpM to ",t2a," HpM / Gene = ", loci))}
      polygon(c(-w,w,w,-w),c(ymin, ymin, ymax, ymax),col="grey", border=NA)
      abline(v=dyad, col="white", lwd=4)
      abline(v=vert, col="white")
      lines((-w:w),Ar[1:width], type=plotT, pch=19, col="yellow",cex=0.01)
      lines((-w:w),Cr[1:width], type=plotT, pch=19, col="green",cex=0.01)
      lines((-w:w),Gr[1:width], type=plotT, pch=19, col="red",cex=0.01)
      lines((-w:w),Tr[1:width], type=plotT, pch=19, col="blue",cex=0.01)
      legend(w*0.85, ymax*0.95, legend=c("A", "T", "C", "G"), col=c("yellow", "blue", "green", "red"), lty=1, cex=1, bg="white")
      
      #Plotting_F+R
      if(loci ==0)  {plot(0,0,xlim=c(-w,w), ylim=c(ymin,ymax), type="n", xlab="Position (Zero is the 5'end)", ylab="Base Frequency", main=paste0(mapS," Mreads = ",round(Mreads,2),"\nCombined_Sites = ",SampledSitesF,"(",SMF*100,"%)"," / Thresh > ",minThreshold, " HpM to ",t2a," HpM / Chroms = ",paste0(chromsC, collapse=" ")))}
      if(loci !=0)  {plot(0,0,xlim=c(-w,w), ylim=c(ymin,ymax), type="n", xlab="Position (Zero is the 5'end)", ylab="Base Frequency", main=paste0(mapS," Mreads = ",round(Mreads,2),"\nCombined_Sites = ",SampledSitesF,"(",SMF*100,"%)"," / Thresh > ",minThreshold, " HpM to ",t2a," HpM / Gene = ", loci))}
      polygon(c(-w,w,w,-w),c(ymin, ymin, ymax, ymax),col="grey", border=NA)
      abline(v=dyad, col="white", lwd=4)
      abline(v=vert, col="white")
      lines((-w:w),(Af[1:width]+Ar[1:width])/2, type=plotT, pch=19, col="yellow",cex=0.01)
      lines((-w:w),(Cf[1:width]+Cr[1:width])/2, type=plotT, pch=19, col="green",cex=0.01)
      lines((-w:w),(Gf[1:width]+Gr[1:width])/2, type=plotT, pch=19, col="red",cex=0.01)
      lines((-w:w),(Tf[1:width]+Tr[1:width])/2, type=plotT, pch=19, col="blue",cex=0.01)
      legend(w*0.85, ymax*0.95, legend=c("A", "T", "C", "G"), col=c("yellow", "blue", "green", "red"), lty=1, cex=1, bg="white")
      
      ##########################################################################################################################################################
      dev.off() #Close connection to PDF writer
      setwd(path.name)
      ##########################################################################################################################################################
    } #Close j loop
  }
}
#' @title chrom.vector.rename
#' @description renames a vector of numbers to a short hand string
#' @examples chrom.vector.rename(c(1:16))
#' @author George Brown
#' @export
vec2str<-function(input){


  input=sort(input)


strin = ""
previousChrome = 0
first=0
chain=F
if (length(input)>1){
for(chrome in input) {

  if (previousChrome == 0) {
    firstChrome = chrome
  }
  if (chrome - previousChrome > 1) {
    # print(chrome)
    # print(previousChrome)
    chain=T
    print(paste0(firstChrome, "-", previousChrome))
    first = first +1
    if (first==1){
    strin=paste0(strin,firstChrome, "-", previousChrome)
    }else{
      if (firstChrome==previousChrome){strin=paste0(strin,"_",firstChrome)}else{
      strin=paste0(strin,"_",firstChrome, "-", previousChrome)}}
    firstChrome = chrome
  } else {
    chain=F
    # strin=paste0(strin,chrome)
  }
  #print(str)
  previousChrome = chrome
}
 if (first==0){
 if (chain ==1){strin=paste0(strin,firstChrome)}else{strin=paste0(firstChrome, "-", previousChrome)}}else{
if (chain ==T){strin=paste0(strin,"_",firstChrome)}else{strin=paste0(strin,"_",firstChrome, "-", previousChrome)}}
}else{
  strin=paste0(input[1])
}
return(strin)
}
