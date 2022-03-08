#' VX
#' @description VariableX smoother which takes a sparse CCseq dataframe or datatable (colnames = c("Chr", "Pos", "Watson", "Crick)), and smooths both the y dimension ("Watson" and "Crick"), and the x dimension ("Pos"). A sliding mean height and mean position of Fsize number of hits is calculated.
#' @param x A sparse CCseq dataframe or datatable (colnames = c("Chr", "Pos", "Watson", "Crick)). E.g the "sae2.0" table generated within mapper.
#' @param Fsize The number of single strand sites (non-zero "Watson" or "Crick" rows) or double-stranded sites ("Total" rows) over which to smooth. This is the number of values that will be averaged, NOT the  distance in bp.
#' @author Matt Neale and Will Gittens
#' @export
VX <- function(x, Fsize = 101, norm.factor = 3){
  f=rep(1,Fsize) #Create a SUM filter. Every valeu in the window of interest will be SUMMED (1 means ADD).
  Fdif=c(1,rep(0,Fsize-2),-1) #Create a filter that is will subtract the first from the last value on the window of interest. This calculates window width
    x$Total <- x$Watson + x$Crick
  b=subset(x,Watson>=1)
  b$Pos <- as.numeric(b$Pos); b$Watson <- as.numeric(b$Watson) # important to avoid integer overflow!
  b<-b[order(b$Pos),]
  b$PosDif=filter(b$Pos,Fdif) #Calculate width of averaging window
  b$WatsonAve=filter(b$Watson,f)/b$PosDif #Divide total Watson signal in window by the window width
  b$PosAve<-b$Watson*b$Pos  #Calculate x-weighted midpoint of averaging window, note that this essentially gives the mean position (rather than median). This is intentional.
  b$PosAve<-filter(b$PosAve,f) #Calculate x-weighted midpoint of averaging window, note that this essentially gives the mean position (rather than median). This is intentional.
  b$PosAve<-b$PosAve/filter(b$Watson,f)

  d=subset(x, Crick>=1)
  d$Pos <- as.numeric(d$Pos); d$Crick <- as.numeric(d$Crick)
  d<-d[order(d$Pos),]
  d$PosDif=filter(d$Pos,Fdif) #Calculate width of averaging window
  d$CrickAve=filter(d$Crick,f)/d$PosDif #Divide total Watson signal in window by the window width
  d$PosAve<-d$Crick*d$Pos   # These 3 lines Calculate midpoint of averaging window, note that this essentially gives the mean position (rather than median). This is intentional.
  d$PosAve<-filter(d$PosAve,f)
  d$PosAve<-d$PosAve/filter(d$Crick,f)

  e=subset(x, Total>=1)
  e$Pos <- as.numeric(e$Pos); e$Total <- as.numeric(e$Total)
  e<-e[order(e$Pos),]
  e$PosDif=filter(e$Pos,Fdif)
  e$Total<-e$Watson+e$Crick#Calculate width of averaging window
  e$TotalAve=filter(e$Total,f)/e$PosDif #Divide total Total signal in window by the window width
  e$PosAve<-e$Total*e$Pos  #Calculate x-weighted midpoint of averaging window, note that this essentially gives the mean position (rather than median). This is intentional.
  e$PosAve<-filter(e$PosAve,f)
  e$PosAve<-e$PosAve/filter(e$Total,f)

  b$WatsonAve <- as.numeric(b$WatsonAve) # COERCE from TIME-SERIES TO NUMERIC (IMPORTANT FOR PLOTTING)
  b$PosAve <- as.numeric(b$PosAve)
  d$CrickAve <- as.numeric(d$CrickAve)
  d$PosAve <- as.numeric(d$PosAve)
  e$TotalAve <- as.numeric(e$TotalAve)
  e$PosAve <- as.numeric(e$PosAve)
  e[is.na(e)] <- 0L; b[is.na(b)] <- 0L; d[is.na(d)] <- 0L; # convert NA to 0
  b <- b[seq(((Fsize-1)/2)+1, (nrow(b)-(Fsize-1)/2)),];   d <- d[seq(((Fsize-1)/2)+1, (nrow(d)-(Fsize-1)/2)),];   e <- e[seq(((Fsize-1)/2)+1, (nrow(e)-(Fsize-1)/2)),]  #top and tail the data

  b$WatsonAveNorm <- b$WatsonAve * norm.factor
  d$CrickAveNorm <-  d$CrickAve * norm.factor
  e$TotalAveNorm <-  (e$TotalAve * norm.factor)/2
  VX.out <- list(b,d,e)
  return(VX.out)
}
#' VX.raw
#' @description VariableX smoother which takes a sparse CCseq dataframe or datatable (colnames = c("Chr", "Pos", "Watson", "Crick)), and smooths both the y dimension ("Watson" and "Crick"), and the x dimension ("Pos"). A sliding mean height and mean position of Fsize number of hits is calculated.
#' @param x A sparse CCseq dataframe or datatable (colnames = c("Chr", "Pos", "Watson", "Crick)). E.g the "sae2.0" table generated within mapper.
#' @param Fsize The number of single strand sites (non-zero "Watson" or "Crick" rows) or double-stranded sites ("Total" rows) over which to smooth. This is the number of values that will be averaged, NOT the  distance in bp.
#' @author George Brown
#' @export
VX.raw  <- function(df,x,y, Fsize = 101, norm.factor = 3,keepdf=TRUE){
  f=rep(1,Fsize) #Create a SUM filter. Every valeu in the window of interest will be SUMMED (1 means ADD).
  Fdif=c(1,rep(0,Fsize-2),-1) #Create a filter that is will subtract the first from the last value on the window of interest. This calculates window width
  # x$Total <- x$Watson + x$Crick
  # x=c(1:1000)
  # y=rnorm(1000)
  # s[>.2]
  df[[x]]
  
  # x = rnorm(100, mean = 0, sd = 1)
  # b=data.frame(x,y)
  df=subset(df,df[[y]]!=0)
  df[[x]] <- as.numeric(df[[x]]); df[[y]] <- as.numeric(df[[y]]) # important to avoid integer overflow!
  df<-df[order(df[[x]]),]
  df[[paste0(x,"Dif")]]=stats::filter(df[[x]],Fdif) #Calculate width of averaging window
  df[[paste0(y,"Ave")]]=stats::filter(df[[y]],f)/df[[paste0(x,"Dif")]] #Divide total Watson signal in window by the window width
  df[[paste0(x,"Ave")]]<-df[[x]]*df[[y]]  #Calculate x-weighted midpoint of averaging window, note that this essentially gives the mean position (rather than median). This is intentional.
  df[[paste0(x,"Ave")]]<-stats::filter(df[[paste0(x,"Ave")]],f) #Calculate x-weighted midpoint of averaging window, note that this essentially gives the mean position (rather than median). This is intentional.
  df[[paste0(x,"VarX")]]<-df[[paste0(x,"Ave")]]/stats::filter(df[[y]],f)
  df <- df[seq(((Fsize-1)/2)+1, (nrow(df)-(Fsize-1)/2)),]
  # names(df)[names(df) == paste0(x,"Ave")] <- paste0(x,"VarX")
  df[[paste0(y,"VarX")]] <- df[[paste0(y,"Ave")]] * norm.factor
  drops <- c(paste0(x,"Dif"),paste0(y,"Ave"),paste0(x,"Ave"))
  df=df[ , !(names(df) %in% drops)]
  output=data.frame(df[[paste0(x,"VarX")]],df[[paste0(y,"VarX")]])
  # df=rename(df,paste0(x,"Ave")=,"b")
  colnames(output) <- c(x, y)
  if (keepdf == TRUE){return(df)
  }else{return(output)}
}
