#' CCSeqbias2
#' @description generates a sequence bias plot. This version is Will's work in progress, to make a faster, and more generic script.
#' @param path.name point to a folder containing fullmaps
#' @param minThresholdv A vector containing the min threshold
#' @param maxThreshold A value containing the max
#' @param w Window width in bp
#' @param rDNA include rDNA
#' @param chroms A Vector of chromosomes to process.
#' @param ylims min max of ylims. E.g c(0.05,0.6)
#' @param dyad Dyad to plot center lines from (0.5 for Spo11, 1.5 for Top2)
#' @param genome location to genome.fa
#' @param loci this can be set to a specific genename. The sequence bias will be calcuated only for sites within the genebody of this gene. HpM thresholding occurs before this, but you will probably still need to modify min max thresholds.
#' @param plotT  #PlotType: "l", "p", "o"
#' @author Will Gittens, George Brown, Matt Neale
#' @export
CCSeqbias2 <-function(path.name=getwd(),minThresholdv=c(0,1,5),maxThreshold="MAX",w=20,rDNA=F,chroms=c(1:16), ylims = c(0.05,0.6), dyad=1.5,genome="Cer3H4L2", plotT="o"){
  library(seqinr)
  seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"))
  #Read in CC-seq data
  parent.directory <- dirname(path.name)
  CCList(path.name)


  if(genome == "Cer3H4L2"){
    if(!rDNA){
      DSBList <- lapply(DSBList, function(x){
        x <- x[!(x$Chr == 12 & x$Pos >= 451000 & x$Pos <= 469000), ]
        return(x)
      })
    }
    }

  if(genome == "pombase220208"){
    if(!rDNA){
      DSBList <- lapply(DSBList, function(x){
        x <- x[!(x$Chr == 3 & x$Pos <= 24600 & x$Pos <= 2439550), ]
        return(x)
      })
       }
     }


# CCList("/Users/wg45/Dropbox/Lab_NGS/CC-seq/CC-seq_Saved/FullMaps/Raw/George/Test")

# Read in entire genome as list of character vectors. Only 97.4 mb of RAM used for the entire genome encoded this way!
  library(stringr)
  if (genome == "Cer3H4L2"){
  }else if(genome == "W303"){
    gen=CCAnnotate(W303_genome)
  }else if (genome == "Cer3H4L2MATa"){
    gen=CCAnnotate(Cer3H4L2MATa_genome)
  }else if (genome == "pombase220208"){
    gen  <- read.fasta("/Users/wg45/Dropbox/Work/termMapper_config/pombase220208/pombase220208.fa")
  }else{
    gen=read_fasta(genome)
  }

# gen <- read.fasta(file = "/Users/wg45/Dropbox/Work/termMapper_config/Cer3H4L2/Cer3H4L2.fa", seqtype = "DNA")
script="SeqBiasR_v03"

# threshold stuff
t2a = maxThreshold
if (maxThreshold=="MAX"|| 1000000){maxThreshold=1000000;t2a="MAX"}

for (j in minThresholdv){
  minThreshold = j

###Watson strand below
for (k in 1:length(DSBList)){
DSBList.2 <- split(DSBList[[k]], DSBList[[k]]$Chr) # subset the first sample and split it by chromosome
gen2 <- gen[chroms] # filter to include only desired chromosomes
DSBList.2 <- DSBList.2[chroms] # filter to include only desired chromosomes
# browser()
maxsitesF <- sum(sapply(DSBList.2, function(x){
nrow(x[x$Watson > 0,])
})
)

DSBList.F <- lapply(DSBList.2, function(x){
  x <- x[x$Watson>j & x$Watson <= maxThreshold,] # throw away rows below threshold
  return(x)
})

SampledSitesF <- sum(sapply(DSBList.F, nrow))
# browser()
# for each chromosome, generate a matrix where each row is a window centred on each watson hit
out <- lapply(DSBList.F, function(x){
  out <- t(seq2(x$Pos-w, (x$Pos+w)))
  out <- out[!rowSums(out <= 0), ] # remove windows that extend off lefthand side of the chromosome
  return(out)
})

#subset the genome by these coordinates, to rapidly convert the coordinate windows into sequence windows
out2 <- mapply(function(x,y) { # mapply the function chromosome-by-chromosome, because both genome and CC coordinates are organized in chromosome-separated lists
  z <- apply(y, 2, function(q){ # this applies the indexing function to the CC coordinates in a row-wise fashion. Not essential, but probably faster than column-wise.
    p <- x[q] # simple indexing: subset the DNA sequence on that chromosome by the CC coordinates
    return(p)
  })
  return(z)
}, gen2, out, SIMPLIFY = F) # these are x and y in the mapply function above
out4 <- do.call("rbind", out2) # bind all the chromosomes back together
out4 <- apply(out4, 2, table) # use super fast table function to quickly count number of bases at each position
# browser()
if(class(out4) == "matrix") {if(ncol(out4) > nrow(out4)){ out4 <- t(out4)}}

# browser()
# some posiitons may have some some "n" nucleotides and some may have none. This part makes sure that those that don't have any still get an "n" column
if(class(out4) == "matrix"){
  out4 <- as.data.frame(out4)
  if(!("n" %in% colnames(out4))) {
    out4[,"n"] <- 0
  }
}
if(class(out4) == "list"){
out4 <- lapply(out4, function(x){
  if(!("n" %in% names(x))) {
    x["n"] <- 0
  }
  x <- x[c("a", "c", "g", "t", "n")]
  return(x)
  }
)
out4 <- do.call("rbind", out4) # bind all the chromosomes back together
}
out4 <- out4/rowSums(out4) # divide by total nucleotides to give fraction at each position
row.names(out4) <- (-((nrow(out4)-1)/2):((nrow(out4)-1)/2))
outf <- out4
#############

###Crick strand below
maxsitesR <- sum(sapply(DSBList.2, function(x){
  nrow(x[x$Crick > 0,])
})
)
DSBList.R <- lapply(DSBList.2, function(x){
  x <- x[x$Crick>j & x$Crick <= maxThreshold,] # throw away rows below threshold
  return(x)
})

SampledSitesR <- sum(sapply(DSBList.R, nrow))

# for each chromosome, generate a matrix where each row is a window centred on each watson hit
out <- lapply(DSBList.R, function(x){
  out <- t(seq2(x$Pos-w, (x$Pos+w)))
  out <- out[!rowSums(out <= 0), ] # remove windows that extend off lefthand side of the chromosome
  return(out)
})

#subset the genome by these coordinates, to rapidly convert the coordinate windows into sequence windows
out2 <- mapply(function(x,y) { # mapply the function chromosome-by-chromosome, because both genome and CC coordinates are organized in chromosome-separated lists
  z <- apply(y, 2, function(q){ # this applies the indexing function to the CC coordinates in a row-wise fashion. Not essential, but probably faster than column-wise.
    p <- x[q] # simple indexing: subset the DNA sequence on that chromosome by the CC coordinates
    return(p)
  })
  return(z)
}, gen2, out, SIMPLIFY = F) # these are x and y in the mapply function above
out4 <- do.call("rbind", out2) # bind all the chromosomes back together
out4 <- apply(out4, 2, table) # use super fast table function to quickly count number of bases at each position
if(class(out4) == "matrix") {if(ncol(out4) > nrow(out4)){ out4 <- t(out4)}}

# some posiitons may have some some "n" nucleotides and some may have none. This part makes sure that those that don't have any still get an "n" column
if(class(out4) == "matrix"){
  out4 <- as.data.frame(out4)
if(!("n" %in% colnames(out4))) {
  out4[,"n"] <- 0
}
}

if(class(out4) == "list"){
out4 <- lapply(out4, function(x){
  if(!("n" %in% names(x))) {
    x["n"] <- 0
  }
  x <- x[c("a", "c", "g", "t", "n")]
  return(x)
}
)
out4 <- do.call("rbind", out4) # bind all the chromosomes back together
}
out4 <- out4/rowSums(out4) # divide by total nucleotides to give fraction at each position
row.names(out4) <- (-((nrow(out4)-1)/2):((nrow(out4)-1)/2))
outr <- out4[rev(1:nrow(out4)),]
row.names(outr) <- row.names(out4)
#############



## plotting module from Matt/George

SMF=round(SampledSitesF/maxsitesF,2)
SMR=round(SampledSitesR/maxsitesR,2)
SMT=round((SampledSitesR+SampledSitesF)/(maxsitesR+maxsitesF),2)

vert=seq(from =-round(w,-1), to=round(w,-1), by=5)+dyad #Add vertical grid lines centred on dyad

#########################################################################################################################################################
setwd(parent.directory)
ifelse(!dir.exists(file.path(parent.directory, "seqBias")), dir.create(file.path(parent.directory, "seqBias")), FALSE); setwd(file.path(parent.directory, "seqBias")); directory2 <- getwd()
ifelse(!dir.exists(file.path(directory2, DSBListNames[k])), dir.create(file.path(directory2, DSBListNames[k])), FALSE); setwd(file.path(directory2, DSBListNames[k])); directory2 <- getwd()
wd = getwd();
# if(missing(loci)) {
  # out = paste(wd,"/","seqBias","/",DSBListNames[k],"/",script,"_",DSBListNames[k],"_rDNA=",rDNA,"_±",w,"_T>",minThreshold, "to",t2a,"_HpM_dyad=",dyad,"_C_",paste(chroms, sep = "_", collapse = "_"),".pdf",sep="")
out = paste(wd,"/",script,"_",DSBListNames[k],"_rDNA=",rDNA,"_±",w,"_T>",minThreshold, "to",t2a,"_HpM_dyad=",dyad,".pdf",sep="")
  # }
# if(!missing(loci)) {out = paste(wd,"/","seqBias","/",DSBListNames[k],"/",script,"_",DSBListNames[k],"_rDNA=",rDNA,"_±",w,"_T>",minThreshold, "to",t2a,"_HpM_dyad=",dyad,"_",loci, "only.pdf",sep="")}
; pdf(file=out, width=9,height=18) # Plot PDF
layout(matrix(c(1,1,1,2,2,2,3,3,3),9, 1, byrow = T))

#Plotting_F
mreads.1 <- round(Mreads[k],2)[[1]]
DSBListNames[k]
ymin <- ylims[1]
ymax <- ylims[2]
# plot(0,0,xlim=c(-w,w), ylim=ylims, type="n")
# if(missing(loci)) {
  plot(0,0,xlim=c(-w,w), ylim=c(ymin,ymax), type="n", xlab="Position (Zero is the 5'end)", ylab="Base Frequency", main=paste0(DSBListNames[k]," Mreads = ",mreads.1,"\nFranklin_Sites = ",SampledSitesF,"(",SMF*100,"%)"," / Thresh > ",minThreshold, " HpM to ",t2a," HpM / Chroms = ", paste0(chroms, collapse = "-")))
  # }
# if(!missing(loci)) {plot(0,0,xlim=c(-w,w), ylim=c(ymin,ymax), type="n", xlab="Position (Zero is the 5'end)", ylab="Base Frequency", main=paste0(DSBListNames[k]," Mreads = ",round(Mreads[k],2),"\nFranklin_Sites = ",SampledSitesF,"(",SMF*100,"%)"," / Thresh > ",minThreshold, " HpM to ",t2a," HpM / Gene = ", loci))}
polygon(c(-w,w,w,-w),c(ymin, ymin, ymax, ymax),col="grey", border=NA)
abline(v=dyad, col="white", lwd=4)
abline(v=vert, col="white")
lines(row.names(outf), outf[,"a"], type=plotT, pch=19, col="yellow", lwd=4)
lines(row.names(outf), outf[,"c"], type=plotT, pch=19, col="green", lwd=4)
lines(row.names(outf), outf[,"g"], type=plotT, pch=19, col="red", lwd=4)
lines(row.names(outf), outf[,"t"], type=plotT, pch=19, col="blue", lwd=4)
legend(w*0.85, ymax*0.95, legend=c("A", "T", "C", "G"), col=c("yellow", "blue", "green", "red"), lty=1, cex=1, bg="white")

#Plotting_R

# plot(0,0,xlim=c(-w,w), ylim=ylims, type="n")
# if(missing(loci)) {
  plot(0,0,xlim=c(-w,w), ylim=c(ymin,ymax), type="n", xlab="Position (Zero is the 5'end)", ylab="Base Frequency", main=paste0(DSBListNames[k]," Mreads = ",mreads.1,"\nRosalind_Sites = ",SampledSitesF,"(",SMF*100,"%)"," / Thresh > ",minThreshold, " HpM to ",t2a," HpM / Chroms = ", paste0(chroms, collapse = "-")))
  # }
# if(!missing(loci)) {plot(0,0,xlim=c(-w,w), ylim=c(ymin,ymax), type="n", xlab="Position (Zero is the 5'end)", ylab="Base Frequency", main=paste0(DSBListNames[k]," Mreads = ",round(Mreads[k],2),"\nRosalind_Sites = ",SampledSitesF,"(",SMF*100,"%)"," / Thresh > ",minThreshold, " HpM to ",t2a," HpM / Gene = ", loci))}
polygon(c(-w,w,w,-w),c(ymin, ymin, ymax, ymax),col="grey", border=NA)
abline(v=dyad, col="white", lwd=4)
abline(v=vert, col="white")
lines(row.names(outr), outr[,"t"], type=plotT, pch=19, col="yellow", lwd=4) # note the bases are complementary here relative to the above strand
lines(row.names(outr), outr[,"g"], type=plotT, pch=19, col="green", lwd=4)
lines(row.names(outr), outr[,"c"], type=plotT, pch=19, col="red", lwd=4)
lines(row.names(outr), outr[,"a"], type=plotT, pch=19, col="blue", lwd=4)
legend(w*0.85, ymax*0.95, legend=c("A", "T", "C", "G"), col=c("yellow", "blue", "green", "red"), lty=1, cex=1, bg="white")

#Plotting_F+R
# plot(0,0,xlim=c(-w,w), ylim=ylims, type="n")
# if(missing(loci)) {
  plot(0,0,xlim=c(-w,w), ylim=c(ymin,ymax), type="n", xlab="Position (Zero is the 5'end)", ylab="Base Frequency", main=paste0(DSBListNames[k]," Mreads = ",mreads.1,"\nCombined_Sites = ",SampledSitesF,"(",SMF*100,"%)"," / Thresh > ",minThreshold, " HpM to ",t2a," HpM / Chroms = ", paste0(chroms, collapse = "-")))
  # }
# if(!missing(loci)) {plot(0,0,xlim=c(-w,w), ylim=c(ymin,ymax), type="n", xlab="Position (Zero is the 5'end)", ylab="Base Frequency", main=paste0(DSBListNames[k]," Mreads = ",round(Mreads[k],2),"\nCombined_Sites = ",SampledSitesF,"(",SMF*100,"%)"," / Thresh > ",minThreshold, " HpM to ",t2a," HpM / Gene = ", loci))}
polygon(c(-w,w,w,-w),c(ymin, ymin, ymax, ymax),col="grey", border=NA)
abline(v=dyad, col="white", lwd=4)
abline(v=vert, col="white")
lines(row.names(outf), (outf[,"a"] + outr[,"t"])/2, type=plotT, pch=19, col="yellow", lwd=4)
lines(row.names(outf), (outf[,"c"] + outr[,"g"])/2, type=plotT, pch=19, col="green", lwd=4)
lines(row.names(outf), (outf[,"g"] + outr[,"c"])/2, type=plotT, pch=19, col="red", lwd=4)
lines(row.names(outf), (outf[,"t"] + outr[,"a"])/2, type=plotT, pch=19, col="blue", lwd=4)
legend(w*0.85, ymax*0.95, legend=c("A", "T", "C", "G"), col=c("yellow", "blue", "green", "red"), lty=1, cex=1, bg="white")

##########################################################################################################################################################
dev.off() #Close connection to PDF writer
}
}
}
