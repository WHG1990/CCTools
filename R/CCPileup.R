#' @title CCPileup
#' @description Maps CCs in a specific chromosome region, alongside various annotations.
#' @param in.mode Reading in .RDS file or a folder of FullMaps? Choose "rds" or "fullmap".  Defaults to "rds".
#' @param path.name A character string defining the RDS file or folder containing FullMaps.
#' @param loci The path to an RDS file or .txt file containing a table of loci in column format: c("Chr", "Pos", "Strand"). Alternatively can use inbuilt loci such as "TSS" (From AllElementsDUB), or "CTSS" (From YeastTSS), "CEN" (Centromeres). More info on the generation of these sets is available from Will.
#' @param exp.name Manually specify a name for the experiment.
#' @param genome Which genome has the data been aligned to? "Cer3H4L2", "W303" or "hg19". Defaults to "Cer3H4L2".
#' @param phase Which phase is it? "meiosis" or "vegetative"
#' @param samples Which samples do you want to plot? Provide a numeric vector of DSBList indexes. Defaults to all samples.
#' @param window.w How wide a region to pileup in bp?
#' @param samples Which samples to plot, as a vector of indexes in the DSBList.
#' @param combine Do you want to plot all samples on the same plot? If so, it is recommended to set strat.mode to "none", or the plots will be very busy. Defaults to FALSE
#' @param strat.mode How do you want to stratify the plots? "none" will leave the Pileups unstratified. "manual" will point the stratifier to a 4th numeric column (column name must be "Value") in the user-defined loci table. "expression" or "length" can be used in combination with the inbuilt TSS or CTSS loci tables, and refer to microarray gene expression, and gene length, respectively.
#' @param strat.levels How many strata do you want to use? Defaults to 4.
#' @param out.mode "pdf" or "png" file format for output.
#' @param ylims A single number specifying the maximum limit of the y axis, for both Watson and Crick. OR "auto", to automatically define for each sample.
#' @param win Hanning window width for smoothing.
#' @param plot.dimensions Specify the plot size. Can be any size from the "A sizing" scale (e.g "A1"). However this is not implemented properly yet. Defaults to "wide", which is fine.
#' @param plot.w Specify the plot range in bp. Defaults to the pileup window.w.
#' @param mirror.mode Do you want to average the plot either side of the loci? This generates a symmetrical plot, which is useful when considering loci with no inherent directionality (e.g centromeres). Defaults to FALSE.
#' @param strand.sep Do you want to plot the Top and Bottom strands (relative to locus) separately? Defaults to FALSE
#' @param positive.strand Do you want to plot both Top and Bottom strands on the same positive axis? If so, TRUE. If you want to plot Top on positive and bottom on negative, then change to FALSE.
#' @param cache.search Do you want to search in the folders surrounding the path.name for a matching saved PileUp? Defaults to true.
#' @param normalize Do you want to normalize every pileup over the range of window.w? (NOTE: IT NORMALISES OVER window.w , not plot.w!)
#' @examples
#' @author Will Gittens
#' @import data.table
#' @export
CCPileup <- function(path.name, exp.name, in.mode = "fullmap",  loci = "TSS", genome = "Cer3H4L2", phase = "meiosis", samples, window.w = 5000, combine = F, strat.mode = "none", strat.levels = 4, out.mode = "png", ylims = "auto", win = 50, plot.dimensions = "wide", mirror.mode = F, strand.sep = F, plot.w, positive.strands = T, normalize = T, cache.search = T, names, col, rDNA = F, mtDNA = F) {

  library("data.table")
  library("e1071")
  library("scales")
  library("RColorBrewer")

  chroms=1:16 # Or just a single chromosome. Hash this line out if you want to plot all of them specified above.
  width1=window.w/2 # bp upstream and downstream of motif
  padding=500
  width=width1+padding
  ################# data loading

  parent.directory <- dirname(path.name)
  ifelse(!dir.exists(file.path(parent.directory, "PileUp")), dir.create(file.path(parent.directory, "PileUp")), FALSE); setwd(file.path(parent.directory, "PileUp")); directory2 <- getwd()
  ifelse(!dir.exists(file.path(directory2, "/",exp.name)), dir.create(file.path(directory2, "/",exp.name)), FALSE); setwd(file.path(directory2, "/", exp.name)); directory2 <- getwd()
  files <- list.files(pattern = "PileUp")

  if(in.mode == "rds") {
    load(path.name)
    exp.name <- substr(basename(path.name), 0, nchar(basename(path.name))-4 )
  }

  if(in.mode == "fullmap") {
    CCList(path.name)
    exp.name <- exp.name
  }

  # samples <- readline("Which samples to plot?")
  # samples <- as.vector(samples)

  # names(DSBListNames) <- 1:length(DSBListNames)
  # print(DSBListNames)
  #   question <- paste0("Which samples do you want to plot? Answer as a numeric vector.")
  #   choices <- 1:length(DSBListNames)
  #   samples <- menu(choices, graphics = FALSE, question)

  if(missing(samples)){samples <- 1:length(DSBList)}
  if(samples[1] == "all"){samples <- 1:length(DSBList)}
  if(!missing(names)){DSBListNames <- names}
  if(missing(plot.w)) {plot.w = window.w}

  # DSBList <- DSBList[samples]; DSBListNames <- DSBListNames[samples]
  DSBList <- lapply(DSBList, data.table)
  DSBList<- lapply(DSBList, setkey, Chr) # set the keys as Chr and Pos (data.table)
  Mreads <- sapply(DSBList, function(x) {sum(x$Watson + x$Crick)/1000000})


  #user specified loci
  if(strat.mode == "none"){ strat.levels = 1}
  if(substr(loci, nchar(loci) - 3, nchar(loci)) == ".rds") {
    motifs <- load(loci)
  }
  if (substr(loci, nchar(loci) - 3, nchar(loci)) == ".txt") {
    motifs <- fread(loci)
  }


  if (genome == "Cer3H4L2") {

    #remove rDNA and mtDNA

    if(!rDNA){
      DSBList <- lapply(DSBList, function(x){
        x <- x[!(x$Chr == 12 & x$Pos >= 451000 & x$Pos <= 469000),]
        return(x)
      })
    }

    if(!mtDNA){
      DSBList <- lapply(DSBList, function(x){
        x <- x[x$Chr != 17,]
        return(x)
      })
    }


  #inbuilt TSS loci
  if (loci == "TSS") {
    if (phase == "vegetative") {
      motifs <-
        CCAnnotate(Cer3_TSS_positions_forwhich_GSE36958_vegetative_available)
    }
    if (phase == "meiosis") {
      motifs <-
        CCAnnotate(Cer3_TSS_positions_forwhich_GSE36958_meiosis_available)
    }
    if(strat.mode == "manual") { warning("Manual stratification mode is not implemented for inbuilt loci, defaulting 'none'")
      strat.mode = "none"}
  }

  #inbuilt TTS loci
  if (loci == "TTS") {
    if (phase == "vegetative") {
      motifs <-
        CCAnnotate(Cer3_TSS_positions_forwhich_GSE36958_vegetative_available)
      motifs$Pos <- motifs$Pos2
    }
    if (phase == "meiosis") {
      motifs <-
        CCAnnotate(Cer3_TSS_positions_forwhich_GSE36958_meiosis_available)
      motifs$Pos <- motifs$Pos2
    }
    if(strat.mode == "manual") { warning("Manual stratification mode is not implemented for inbuilt loci, defaulting 'none'")
      strat.mode = "none"}
  }

  #inbuilt CTSS loci
  if (loci == "CTSS") {
    if (phase == "vegetative") {
      motifs <-
        CCAnnotate(Cer3H4L2_CTSS_positions_forwhich_GSE36958_vegetative_available)
    }
    if (phase == "meiosis") {
      motifs <-
        CCAnnotate(Cer3H4L2_CTSS_positions_forwhich_GSE36958_meiosis_available)
    }
    if(strat.mode == "manual") { warning("Manual stratification mode is not implemented for inbuilt loci, defaulting 'none'")
      strat.mode = "none"
    }
  }

  #inbuilt CEN loci
  if (loci == "CEN") {
    motifs <- CCAnnotate(Cer3H4L2_centro)
    motifs$Pos <- ceiling((motifs$chromStart+motifs$chromEnd)/2)
    motifs$Strand <- rep("+", nrow(motifs))
    motifs <- motifs[,c(1,5,6,4)]
    colnames(motifs) <- c("Chr", "Pos", "Strand", "Value")
    strat.mode = "none"
  }

  #inbuilt Nucleosome loci
  if (loci == "nucleosome") {
    motifs <- CCAnnotate(GSE36063_nucleoChemUnique)
    colnames(motifs) <- c("Chr", "Pos", "Strand", "Value")

  }

  #inbuilt ARS loci
  if (loci == "ARS") {
    motifs <- CCAnnotate(ARS_consensus_noChr12)
    motifs$Value <- rep(1, nrow(motifs))
    colnames(motifs) <- c("Chr", "Pos", "Strand", "Value")

  }

  #inbuilt HS loci
  if (loci == "HS_start") {
    motifs <- CCAnnotate(PanHS)
    motifs <- motifs[, c("CHROM", "midpoint", "HITS","LENGTH")]
    colnames(motifs) <- c("Chr", "Pos", "Value","Width")
    motifs$Strand <- rep("+", nrow(motifs))
    motifs <- motifs[, c(1, 2, 5, 4,3)]
  }
  if (loci == "HS_mid") {
    motifs <- CCAnnotate(PanHS)
    motifs <- motifs[, c("CHROM", "midpoint", "HITS","LENGTH")]
    colnames(motifs) <- c("Chr", "Pos", "Value","Width")
    motifs$Strand <- rep("+", nrow(motifs))
    motifs <- motifs[, c(1, 2, 5, 4,3)]
    ot<<-motifs
  }


  }

  #######

  if(genome == "W303"){
    motifs <- CCAnnotate(W303_centro)
    motifs$Pos <- ceiling((motifs$chromStart+motifs$chromEnd)/2)
    motifs$Strand <- rep("+", nrow(motifs))
    motifs <- motifs[,c(1,5,6,4)]
    colnames(motifs) <- c("Chr", "Pos", "Strand", "Value")
    strat.mode = "none"
  }

  #######

  if (genome == "pombase220208") {

    if(!rDNA){
      DSBList <- lapply(DSBList, function(x){
        x <- x[!(x$Chr == 3 & x$Pos <= 24600),]
        x <- x[!(x$Chr == 3 & x$Pos >= 2439550),]
        return(x)
      })
    }

    if(!mtDNA){
      DSBList <- lapply(DSBList, function(x){
        x <- x[x$Chr != 4,]
        return(x)
      })
    }


    #inbuilt TSS loci
    if (loci == "TSS") {
      if (phase == "vegetative") {
         motifs <- CCAnnotate(pombase220208_TSS_TTS)
         motifs$Width <- motifs$stop-motifs$start

      }
      # if (phase == "meiosis") {
      #   motifs <-
      #     CCAnnotate(Cer3_TSS_positions_forwhich_GSE36958_meiosis_available)
      # }
      if(strat.mode == "manual") { warning("Manual stratification mode is not implemented for inbuilt loci, defaulting 'none'")
        strat.mode = "none"}
    }

    #inbuilt TTS loci
    if (loci == "TTS") {
      if (phase == "vegetative") {
        motifs <- CCAnnotate(pombase220208_TSS_TTS)
        motifs$Pos <- motifs$Pos2
        motifs$Width <- motifs$stop-motifs$start
      }
      # if (phase == "meiosis") {
      #   motifs <-
      #     CCAnnotate(Cer3_TSS_positions_forwhich_GSE36958_meiosis_available)
      #   motifs$Pos <- motifs$Pos2
      # }
      if(strat.mode == "manual") { warning("Manual stratification mode is not implemented for inbuilt loci, defaulting 'none'")
        strat.mode = "none"}
    }

    # #inbuilt CTSS loci
    # if (loci == "CTSS") {
    #   if (phase == "vegetative") {
    #     motifs <-
    #       CCAnnotate(Cer3H4L2_CTSS_positions_forwhich_GSE36958_vegetative_available)
    #   }
    #   if (phase == "meiosis") {
    #     motifs <-
    #       CCAnnotate(Cer3H4L2_CTSS_positions_forwhich_GSE36958_meiosis_available)
    #   }
    #   if(strat.mode == "manual") { warning("Manual stratification mode is not implemented for inbuilt loci, defaulting 'none'")
    #     strat.mode = "none"
    #   }
    # }

    #inbuilt CEN loci
    if (loci == "CEN") {
        motifs <- CCAnnotate(pombase220208_centro)
        motifs$Pos <- ceiling((motifs$chromStart+motifs$chromEnd)/2)
        motifs$Strand <- rep("+", nrow(motifs))
        motifs <- motifs[,c(1,5,6,4)]
        colnames(motifs) <- c("Chr", "Pos", "Strand", "Value")
        strat.mode = "none"
    }

    # #inbuilt Nucleosome loci
    # if (loci == "nucleosome") {
    #   motifs <- CCAnnotate(GSE36063_nucleoChemUnique)
    #   colnames(motifs) <- c("Chr", "Pos", "Strand", "Value")
    #
    # }

    #inbuilt ARS loci
    if (loci == "ARS") {
      motifs <- fread("/Users/wg45/Dropbox/Work/Top2Seq/Genomes_and_Annotations/Pombe/Features/pombase_origins.tsv")
      motifs$Pos <- ceiling(rowMeans(motifs[,c("start", "end")]))
      motifs$Chr <- motifs$chr
      motifs$Strand <- "+"
      motifs$Value <- 1
      motifs <- motifs[motifs$status != "Dubious",]
      motifs <- motifs[,c("Chr", "Pos", "Strand", "Value")]
          }

  #   #inbuilt HS loci
  #   if (loci == "HS_start") {
  #     motifs <- CCAnnotate(PanHS)
  #     motifs <- motifs[, c("CHROM", "midpoint", "HITS","LENGTH")]
  #     colnames(motifs) <- c("Chr", "Pos", "Value","Width")
  #     motifs$Strand <- rep("+", nrow(motifs))
  #     motifs <- motifs[, c(1, 2, 5, 4,3)]
  #   }
  #   if (loci == "HS_mid") {
  #     motifs <- CCAnnotate(PanHS)
  #     motifs <- motifs[, c("CHROM", "midpoint", "HITS","LENGTH")]
  #     colnames(motifs) <- c("Chr", "Pos", "Value","Width")
  #     motifs$Strand <- rep("+", nrow(motifs))
  #     motifs <- motifs[, c(1, 2, 5, 4,3)]
  #     ot<<-motifs
  #   }
  #   if (loci == "HS_end") {
  #     motifs <- CCAnnotate(PanHS)
  #     motifs <- motifs[, c("CHROM", "midpoint", "HITS","LENGTH")]
  #     colnames(motifs) <- c("Chr", "Pos", "Value","Width")
  #     motifs$Strand <- rep("+", nrow(motifs))
  #     motifs <- motifs[, c(1, 2, 5, 4,3)]
  #   }
  # }


  #inbuilt Mnase loci. NEED TO DO SOME PEAK CALLING FIRST
  # if (loci == "MNase") {
  #   if (phase == "vegetative") {
  #     motifs <-
  #       CCAnnotate(MNase_Veg_GSM1849297_GSM1849298_R1R2)
  #   }
  #   if (phase == "meiosis") {
  #     motifs <-
  #       CCAnnotate(MNase_meiosis_GSM1424408)
  #   }
  #   if(strat.mode == "manual") { warning("Manual stratification mode is not implemented for inbuilt loci, defaulting 'none'")
  #     strat.mode = "none"}
  # }
  }

  if((strat.mode == "none") || (strat.levels == 1)) {
    motifs <- motifs[, c("Chr", "Pos", "Strand")]
    motifs$Value <- rep(1,nrow(motifs))
    motifs <- split(motifs, dplyr::ntile(motifs$Value, 1))
  } else {

    if ((class(motifs) != "list") && (strat.mode == "manual")) {
      if (all(colnames(motifs)[1:3] != c("Chr", "Pos", "Strand"))) {
        stop(
          "The loci table is not formated correctly. The first three column names must be Chr, Pos and Strand. The optional fourth column name must be Value"
        )
      }
      motifs <- split(motifs, dplyr::ntile(motifs$Value, strat.levels))

    }

    if ((class(motifs) != "list") && (strat.mode == "expression")) {
      motifs <- motifs[, c("Chr", "Pos", "Strand", "ExpressionValue")]
      colnames(motifs) <- c("Chr", "Pos", "Strand", "Value")
      motifs <- split(motifs, dplyr::ntile(motifs$Value, strat.levels))
    }

    if ((class(motifs) != "list") && (strat.mode == "length")) {
      if (loci %in% c("CTSS", "CEN")) {
        warning("Length stratification mode is not implemented for CTSS or CEN, defaulting 'none'")
        motifs <- motifs[, c("Chr", "Pos", "Strand")]
        motifs$Value <- rep(1,nrow(motifs))
        motifs <- split(motifs, dplyr::ntile(motifs$Value, 1))
      } else {
        motifs <- motifs[, c("Chr", "Pos", "Strand", "Width")]
        colnames(motifs) <- c("Chr", "Pos", "Strand", "Value")
        motifs <- split(motifs, dplyr::ntile(motifs$Value, strat.levels))
      }
    }
  }

  if(!(loci %in% c("TSS", "TTS", "CTSS", "CEN", "ARS"))) {
    loci.name <- substr(basename(loci), 0, nchar(basename(loci))-4)
  } else {loci.name = loci}

  ifelse(!dir.exists(file.path(directory2, "/",loci.name)), dir.create(file.path(directory2, "/",loci.name)), FALSE); setwd(file.path(directory2, "/", loci.name)); directory2 <- getwd()

  if(strat.mode != "none") {loci.name = paste0(strat.mode, "_stratified", strat.levels, "_", loci.name)}
  if(strat.mode == "none") {loci.name = paste0("unstratified_", loci.name)}
  filename <- paste0("PileUp_",exp.name, "_on_", loci.name,"_", window.w/1000,"KB")
  setwd(parent.directory)
  setwd(directory2)

  #   CCSearch <- function(path.name, file.name){
  #   if(missing(path.name)){path.name <- getwd()}
  #   path.vec <- strsplit(path.name, "/")[[1]]
  #   path.vec <- path.vec[2:length(path.vec)]
  #   path.vec2 <- path.vec
  #   for (i in 1:length(path.vec)){
  #     temp <- c(path.vec[1:i])
  #     temp <- paste0(temp, sep = "/", collapse = "")
  #     temp <- paste0("/", temp)
  #     path.vec2[i] <- temp
  #   }
  #   path.vec <- path.vec2
  #   path.vec <- path.vec[(ceiling(length(path.vec)/2)):length(path.vec)]
  #   path.vec <- rev(path.vec)
  #   for (i in 1:length(path.vec)){
  #     str <- paste0('list.files(path = path.vec[i], pattern = "', file.name, '", recursive = T)')
  #     out <<- eval(parse(text = str))
  #     if (length(out) > 0) { break }
  #   }
  #   return(paste0(path.vec[i],out))
  # }


filename2 <- paste0(filename, ".rds")
  cached.file <- character(length = 0L)
  if(cache.search == T){
    cached.file <- CCSearch(path.name = directory2, file.name = filename2)
    if(length(cached.file) > 0){
            question <- paste0("Cached PileUp(s) with this experiment name were found. Which do you want to use? Alternatively, rerun function with cache.search = F")
      answer <- menu(cached.file, graphics = F, question)
      load(cached.file[answer])
    } else {
      print("No cached PileUp found, running de novo.")
    }
  }
    if(cache.search == F | length(cached.file) == 0) {

  motifs <- lapply(motifs, data.table)
  lapply(motifs, setkey, Strand)
  motifsW<-lapply(motifs, subset, Strand=="+")
  motifsC<-lapply(motifs, subset, Strand=="-")
  lapply(motifsW, setkey, Chr); lapply(motifsC, setkey, Chr)

  #################
  PileUp_list <- rep(list(list()), length(DSBList))
  for (p in 1:length(DSBList)){
    for (q in 1:length(motifs)) {

      #Initialise PileUp dataframe
      PileUp=data.frame(NULL)
      PileUp[1:((width*2)+1),"Pos"]=(1:((2*width)+1))
      PileUp$WatsonW=0 # Watson MOTIFs watson hits
      PileUp$CrickW=0 # Watson MOTIFs crick hits
      PileUp$WatsonC=0 # Crick MOTIFs watson hits
      PileUp$CrickC=0 # Crick MOTIFs crick hits

      temp=NULL
      #### Start of Chromosome loop

      for (chromo in chroms) { # Specify which chromosomes to process
        df.1 <- DSBList[[p]][J(chromo)]
        motifsWtemp <- motifsW[[q]][J(chromo)]
        motifsCtemp <- motifsC[[q]][J(chromo)]
        setkey(df.1, Pos); setkey(motifsWtemp,NULL);  setkey(motifsCtemp,NULL)

        if(!is.na(motifsWtemp[1,"Pos"])){
          for (i in 1:nrow(motifsWtemp)) {
            temp <- df.1[J((as.numeric(motifsWtemp[i,"Pos"])-width):(as.numeric(motifsWtemp[i,"Pos"])+width)), nomatch=0L]
            temp$Pos=temp$Pos-motifsWtemp[[i,"Pos"]]+width+1
            PileUp[temp$Pos,"WatsonW"]=PileUp[temp$Pos,"WatsonW"]+temp[,"Watson"]
            PileUp[temp$Pos,"CrickW"]=PileUp[temp$Pos,"CrickW"]+temp[,"Crick"]
          }
        }
        if(!is.na(motifsCtemp[1,"Pos"])){
          for (i in 1:nrow(motifsCtemp)) {
            temp <- df.1[J((as.numeric(motifsCtemp[i,"Pos"])-width):(as.numeric(motifsCtemp[i,"Pos"])+width)), nomatch=0L]
            temp$Pos=temp$Pos-motifsCtemp[[i,"Pos"]]+width+1
            PileUp[temp$Pos,"WatsonC"]=PileUp[temp$Pos,"WatsonC"]+temp[,"Watson"]
            PileUp[temp$Pos,"CrickC"]=PileUp[temp$Pos,"CrickC"]+temp[,"Crick"]
          }
        }
      }

      # This section combines the Watson and Crick hits, reversing the order of the crick gene hits and adding the Watsons to the Crick hits and Crick to Watsons (i.e. reverse complements the data for the Crick MOTIFs).
      # This allows both to be combined in the plot
      PileUp$revWc=rev(PileUp$WatsonC) #Crick MOTIFs with Watson hits
      PileUp$revCc=rev(PileUp$CrickC) #Crick MOTIFs with Crick hits
      PileUp$WatsonTotal=PileUp$WatsonW+PileUp$revWc
      PileUp$CrickTotal=PileUp$CrickW+PileUp$revCc
      PileUp$TopTotal = PileUp$WatsonW + PileUp$revCc
      PileUp$BottomTotal = PileUp$CrickW + PileUp$revWc
      PileUp$WCTotal=(PileUp$WatsonTotal+PileUp$CrickTotal) #Create total column
      PileUp$Pos<-PileUp$Pos-(width+1)
      PileUp$WC.HpM <- PileUp$WCTotal/(Mreads[p]*(nrow(motifsW[[q]]) + nrow(motifsW[[q]])))
      PileUp$Top.HpM <- PileUp$TopTotal/ (Mreads[p]*(nrow(motifsW[[q]]) + nrow(motifsW[[q]])))
      PileUp$Bottom.HpM <- PileUp$BottomTotal/ (Mreads[p]*(nrow(motifsW[[q]]) + nrow(motifsW[[q]])))
      PileUp_list[[p]][[q]] <- PileUp
      # write.table(PileUp, file = paste0("PileUp_datatable_", (2*width1)/1000,"KB_", DSBListNames[p], "_", q,".txt"), row.names=FALSE, sep="\t")
    }
  }

  save(PileUp_list, file = paste0(filename, ".rds"))

    }
  ##########################################################################################################################################################
  # this bit calculates genome size and expected density
  # directory1 <- getwd(); setwd("/Users/wg45/Dropbox/Work/Top2Seq/WD")
  # load("Cer3H4L2_exclusionlist.rds")
  # setwd(parent.directory)
  # excluded.size <- 0
  # expected.density <- 1000000/((12100000-excluded.size)/1)

  recurse <- function (L, f) {
    if (inherits(L, "data.frame")) f(L)
    else lapply(L, recurse, f)
  }

  data <- PileUp_list
  # if(mirror.mode){
  #   data <- recurse(data, function(x){
  #     x$WC.HpM <- (x$WC.HpM+ (rev(x$WC.HpM)))/2
  #     x$Top.HpM <- (x$Top.HpM+ (rev(x$Top.HpM)))/2
  #     x$Bottom.HpM <- (x$Bottom.HpM+ (rev(x$Bottom.HpM)))/2
  #     return(x)
  #   })
  # }

  hw=hanning.window(win) #create hanning window (require package e1071 to be loaded)
  scalar=1/sum(hw) #scalar [1]

  # if(normalize){
  #   data <- recurse(data, function(x){
  #     x$WC.HpM <- x$WC.HpM/sum(x$WC.HpM)
  #     x$Top.HpM <-  x$Top.HpM/sum(x$Top.HpM)
  #     x$Bottom.HpM <-  x$Bottom.HpM/sum(x$Bottom.HpM)
  #     return(x)
  #   })
  # }


  data <- recurse(data, function(x){
    x$Smooth <- (stats::filter(x$WC.HpM, hw))*scalar
    x$Top.Smooth <- (stats::filter(x$Top.HpM, hw))*scalar
    x$Bottom.Smooth <- (stats::filter(x$Bottom.HpM, hw))*scalar
    return(x)
  })

  if(mirror.mode){
    data <- recurse(data, function(x){
      x$Smooth <- (x$Smooth+ (rev(x$Smooth)))/2
      x$Smooth <- (x$Smooth+ (rev(x$Smooth)))/2
      x$Smooth <- (x$Smooth+ (rev(x$Smooth)))/2
      return(x)
    })
  }

  data <- recurse(data, function(x){
    x$Smooth[is.na(x$Smooth)] <- 0
    x$Top.Smooth[is.na(x$Top.Smooth)] <- 0
    x$Bottom.Smooth[is.na(x$Bottom.Smooth)] <- 0
    return(x)
  })

  if(normalize){
  data <- recurse(data, function(x){
    x$Smooth <- x$Smooth/sum(x$Smooth)
    x$Top.Smooth <-  x$Top.Smooth/sum(x$Top.Smooth)
    x$Bottom.Smooth <-  x$Bottom.Smooth/sum(x$Bottom.Smooth)
    return(x)
  })
  }


  # plot.name <- paste0(DSBListNames[samples], sep = "_&_", collapse = "")
min2 <- function(x) {min(x[x > 0])}

  if(!strand.sep){
    function1 <- function(x) {
      max(x$Smooth[!is.na(x$Smooth)])
    }
    ylim2<-1.1*(max(unlist(recurse(data, function1))))

    function1 <- function(x) {
      min2(x$Smooth[!is.na(x$Smooth)])
    }
    ylim1<-0.9*(min2(unlist(recurse(data, function1))))

    if (ylims[1] != "auto"){
      if(length(ylims) == 2) {ylim1 <- ylims[1]; ylim2 <- ylims[2]}
      if(length(ylims) == 1) {ylim1 <- ylims[1]; ylim2 <- ylims[1]}
    }

        if(!combine){
      plot.colours <- rev(RColorBrewer::brewer.pal(4, "Spectral")); myColor.name <- "Spectral"
      if(!missing(col)){plot.colours <- col}
      plot.colours <- colorRampPalette(plot.colours)(strat.levels)
      for (i in samples) {
        {pdf(width = 11, height = 8, file = paste0("PileUp_",DSBListNames[i], "_on_", loci.name, "_", plot.w/1000, "KB_HW" , win, "bp_ylim", round(ylim2,2),"_Mirror=", mirror.mode,"_Norm=", normalize, "_StrandSep=", strand.sep, ".pdf"))

          plot(data[[i]][[1]]$Pos, yaxs = "i", xaxs = "i", ylim = c(ylim1, ylim2), xlim = c(-plot.w/2,plot.w/2), type = "n", ylab = paste0("Aggregated CC-seq signal (",win, " bp smoothed) / HpM"), xlab = "Position relative to feature midpoint/ bp")

          for(q in 1:strat.levels){
            lines(data[[i]][[q]]$Pos, data[[i]][[q]]$Smooth, lwd = 2, type = "l", col = plot.colours[q])
            if(strat.levels > 1){
              if(strat.mode == "manual") text((-plot.w/2), ylim2-(q*(ylim2/20)), labels = paste0(ordinal(q)," Quartile of Value"), col = plot.colours[q], pos = 4)
              if(strat.mode == "expression") text((-plot.w/2), ylim2-(q*(ylim2/20)), labels = paste0(ordinal(q)," Quartile of Expression"), col = plot.colours[q], pos = 4)
              if(strat.mode == "length") text((-plot.w/2), ylim2-(q*(ylim2/20)), labels = paste0(ordinal(q)," Quartile of Length"), col = plot.colours[q], pos = 4)
            } else {
              text((-plot.w/2), ylim2-(q*(ylim2/20)), labels = DSBListNames[i], col = plot.colours[q], pos = 4)
            }
          }
          dev.off()}
      }
    }

    if(combine){
      plot.colours <- rev(RColorBrewer::brewer.pal(4, "Spectral")); myColor.name <- "Spectral"
      if(!missing(col)){plot.colours <- col}
      plot.colours <- colorRampPalette(plot.colours)(length(samples))
      {pdf(width = 11, height = 8, file = paste0(filename, "_HW" , win, "bp_ylim", round(ylim2,2),"_Mirror=", mirror.mode,"_Norm=", normalize, "_StrandSep=", strand.sep, ".pdf"))
        plot(data[[1]][[1]]$Pos, yaxs = "i", xaxs = "i", ylim = c(ylim1, ylim2), xlim = c(-plot.w/2,plot.w/2), type = "n", ylab = paste0("Aggregated CC-seq signal (",win, " bp smoothed) / HpM"), xlab = "Position relative to feature midpoint/ bp")

        for (i in samples) {
          for(q in 1:strat.levels){
            lines(data[[i]][[q]]$Pos, data[[i]][[q]]$Smooth, lwd = 2, type = "l", col = plot.colours[match(i, samples)])
            text((-plot.w/2), ylim2-((match(i, samples))*(ylim2/20)), labels = DSBListNames[i], col = plot.colours[match(i, samples)], pos = 4)
          }

        }
        dev.off()}
    }
  }


  if(strand.sep){
    if(positive.strands){
    function2 <- function(x) {
      max(max(x$Top.Smooth[!is.na(x$Top.Smooth)]), max(x$Bottom.Smooth[!is.na(x$Top.Smooth)]))
    }
    ylim2<-1.1*(max(unlist(recurse(data, function2))))

    function2 <- function(x) {
      min(min2(x$Top.Smooth[!is.na(x$Top.Smooth)]), min2(x$Bottom.Smooth[!is.na(x$Top.Smooth)]))
    }
    ylim1<-0.9*(min2(unlist(recurse(data, function2))))
        if (ylims[1] != "auto"){
          if(length(ylims) == 2) {ylim1 <- ylims[1]; ylim2 <- ylims[2]}
          if(length(ylims) == 1) {ylim1 <- ylims[1]; ylim2 <- ylims[1]}
          }

    if(!combine){
      plot.colours <- rev(RColorBrewer::brewer.pal(4, "Spectral")); myColor.name <- "Spectral"
      if(!missing(col)){plot.colours <- col}
      plot.colours <- colorRampPalette(plot.colours)(strat.levels)
      for (i in samples) {
        {pdf(width = 11, height = 8, file = paste0("PileUp_",DSBListNames[i], "_on_", loci.name, "_", plot.w/1000, "KB_HW" , win, "bp_ylim", round(ylim2,2),"_Mirror=", mirror.mode,"_Norm=", normalize, "_StrandSep=", strand.sep,".pdf"))

          plot(data[[i]][[1]]$Pos, yaxs = "i", xaxs = "i", ylim = c(ylim1, ylim2), xlim = c(-plot.w/2,plot.w/2), type = "n", ylab = paste0("Aggregated CC-seq signal (",win, " bp smoothed) / HpM"), xlab = "Position relative to feature midpoint/ bp")

          text((plot.w/2), ylim2-(ylim2/20), labels = "Top", col = rgb(1,0,0,0.5), pos = 2)
          text((plot.w/2), (ylim2-2*(ylim2/20)), labels = "Bottom", col = rgb(0,0,1,0.5), pos = 2)
          for(q in 1:strat.levels){
            lines(data[[i]][[q]]$Pos, data[[i]][[q]]$Top.Smooth, lwd = 2, type = "l", col = rgb(1,0,0,0.5))
            lines(data[[i]][[q]]$Pos, data[[i]][[q]]$Bottom.Smooth, lwd = 2, type = "l", col = rgb(0,0,1,0.5))
            if(strat.levels > 1){
              if(strat.mode == "manual") text((-plot.w/2), ylim2-(q*(ylim2/20)), labels = paste0(ordinal(q)," Quartile of Value"), col = plot.colours[q], pos = 4)
              if(strat.mode == "expression") text((-plot.w/2), ylim2-(q*(ylim2/20)), labels = paste0(ordinal(q)," Quartile of Expression"), col = plot.colours[q], pos = 4)
              if(strat.mode == "length") text((-plot.w/2), ylim2-(q*(ylim2/20)), labels = paste0(ordinal(q)," Quartile of Length"), col = plot.colours[q], pos = 4)
            } else {
              text((-plot.w/2), ylim2-(q*(ylim2/20)), labels = DSBListNames[i], col = plot.colours[q], pos = 4)
            }
          }
          dev.off()}
      }
    }


    if(combine){
      plot.colours <- rev(RColorBrewer::brewer.pal(4, "Spectral")); myColor.name <- "Spectral"
      if(!missing(col)){plot.colours <- col}
      plot.colours <- colorRampPalette(plot.colours)(length(samples))
      {pdf(width = 11, height = 8, file = paste0(filename, "_HW" , win, "bp_ylim", round(ylim2,2),"_Mirror=", mirror.mode,"_Norm=", normalize, "_StrandSep=", strand.sep,".pdf"))
        plot(data[[1]][[1]]$Pos, yaxs = "i", xaxs = "i", ylim = c(ylim1, ylim2), xlim = c(-plot.w/2,plot.w/2), type = "n", ylab = paste0("Aggregated CC-seq signal (",win, " bp smoothed) / HpM"), xlab = "Position relative to feature midpoint/ bp")
        abline(h = 0, col = "grey")
        text((plot.w/2), ylim2-(ylim2/20), labels = "Top", col = rgb(1,0,0,0.5), pos = 2)
        text((plot.w/2), (ylim2-2*(ylim2/20)), labels = "Bottom", col = rgb(0,0,1,0.5), pos = 2)
        for (i in samples) {
          for(q in 1:strat.levels){
            lines(data[[i]][[q]]$Pos, data[[i]][[q]]$Top.Smooth, lwd = 2, type = "l", col = plot.colours[match(i, samples)])
            lines(data[[i]][[q]]$Pos, data[[i]][[q]]$Bottom.Smooth, lwd = 2, type = "l", col = plot.colours[match(i, samples)])
            text((-plot.w/2), ylim2-((match(i, samples))*(ylim2/20)), labels = DSBListNames[i], col = plot.colours[match(i, samples)], pos = 4)
          }

        }
        dev.off()}
    }
    } else {

    function2 <- function(x) {
      max(max(x$Top.Smooth[!is.na(x$Top.Smooth)]), max(x$Bottom.Smooth[!is.na(x$Top.Smooth)]))
    }
    ylim1<--1.1*(max(unlist(recurse(data, function2))))
    ylim2<-1.1*(max(unlist(recurse(data, function2))))
        if (ylims[1] != "auto"){
          if(length(ylims) == 2) {ylim1 <- ylims[1]; ylim2 <- ylims[2]}
          if(length(ylims) == 1) {ylim1 <- ylims[1]; ylim2 <- ylims[1]}
          }
    if(!combine){
      plot.colours <- rev(RColorBrewer::brewer.pal(4, "Spectral")); myColor.name <- "Spectral"
      if(!missing(col)){plot.colours <- col}
      plot.colours <- colorRampPalette(plot.colours)(strat.levels)
      for (i in samples) {
        {pdf(width = 11, height = 8, file = paste0("PileUp_",DSBListNames[i], "_on_", loci.name, "_", plot.w/1000, "KB_HW" , win, "bp_ylim", round(ylim2,2),"_Mirror=", mirror.mode,"_Norm=", normalize, "_StrandSep=", strand.sep,".pdf"))

          plot(data[[i]][[1]]$Pos, yaxs = "i", xaxs = "i", ylim = c(ylim1, ylim2), xlim = c(-plot.w/2,plot.w/2), type = "n", ylab = paste0("Aggregated CC-seq signal (",win, " bp smoothed) / HpM"), xlab = "Position relative to feature midpoint/ bp")
          abline(h = 0, col = "grey")
          text((plot.w/2), ylim2-(ylim2/20), labels = "Top", col = "black", pos = 2)
          text((plot.w/2), -(ylim2-(ylim2/20)), labels = "Bottom", col = "black", pos = 2)
          for(q in 1:strat.levels){
            lines(data[[i]][[q]]$Pos, data[[i]][[q]]$Top.Smooth, lwd = 2, type = "l", col = plot.colours[q])
            lines(data[[i]][[q]]$Pos, -data[[i]][[q]]$Bottom.Smooth, lwd = 2, type = "l", col = plot.colours[q])
            if(strat.levels > 1){
              if(strat.mode == "manual") text((-plot.w/2), ylim2-(q*(ylim2/20)), labels = paste0(ordinal(q)," Quartile of Value"), col = plot.colours[q], pos = 4)
              if(strat.mode == "expression") text((-plot.w/2), ylim2-(q*(ylim2/20)), labels = paste0(ordinal(q)," Quartile of Expression"), col = plot.colours[q], pos = 4)
              if(strat.mode == "length") text((-plot.w/2), ylim2-(q*(ylim2/20)), labels = paste0(ordinal(q)," Quartile of Length"), col = plot.colours[q], pos = 4)
            } else {
              text((-plot.w/2), ylim2-(q*(ylim2/20)), labels = DSBListNames[i], col = plot.colours[q], pos = 4)
            }
          }
          dev.off()}
      }
    }

    if(combine){
      plot.colours <- rev(RColorBrewer::brewer.pal(4, "Spectral")); myColor.name <- "Spectral"
      if(!missing(col)){plot.colours <- col}
      plot.colours <- colorRampPalette(plot.colours)(length(samples))
      {pdf(width = 11, height = 8, file = paste0(filename, "_HW" , win, "bp_ylim", round(ylim2,2),"_Mirror=", mirror.mode,"_Norm=", normalize, "_StrandSep=", strand.sep,".pdf"))
        plot(data[[1]][[1]]$Pos, yaxs = "i", xaxs = "i", ylim = c(ylim1, ylim2), xlim = c(-plot.w/2,plot.w/2), type = "n", ylab = paste0("Aggregated CC-seq signal (",win, " bp smoothed) / HpM"), xlab = "Position relative to feature midpoint/ bp")
        abline(h = 0, col = "grey")
        text((plot.w/2), ylim2-(ylim2/20), labels = "Top", col = "black", pos = 2)
        text((plot.w/2), -(ylim2-(ylim2/20)), labels = "Bottom", col = "black", pos = 2)
        for (i in samples) {
          for(q in 1:strat.levels){
            lines(data[[i]][[q]]$Pos, data[[i]][[q]]$Top.Smooth, lwd = 2, type = "l", col = plot.colours[match(i, samples)])
            lines(data[[i]][[q]]$Pos, -data[[i]][[q]]$Bottom.Smooth, lwd = 2, type = "l", col = plot.colours[match(i, samples)])
            text((-plot.w/2), ylim2-(i*(ylim2/20)), labels = DSBListNames[i], col = plot.colours[match(i, samples)], pos = 4)
          }

        }
        dev.off()}
    }
    }
  }
  # ENSURE THAT COMBINE IS ONLY ALLOWED IF STRAT.LEVELS == 1
  # lines(rep(0,2), c(0,20), lty = 2)
  # lines(c(-10000,10000), rep(expected.density,2), lty = 3)

}
