# How to use some of Will's core CC-seq functions.

# Author: Will Gittens 2021

# Here I provide use-able functions,  with commonly-used arguments explained as written in the documentation.
# If you identify problems with a function, please can you raise and issue on GitHub at https://github.com/gb305/Rpackages/issues
# Please put the function name in the Issue title, and include a copy-paste of the failing function call. Thanks.

#install and load packages
install.packages("devtools")
install.packages("git2r")
library("devtools")
library("git2r")
library(CCTools)


f.dir <- "PATH/TO/FOLDER/CONTAINING/FULLMAPS"


# CCBin
CCBin(path.name = f.dir, # A character string defining the RDS file or folder containing FullMaps.
      exp.name = "exp1", # Manually specify a name for the experiment.
      genome = "Cer3H4L2",
      bin.width = 100,
      os = 3, # Do you want to offset the Crick strand relative to the Watson strand prior to binning? Choose 3 for Top2, or 1 for Spo11. Choose 0 for no ofsetting. Defaults to 3.
      out.mode = "full", # "sparse" or "full" format. Defaults to "full". Defaults to "full"
      out.mode2 = "strand", # Do you want to sum the Watson and Crick strands ("total") or report them seperately ("strand")? Defaults to "strand".
      plot.mode = "png", # 	Do you want to plot the binned chromosome maps as "pdf" of "png"? Choose "none" for no plots
      ylims = 4) #  A single number specifying the maximum limit of the y axis, for both Watson and Crick. OR "auto", to automatically define for the entire set.

?CCBin # run this to see documentation in righthand Help box, with details on further options that you may want to use.

# CCMap

CCMap(path.name = f.dir, # A character string defining the RDS file or folder containing FullMaps.
      exp.name = "exp1", # Manually specify a name for the experiment.
      genome = "Cer3H4L2",
      window.w = 500, # How wide a region to plot in bp?
      loci = "manual", # Do you want to manually specify a region using Chr and Pos ("manual"), a specific single gene ("GENENAME"), or multiple maps centred on a types of features in the AllElementsDUB table (e.g "gene", "CEN", "ARS_consensus_sequence". For more options run CCFeatures() ). You can also provide a vector where the second element is "start", "mid", or "end", to centre the pileup on , for example, the end of a gene.
      chrom. = 3, # only used if loci = "manual"
      pos = 2439550, # only used if loci = "manual"
      os = 0, # Do you want to offset the Crick strand relative to the Watson strand prior to binning? Choose 3 for Top2, or 1 for Spo11. Choose 0 for no ofsetting. Defaults to 3.
      samples = "all", # Do you just want to plot specific samples from the folder? If so, provide a numeric vector of the sample indexes, assuming that they are sorted by name (like in Finder). Usually easier to just remove fullmaps from the folder.
      ylims = 50, #  A single number specifying the maximum limit of the y axis, for both Watson and Crick. OR "auto", to automatically define for each sample.
      smooth.mode = "none", # This option controls whether to apply VariableX ("VX") or Hann window ("HW") smoothing to the broad overlay. Deafults to no smoothing.
      Ty.track = F, # Plot Ty elements?
      gene.track = T, # Plot gene track?
      ARS.track = T, # Plot plot ARS consensus sequences?
      MNase.track = F,
      scc1.track = F,
      phase = "vegetative", gene.names = "genename") # "meiosis" or "vegetative". This changes some of the annotation files (e.g MNase data), but is not fully implemented.

?CCMap # run this to see documentation in righthand Help box, with details on further options that you may want to use.

# CCPileup

#IMPORTANT NOTE: In CCPileup, the script will first look locally for a previously generated pileup which matches the exp.name and loci specified in the function call.
# This is to save time, because the pileup takes a long time, whereas the plotting does not. Therefore, if you simply wish to replot after changing plotting parameters
# (e.g zoom in, or change the smoothing, strand split, etc), the script should intelligently realise this. NOTE THOUGH, that it does this based only on a unique combination
# of exp.name, loci, strat.mode, strat.levels and window.w. If you use the same exp.name for different experiments, you will run into problems. I will figure out a better implementation soon.
# This feature can be turned off with cache.search = F.

CCPileup(path.name = f.dir,  # A character string defining the RDS file or folder containing FullMaps.
         exp.name = "Pombe1", # Manually specify a name for the experiment.
         phase = "vegetative",
         genome = "Cer3H4L2",
         loci =  "CEN", # he path to an RDS file or .txt file containing a table of loci in column format: c("Chr", "Pos", "Strand"). Alternatively can use inbuilt loci such as "TSS" (From AllElementsDUB), or "CTSS" (From YeastTSS), "CEN" (Centromeres), ARS (ARS_consensus_sequence excluding Chr12, to avoid rDNA). More info on the generation of these sets is available from Will. Note you can't use any feature type from the AllElementsDUB, like you can with CCMap. This will be implemented soon.
         window.w = 100000, # How wide a region to PileUp in bp?
         plot.w = 100000, # How wide a region to plot in bp? This is distinct from window.w, to avoid repeating the time-consuming pileup unless necessary. If you plot a smaller window (e.g zoom in), it just uses the same PileUp
         combine = T, # Do you want to plot all samples on the same plot? If so, it is recommended to set strat.mode to "none", or the plots will be very busy. Defaults to FALSE
         samples = "all", # Do you just want to plot specific samples from the folder? If so, provide a numeric vector of the sample indexes, assuming that they are sorted by name (like in Finder).
         strand.sep = F, # Do you want to seperate the Watson and Crick strand?
         positive.strands = T, #If you do seperate the Watson and Crick strands, do you want to plot them both on the same positive axis?
         normalize = F, # Do you want to normalize every pileup over the range of window.w. (NOTE: IT NORMALISES OVER window.w , not plot.w!)
         cache.search = F, # Do you want to search in the folders surrounding the path.name for a matching saved PileUp? Defaults to true.
         mirror.mode = F, # Do you want to average the plot either side of the loci? This generates a symmetrical plot, which is useful when considering loci with no inherent directionality (e.g centromeres). Defaults to FALSE.
         strat.mode = "none", # How do you want to stratify the plots? "none" will leave the Pileups unstratified. "manual" will point the stratifier to a 4th numeric column (column name must be "Value") in the user-defined loci table. "expression" or "length" can be used in combination with the inbuilt TSS or CTSS loci tables, and refer to microarray gene expression, and gene length, respectively.
         strat.levels = 4, # If you are straifying, how many strata do you want to use? Defaults to 4.
         ylims = "auto", #  A single number specifying the maximum limit of the y axis, for both Watson and Crick. OR "auto", to automatically define for each sample.
         win = 5000) # Hanning window width for smoothing. Change to 1 for no smoothing.


?CCPileup # run this to see documentation in righthand Help box, with details on further options that you may want to use.


(path.name = f.dir, # A character string defining the RDS file or folder containing FullMaps.
           rDNA = T, #do you want to include the rDNA?
           chroms = 1:3, # which chromosomes to use?
           genome = "Cer3H4L2", #Which genome?
           dyad = 1.5, # where do you man to plot the dyad axis guidelines? 1.5 for Top2, 0.5 for Spo11.
           minThresholdv = c(0, 1, 5),
           maxThreshold = "MAX",)


# CCSeqbias2

# Aggregates DNA base composition around CC-seq hits. This version is Will's work in progress, to make a faster, and more generic script.

CCSeqbias2(path.name = f.dir, # A character string defining the RDS file or folder containing FullMaps.
minThresholdv = c(0, 1, 5), # vector of lower threshold limits
maxThreshold = "MAX", # a vector of upper threshold limits. Or can use "MAX"
w = 20, #Window width in bp to aggregate over
rDNA = F,  #do you want to include the rDNA?
chroms = c(1:16),  # which chromosomes to use?
ylims = c(0.05, 0.6), # ylims
dyad = 1.5,  # where do you man to plot the dyad axis guidelines? 1.5 for Top2, 0.5 for Spo11.
genome = "Cer3H4L2", #Which genome?
plotT = "o") #PlotType: "l", "p", "o" (lines, points or both)


?CCSeqbias2


