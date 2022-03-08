#' @title featureplotter
#' @description Wrapper for Matt's gene track plotter in the original mapper.
#' @param x Featurelist subsetted from the AllElementsDUB
#' @param col The fill colour for the feature arrows.
#' @param lty Line type for border. 1 for solid, 2 for dotted.
#' @param strand Which strand are the features on? can be "+", "-" or "*"
#' @param av Vertical location at which to plot feature.
#' @param tv Vertical location at which to plot text. Default is same as av (plots within feature).
#' @examples
#' @author Matt Neale and Will Gittens
#' @export
featureplotter <- function(x, coln="genename",col = "wheat", lty = 1, strand = "+", av = 50,lab = T, tv = av, pos1 = xl1, pos2 = xl2,cex.size=0.9) {
  xrange <- pos2-pos1
  ahead<-xrange/25 #make arrowhead length proportional to plot range
  ahead[(ahead>500)]<-500 #limit max length to 500
  ahw<-15 #arrowhead width
    if((pos2-pos1)>=100000){ lab = F}
    if(nrow(x) > 0){
    if(strand == "+") {
      for(i in 1:nrow(x)){
        if(ahead < (x[i,"stop"] - x[i,"start"])){
          polygon(c(x[i,"start"], x[i,"stop"]-ahead, x[i,"stop"]-ahead,x[i,"stop"], x[i,"stop"]-ahead, x[i,"stop"]-ahead, x[i,"start"]),c(av+ahw,av+ahw,av+ahw+ahw,av,av-ahw-ahw,av-ahw, av-ahw), col=col, border =NA, lty = lty)
          if(lab){text((x[i,"start"]+x[i,"stop"])/2,tv, font=3, x[i,coln], cex=cex.size)
           if(x[i,coln] == "") {text((x[i,"start"]+x[i,"stop"])/2,tv, font=3, x[i,"sysname"], cex=cex.size)}
}
        } else {
          polygon(c(x[i,"start"], x[i,"stop"], x[i,"start"]),c(av+ahw+ahw,av,av-ahw-ahw), col=col, border ="snow4", lty = lty)
         if(lab){ text((x[i,"start"]+x[i,"stop"])/2,tv, font=3, x[i,coln], cex=cex.size)
           if(x[i,coln] == "") {text((x[i,"start"]+x[i,"stop"])/2,tv, font=3, x[i,"sysname"], cex=cex.size)}
           }
        }
      }
    }
    if(strand == "-") {
      for(i in 1:nrow(x)){
        if(ahead < (x[i,"stop"] - x[i,"start"])){
          polygon(c(x[i,"stop"], x[i,"start"]+ahead, x[i,"start"]+ahead,x[i,"start"], x[i,"start"]+ahead, x[i,"start"]+ahead, x[i,"stop"]),c(av+ahw,av+ahw,av+ahw+ahw,av,av-ahw-ahw,av-ahw, av-ahw), col=col, border =NA, lty = lty)
          if(lab){  text((x[i,"start"]+x[i,"stop"])/2,tv, font=3, x[i,coln], cex=cex.size)
          if(x[i,coln] == "") {text((x[i,"start"]+x[i,"stop"])/2,tv, font=3, x[i,"sysname"], cex=cex.size)}
}
        } else {
          polygon(c(x[i,"stop"], x[i,"start"], x[i,"stop"]),c(av+ahw+ahw,av,av-ahw-ahw), col=col, border =NA, lty = lty)
          if(lab){text((x[i,"start"]+x[i,"stop"])/2,tv, font=3, x[i,coln], cex=cex.size)
            if(x[i,coln] == "") {text((x[i,"start"]+x[i,"stop"])/2,tv, font=3, x[i,"sysname"], cex=cex.size)}
            }
        }
      }
    }
    if(strand == "*") {
      for(i in 1:nrow(x)){
        polygon(c(x[i,"start"],x[i,"stop"], x[i,"stop"], x[i,"start"]),c(av+ahw,av+ahw,av-ahw, av-ahw), col=col, border ="black", lty = lty)
         text((x[i,"start"]+x[i,"stop"])/2,tv, font=3, x[i,coln], cex=cex.size)
          if(x[i,coln] == "") {text((x[i,"start"]+x[i,"stop"])/2,tv, font=3, x[i,"sysname"], cex=cex.size)}

      }
    }
  }
}
