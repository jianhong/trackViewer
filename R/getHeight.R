getHeight <- function(SNPs, ratio.yx, LINEW, GAP, cex, type, scoreMax,
                      level=c("data", "data&labels")){
    level=match.arg(level)
    stack.factors <- unique(as.character(SNPs$stack.factor))
    stack.factors <- sort(stack.factors)
    if(level=="data"){
        switch(type,
               circle={
                   labels.y <- LINEW + # add gaps for labels
                       6.5*GAP*cex + 
                       (scoreMax-0.5) * LINEW * ratio.yx*cex
               },
               pin={
                   if(length(SNPs$score)>0) {
                       this.scores <- ceiling(SNPs$score)
                   }else {
                       this.scores <- .5
                   }
                   this.scores[is.na(this.scores)] <- .5
                   labels.y <- LINEW + 
                       6.5*GAP*cex + 
                       (this.scores-0.5) * LINEW * ratio.yx*cex
               },
               pie={
                   labels.y <- LINEW*max(ratio.yx, 1.2) + 
                       6.5*GAP*cex + 0.5 * LINEW * ratio.yx * cex
               },
               pie.stack={
                   labels.y <- LINEW + 
                       6.5*GAP*cex + 
                       (scoreMax-0.5) * LINEW * ratio.yx*cex
               },
               flag={
                 labels.y <- LINEW + 
                   6.5*GAP*cex + 
                   scoreMax * LINEW * ratio.yx*cex
               })
        labels.y
    }else{
        if(length(SNPs$label.parameter.rot)>0) {
            labels.rot <- SNPs$label.parameter.rot
        }else{
            labels.rot <- 90
        }
        labels.cex <- 1
        if(length(SNPs$label.parameter.gp)>0){
          if(length(SNPs$label.parameter.gp$cex)>0)
            labels.cex <- SNPs$label.parameter.gp$cex[[1]][1]
        }
        labels.length.rate <- labels.cex * max(cospi((labels.rot-90)/180), 0) * ratio.yx
        stringH <- as.numeric(convertY(stringHeight("W"), "npc"))
        
        switch(type,
               circle={
                   if(length(names(SNPs))>0){
                       maxStrHeight <- 
                           max(as.numeric(
                               convertX(stringWidth(names(SNPs)), "npc")
                           ))+stringH
                   }else{
                       maxStrHeight <- 0
                   }
                   maxStrHeight <- maxStrHeight * labels.length.rate
                   ypos <- LINEW + 6.5*GAP*cex + 
                       (scoreMax-0.5) * LINEW * ratio.yx*cex + maxStrHeight
               },
               pin={
                   if(length(names(SNPs))>0){
                       thisStrHeight <- max(as.numeric(
                           convertX(stringWidth(names(SNPs)), "npc")) ) +
                         stringH
                   }else{
                       thisStrHeight <- 0
                   }
                   thisStrHeight <- thisStrHeight * labels.length.rate
                   if(length(SNPs$score)>0){
                       ypos <- 
                           max(LINEW + 
                                   6.5*GAP*cex + 
                                   (SNPs$score-0.5) * LINEW * ratio.yx*cex + 
                                   thisStrHeight)
                   }else{
                       ypos <- max(LINEW*max(ratio.yx, 1.2) + 
                                       6.5*GAP*cex + thisStrHeight)
                   }
               },
               pie={
                   if(length(names(SNPs))>0){
                       maxStrHeight <- 
                         max(as.numeric(
                           convertX(stringWidth(names(SNPs)), "npc")
                         ))+stringH
                   }else{
                       maxStrHeight <- 0
                   }
                   maxStrHeight <- maxStrHeight * labels.length.rate
                   ypos <- LINEW + 
                       6.5*GAP*cex + maxStrHeight
               },
               pie.stack={
                   if(length(names(SNPs))>0){
                       maxStrHeight <- 
                         max(as.numeric(
                           convertX(stringWidth(names(SNPs)), "npc")
                         ))+stringH
                   }else{
                       maxStrHeight <- 0
                   }
                   maxStrHeight <- maxStrHeight * labels.length.rate
                   ypos <- LINEW + 
                       6.5*GAP*cex + maxStrHeight +
                       (scoreMax-0.5) * LINEW * ratio.yx*cex
               },
               flag={
                 if(length(names(SNPs))>0){
                   maxStrHeight <- 
                     max(as.numeric(
                       convertX(stringWidth(names(SNPs)), "npc")
                     ))+stringH
                 }else{
                   maxStrHeight <- 0
                 }
                 maxStrHeight <- maxStrHeight * labels.length.rate
                 ypos <- LINEW + 6.5*GAP*cex + 
                   scoreMax * LINEW * ratio.yx*cex + maxStrHeight
               }
        )
        ypos
    }
}
