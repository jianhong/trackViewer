getHeight <- function(SNPs, ratio.yx, LINEW, GAP, cex, type,
                      level=c("data", "data&labels")){
    level=match.arg(level)
    stack.factors <- unique(as.character(SNPs$stack.factor))
    stack.factors <- sort(stack.factors)
    if(length(SNPs$score)>0) {
        scoreMax <- ceiling(max(c(SNPs$score, 1), na.rm=TRUE))
    }else {
        scoreMax <- 1
    }
    if(type=="pie.stack") scoreMax <- length(stack.factors)
    if(!type %in% c("pie", "pie.stack")){
        if(scoreMax>10) {
            SNPs$score <- 10*SNPs$score/scoreMax 
            scoreMax <- ceiling(max(c(SNPs$score, 1), na.rm=TRUE))
        }
    }
    if(level=="data"){
        switch(type,
               circle={
                   labels.y <- LINEW*max(ratio.yx, 1.2) + 
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
                   labels.y <- LINEW*max(ratio.yx, 1.2) + 
                       6.5*GAP*cex + 
                       (this.scores-0.5) * LINEW * ratio.yx*cex
               },
               pie={
                   labels.y <- LINEW*max(ratio.yx, 1.2) + 
                       6.5*GAP*cex + 0.5 * LINEW * ratio.yx * cex
               },
               pie.stack={
                   labels.y <- LINEW*max(ratio.yx, 1.2) + 
                       6.5*GAP*cex + 
                       (scoreMax-0.5) * LINEW * ratio.yx*cex
               })
        labels.y
    }else{
        if(length(SNPs$label.parameter.rot)>0) {
            labels.rot <- SNPs$label.parameter.rot
        }else{
            labels.rot <- 90
        }
        labels.length.rate <- max(cospi((labels.rot-90)/180), 0)
        
        switch(type,
               circle={
                   if(length(names(SNPs))>0){
                       maxStrHeight <- 
                           max(as.numeric(
                               convertY(stringWidth(names(SNPs)), "npc")
                           ))+LINEW/2
                   }else{
                       maxStrHeight <- 0
                   }
                   maxStrHeight <- maxStrHeight * labels.length.rate
                   ypos <- LINEW*max(ratio.yx, 1.2) + 6.5*GAP*cex + 
                       (scoreMax-0.5) * LINEW * ratio.yx*cex + maxStrHeight*cex
               },
               pin={
                   if(length(names(SNPs))>0){
                       thisStrHeight <- as.numeric(
                           convertY(stringWidth(names(SNPs)), "npc"))+LINEW/2
                   }else{
                       thisStrHeight <- 0
                   }
                   thisStrHeight <- thisStrHeight * labels.length.rate
                   if(length(SNPs$score)>0){
                       ypos <- 
                           max(LINEW*max(ratio.yx, 1.2) + 
                                   6.5*GAP*cex + 
                                   (SNPs$score-0.5) * LINEW * ratio.yx*cex + 
                                   thisStrHeight*cex)
                   }else{
                       ypos <- max(LINEW*max(ratio.yx, 1.2) + 
                                       6.5*GAP*cex + thisStrHeight*cex)
                   }
               },
               pie={
                   if(length(names(SNPs))>0){
                       maxStrHeight <- 
                           max(as.numeric(
                               convertY(stringWidth(names(SNPs)), "npc")
                           ))+LINEW/2
                   }else{
                       maxStrHeight <- 0
                   }
                   maxStrHeight <- maxStrHeight * labels.length.rate
                   ypos <- LINEW*max(ratio.yx, 1.2) + 
                       6.5*GAP*cex + maxStrHeight*cex
               },
               pie.stack={
                   if(length(names(SNPs))>0){
                       maxStrHeight <- 
                           max(as.numeric(
                               convertY(stringWidth(names(SNPs)), "npc")
                           ))+LINEW/2
                   }else{
                       maxStrHeight <- 0
                   }
                   maxStrHeight <- maxStrHeight * labels.length.rate
                   ypos <- LINEW*max(ratio.yx, 1.2) + 
                       6.5*GAP*cex + maxStrHeight*cex +
                       (scoreMax-0.5) * LINEW * ratio.yx*cex
               }
        )
        ypos
    }
}
