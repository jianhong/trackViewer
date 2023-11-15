plotLollipopData <- function(dat, xscale, chr, yaxis, yaxis.gp, ybase=0, side="top", main=TRUE, baselineCol="black", maxHeight, lollipop_style_switch_limit=10){
  if(length(dat)>0){
    dat <- orderedGR(dat)
    xscale.gr <- GRanges(seqnames = chr, ranges=IRanges(min(xscale), max(xscale)))
    dat <- subsetByOverlaps(dat, xscale.gr, ignore.strand=TRUE)
    if(length(dat)<1) return(invisible())
    width(dat) <- 1
    TYPES <- c("circle", "pie", "pin", "pie.stack", "flag")
    type <- if(is.list(dat$type)) dat$type[[1]] else dat$type[1]
    if(length(type)==0) type <- "circle"
    if(!type %in% TYPES) type <- "circle"
    if(type=="pin"){ ## read the pin shape file
      pinpath <- system.file("extdata", "map-pin-red.xml", package="trackViewer")
      pin <- readPicture(pinpath)
    }else{
      pin <- NULL
    }
    cex <- if(is.list(dat$cex)) dat$cex[[1]] else dat$cex[1]
    if(length(cex)==0) cex <- 1
    dashline.col <- if(is.list(dat$dashline.col)) dat$dashline.col[[1]] else dat$dashline.col[1]
    if(length(dashline.col)==0) dashline.col <- "gray80"
    jitter <- if(is.list(dat$jitter)) dat$jitter[[1]] else dat$jitter[1]
    if(length(jitter)==0) jitter <- "node"
    if(!jitter %in% c("node", "label")) jitter <- "node"
    scoreMax0 <- scoreMax <- 
      if(length(dat$score)>0) ceiling(max(c(dat$score, 1), na.rm=TRUE)) else 1
    if(type=="pie.stack") scoreMax <- length(unique(dat$stack.factor))
    if(!type %in% c("pie", "pie.stack")){
      if(scoreMax>lollipop_style_switch_limit) {
        dat$score <- 10*dat$score/scoreMax
        scoreMax <- 10*scoreMax0/scoreMax
      }else{
        scoreMax <- scoreMax0
      }
      scoreType <- 
        if(length(dat$score)>0) all(floor(dat$score)==dat$score) else FALSE
    }else{
      scoreType <- FALSE
    }
    LINEW <- as.numeric(convertX(unit(1, "line"), "npc"))/2
    LINEH <- as.numeric(convertY(unit(1, "line"), "npc"))/2
    ## GAP the gaps between any elements
    GAP <- .2 * LINEH
    ratio.yx <- 1/as.numeric(convertX(unit(1, "snpc"), "npc"))
    feature.height <- 
      if(is.list(dat$feature.height)) dat$feature.height[[1]] else dat$feature.height[1]
    if(length(feature.height)==0) feature.height <- GAP
    if(maxHeight>1) cex <- cex/maxHeight
    scoreMax0 <- scoreMax0/cex
    guideline <- ifelse(side=="top", ybase, ybase+feature.height)
    grid.lines(y=c(guideline, guideline), gp = gpar(col=baselineCol))
    plotLollipops(dat, feature.height=feature.height, bottomHeight=ybase, baseline=0, 
                  type=type, ranges=xscale.gr, yaxis=yaxis, yaxis.gp=yaxis.gp,
                  scoreMax=scoreMax, scoreMax0=scoreMax0, scoreType=scoreType, 
                  LINEW, cex, ratio.yx, GAP, pin, dashline.col,
                  side=side, jitter=jitter, main=main)
  }
  return(invisible())
}