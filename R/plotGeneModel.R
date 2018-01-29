##levels of feature "CDS"   "ncRNA" "utr3"  "utr5"
plotGeneModel <- function(track, xscale, chr){
  if(length(track@dat2)>0){
    feature.height <- 
      if(is.list(track@dat2$feature.height)) track@dat2$feature.height[[1]] else track@dat2$feature.height[1]
    if(length(feature.height)==0) feature.height <- .5
    unit <- feature.height/4
    y <- feature.height/2
  }else{
    unit <- 0.25
    y <- 0.5
  }
    transcript <- track@dat
    col <- track@style@color
    strand <- as.character(strand(transcript))[1]
    unitx <- convertWidth(unit(3, "lines"), 
                          unitTo = "npc",
                          valueOnly = TRUE)
    if(unitx > 0.5) unitx <- 0.5
    unitx <- unitx * abs(xscale[2] - xscale[1]) #abs(xscale[2] - xscale[1]) / 10
    ## plot introns
    introns <- gaps(ranges(transcript))
    introns <- introns[start(introns)>min(end(transcript))
                       &end(introns)<max(start(transcript))]
    plotArrow <- function(ax0, ax1, y, col, vp){
        ay0 <- y + unit/2
        ay1 <- y
        ay2 <- y - unit/2
        grid.segments(ax0, ay0, ax1, ay1, 
                      default.units="native", 
                      gp=gpar(col=col))
        grid.segments(ax1, ay1, ax0, ay2, 
                      default.units="native", 
                      gp=gpar(col=col))
    }
    plotIntrons <- function(strand, a, b, y, xscale, col){
        if(a<xscale[1]) a <- xscale[1]
        if(b>xscale[2]) b <- xscale[2]
        w <- b - a
        times <- floor(5*w/unitx)
        for(i in seq_len(times)){
            if(b > a+i*unitx/5+unitx*3/20){
                if(strand=="+") plotArrow(a+i*unitx/5+unitx/10, 
                                          a+i*unitx/5+unitx*3/20,
                                          y, col)
                else {
                    if(strand=="-") 
                        plotArrow(a+i*unitx/5+unitx*3/20, 
                                  a+i*unitx/5+unitx/10,
                                  y, col)
                }
            }
        }
    }
    if(length(introns)>0){
        for(i in 1:length(introns)){
            if(start(introns)[i]<=xscale[2] && end(introns)[i]>=xscale[1]){
                grid.segments(start(introns)[i], y, 
                              end(introns)[i], y, default.units="native", 
                              gp=gpar(col=col))
                plotIntrons(strand, start(introns[i]), 
                            end(introns[i]),
                            y, xscale, col)
            }
        }
    }

    ## plot utr
    utr <- transcript[transcript$feature %in% c("utr3", "utr5")]
    if(length(utr)>0){
        for(i in 1:length(utr)){
            if(start(utr)[i]<=xscale[2] && end(utr)[i]>=xscale[1]){
                grid.rect(start(utr)[i], y,
                          end(utr)[i]-start(utr)[i]+1, unit, 
                          default.units="native", just=c(0, 0.5),
                          gp=gpar(col=NA, fill=col))
            }
        }
    }
    
    ## plot exons
    exons <- transcript[transcript$feature %in% c("CDS", "ncRNA", "exon")]
    if(length(exons)>0){
        for(i in 1:length(exons)){
            if(start(exons)[i]<=xscale[2] && end(exons)[i]>=xscale[1]){
                grid.rect(start(exons)[i], 
                          y, 
                          end(exons)[i]-start(exons)[i]+1, 
                          2*unit, default.units="native", just=c(0, 0.5),
                          gp=gpar(col=NA, fill=col))
            }
        }
    }
    
    if(length(track@dat2)>0){
      track@dat2 <- orderedGR(track@dat2)
      xscale.gr <- GRanges(seqnames = chr, ranges=IRanges(min(xscale), max(xscale)))
      track@dat2 <- subsetByOverlaps(track@dat2, xscale.gr, ignore.strand=TRUE)
      if(length(track@dat2)<1) return(invisible())
      width(track@dat2) <- 1
      TYPES <- c("circle", "pie", "pin", "pie.stack")
      type <- if(is.list(track@dat2$type)) track@dat2$type[[1]] else track@dat2$type[1]
      if(length(type)==0) type <- "circle"
      if(!type %in% TYPES) type <- "circle"
      if(type=="pin"){ ## read the pin shape file
        pinpath <- system.file("extdata", "map-pin-red.xml", package="trackViewer")
        pin <- readPicture(pinpath)
      }else{
        pin <- NULL
      }
      cex <- if(is.list(track@dat2$cex)) track@dat2$cex[[1]] else track@dat2$cex[1]
      if(length(cex)==0) cex <- 1
      dashline.col <- if(is.list(track@dat2$dashline.col)) track@dat2$dashline.col[[1]] else track@dat2$dashline.col[1]
      if(length(dashline.col)==0) dashline.col <- "gray80"
      jitter <- if(is.list(track@dat2$jitter)) track@dat2$jitter[[1]] else track@dat2$jitter[1]
      if(length(jitter)==0) jitter <- "node"
      if(!jitter %in% c("node", "label")) jitter <- "node"
      scoreMax0 <- scoreMax <- 
        if(length(track@dat2$score)>0) ceiling(max(c(track@dat2$score, 1), na.rm=TRUE)) else 1
      if(type=="pie.stack") scoreMax <- length(unique(track@dat2$stack.factor))
      if(!type %in% c("pie", "pie.stack")){
        if(scoreMax>10) {
          track@dat2$score <- 10*track@dat2$score/scoreMax
          scoreMax <- 10*scoreMax0/scoreMax
        }else{
          scoreMax <- scoreMax0
        }
        scoreType <- 
          if(length(track@dat2$score)>0) all(floor(track@dat2$score)==track@dat2$score) else FALSE
      }else{
        scoreType <- FALSE
      }
      LINEW <- as.numeric(convertX(unit(1, "line"), "npc"))
      LINEH <- as.numeric(convertY(unit(1, "line"), "npc"))
      ## GAP the gaps between any elements
      GAP <- .2 * LINEH
      ratio.yx <- 1/as.numeric(convertX(unit(1, "snpc"), "npc"))
      plotLollipops(track@dat2, feature.height=y+unit, bottomHeight=0, baseline=y, 
                    type=type, ranges=xscale.gr, yaxis=FALSE, 
                    scoreMax=scoreMax, scoreMax0=scoreMax0, scoreType=scoreType, 
                    LINEW, cex, ratio.yx, GAP, pin, dashline.col,
                    side="top", jitter=jitter)
    }
    return(invisible())
}