xysmooth <- function(x2, y2, smooth=5){
  if(is.logical(smooth)) smooth=5
  if(smooth>=length(x2)){
    return(list(x=x2, y=y2))
  }
  ord <- order(x2)
  x2 <- x2[ord]
  y2 <- y2[ord]
  x <- min(x2):max(x2)
  k <- ksmooth(x2, y2, n.points = length(x), x.points = x, bandwidth = 1)
  y <- k$y
  y.num <- which(!is.na(y))
  y[is.na(y)] <- unlist(mapply(function(from, to, length.out){
    s <- seq(from, to, length.out = length.out)
    s[-c(1, length(s))]
    }, y[y.num[-length(y.num)]], y[y.num[-1]], diff(y.num)+1, SIMPLIFY = FALSE))
  k <- ksmooth(x, y, n.points = length(x), x.points = x, bandwidth = smooth)
  id <- which(x %in% x2)
  x <- k$x[id]
  y <- k$y[id]
  id <- which(y2==0)
  id1 <- which(diff(id)==1)
  id1 <- sort(unique(c(id1, id1+1)))
  id <- id[id1]
  id1 <- which(diff(x[id])>(diff(range(x))+1)/100) ## resolution set to 100
  id1 <- sort(unique(c(id1, id1+1)))
  id <- id[id1]
  y[id] <- 0
  list(x=x, y=y)
}
plotDataTrack <- function(.dat, chr, strand, scale, color, yscale, smooth=FALSE){
    names(.dat) <- NULL
    mcols(.dat) <- mcols(.dat)[, "score"]
    colnames(mcols(.dat)) <- "score"
    if(missing(yscale)) yscale <- c(0, 1)
    if(length(.dat)>0){##subset if length .dat > 1000
        trackViewerCache <- list(step=1000)
        if(length(.dat)>trackViewerCache$step){
            .width <- which(width(.dat) > diff(scale)[1]/50)
            .idx <- sample(1:length(.dat), trackViewerCache$step)
            .idx <- unique(c(.idx, .width))
            .idx <- .idx[order(.idx)]
            .dat1 <- .dat[-.idx]
            .datM <- reduce(.dat1, with.revmap=TRUE)
            .map <- .datM$revmap
            .gp <- rep(1:length(.map), sapply(.map, length))
            .datM$score <- floor(sapply(split(.dat1$score[unlist(.map)], 
                                              .gp), mean))
            .datM$revmap <- NULL
            .dat <- c(.dat[.idx], .datM)
            .dat <- orderedGR(.dat)
        }
        .dat <- c(.dat, GRanges(seqnames=chr, 
                                IRanges(start=scale, end=scale), 
                                strand=strand, score=c(0,0)))
    }else{
        .dat <- GRanges(seqnames=chr, 
                        IRanges(start=scale, end=scale), 
                        strand=strand, score=c(0,0))
    }
    .data <- split(.dat, strand(.dat))
    if(strand!="*") .data <- .data[names(.data)==strand]
    xt <- c()
    yt <- c()
    for(i in 1:length(.data)){
        if(length(.data[[i]])>0){
            .dat <- .data[[i]]
            .dat <- condenceGRs(.dat) ## can not handle +/- strands in the same track
            .dat$type <- "D"
            .gap <- gaps(.dat) ## insert 0
            .gap <- .gap[start(.gap)>=scale[1] & end(.gap)<=scale[2]]
            if(length(.gap)>0){
                mcols(.gap)$score <- 0
                .gap$type <- "G"
                .dat <- c(.dat, .gap)
            }
            .dat <- orderedGR(.dat)
            .dat <- as.data.frame(.dat)
            ##expand the gap by 1
            .dat[.dat$type=="G", "start"] <- .dat[.dat$type=="G", "start"] - .5
            .dat[.dat$type=="G", "end"] <- .dat[.dat$type=="G", "end"] + .5
            x <- as.numeric(t(.dat[,c("start", "end")]))
            y <- as.numeric(rep(.dat[,"score"], each=2))
            x2 <- c(min(x), x, max(x))
            y2 <- c(0, y, 0)
            ## do smooth
            if(smooth){
              xy.smoothed <-xysmooth(x2, y2, smooth)
              x2 <- xy.smoothed$x
              y2 <- xy.smoothed$y
            }
            grid.polygon(x2, y2, default.units="native", 
                         gp=gpar(col=NA, fill=color))
            yscale <- range(yscale)
            if(yscale[1]<=0 && yscale[2]>=0){
              grid.lines(x=scale, y=0, default.units="native",
                       gp=gpar(col=color))
            }else{
              dy <- abs(yscale - 0)
              dy <- yscale[order(dy)][1]
              grid.lines(x=scale, y=dy, default.units="native",
                         gp=gpar(col="red"))
            }
            xt <- c(xt, x)
            yt <- c(yt, y)
        }
    }
    return(list(x=xt, y=yt))
}
plotTrack <- function(name, track, curViewStyle, curYpos,
                      yscale, height, xlim, chr, strand,
                      operator, wavyLine, smooth=FALSE){
    style <- track@style
    yHeightBottom <- yHeightTop <- 0.01
    
    if(curViewStyle@flip){ 
        xscale <- rev(xlim)
    }else{
        xscale <- xlim
    }
    if(length(style@marginBottom)>0){
        yHeightBottom <- style@marginBottom
    }
    if(length(style@marginTop)>0){
        yHeightTop <- style@marginTop
    }
    pushViewport(viewport(x=0, y=curYpos, 
                          height=height,
                          width=1, 
                          just=c(0,0))) ## vp1
    ##put ylab
    if(track@type %in% c("transcript", "gene") && style@ylabpos %in% c("upstream", "downstream")){
        if(length(findOverlaps(ranges(range(track@dat)), IRanges(xscale[1], xscale[2])))>0){
          putGeneYlab(curViewStyle, style, name, height, xscale,
                      range(track@dat), length(track@dat2)>0)
        } 
    }else{
        if(style@ylabpos=="upstream") style@ylabpos <- "left"
        if(style@ylabpos=="downstream") style@ylabpos <- "right"
        putYlab(curViewStyle, style, name, yHeightBottom, yHeightTop, height, yscale)
    }
    
    pushViewport(viewport(x=0, y=yHeightBottom, 
                          height=1 - yHeightTop - yHeightBottom,
                          width=1, 
                          just=c(0,0))) ## vp2
    xy <- list()
    if(track@type %in% c("data", "lollipopData", "interactionData")){
        if(track@type=="data") {
          ##plot yaxis
          drawYaxis(yscale, style@yaxis, curViewStyle)
          pushViewport(viewport(x=curViewStyle@margin[2], y=0, 
                                height=1, 
                                width=1-curViewStyle@margin[2]-curViewStyle@margin[4], 
                                clip="on",
                                just=c(0,0), 
                                xscale=xscale, 
                                yscale=yscale))
          ##grid.clip()
          ##for dat
          xy1 <- plotDataTrack(track@dat, chr, strand, xlim, style@color[1], yscale=yscale, smooth=smooth)
          xy2 <- list(x=numeric(length=0L), y=numeric(length=0L))
          ##for dat2
          if(length(track@dat2)>0){
            if(is.null(operator)) track@dat2$score <- -1 * track@dat2$score ##convert to negtive value
            xy2 <- plotDataTrack(track@dat2, chr, strand, xlim, style@color[2], yscale=yscale, smooth=smooth)
          }
          xy <- list(x=c(xy1$x, xy2$x), y=c(xy1$y, xy2$y))
          xy.id <- order(xy$x)
          xy <- list(x=xy$x[xy.id], y=xy$y[xy.id])
          ##plot xscale?
          if(style@xscale@draw){
            if(length(style@xscale@label)==0){
              scale <- hrScale(xlim)
              xpos <- locateScale(xy$x, abs(xy$y), max(abs(yscale)), scale)
              ypos <- yscale[abs(yscale)==max(abs(yscale))] * 0.7
              if(length(ypos)>1) ypos <- max(ypos)
              style@xscale@from <- new("pos", x=xpos[1], y=ypos, unit="native")
              style@xscale@to <- new("pos", x=xpos[2], y=ypos, unit="native")
              style@xscale@label <- paste(2*scale[[1]], scale[[2]])
            }
            drawXscale(style@xscale)
          }
        }else{
          if(track@type=="lollipopData"){
            pushViewport(viewport(x=curViewStyle@margin[2], y=0, 
                                  height=1, 
                                  width=1-curViewStyle@margin[2]-curViewStyle@margin[4], 
                                  clip="on",
                                  just=c(0,0), 
                                  xscale=xscale))
            ybase <- ifelse(length(track@dat2)>0, .5, 0)
            
            LINEW <- as.numeric(convertX(unit(1, "line"), "npc"))
            LINEH <- as.numeric(convertY(unit(1, "line"), "npc"))
            ## GAP the gaps between any elements
            GAP <- .2 * LINEH
            ratio.yx <- 1/as.numeric(convertX(unit(1, "snpc"), "npc"))
            getMaxHeight <- function(lollipopData){
              if(length(lollipopData)==0) return(0)
              TYPES <- c("circle", "pie", "pin", "pie.stack", "flag")
              type <- if(is.list(lollipopData$type)) lollipopData$type[[1]] else lollipopData$type[1]
              if(length(type)==0) type <- "circle"
              if(!type %in% TYPES) type <- "circle"
              cex <- if(is.list(lollipopData$cex)) lollipopData$cex[[1]] else lollipopData$cex[1]
              if(length(cex)==0) cex <- 1
              scoreMax0 <- scoreMax <- 
                if(length(lollipopData$score)>0) ceiling(max(c(lollipopData$score, 1), na.rm=TRUE)) else 1
              if(type=="pie.stack") scoreMax <- length(unique(lollipopData$stack.factor))
              if(!type %in% c("pie", "pie.stack")){
                if(scoreMax>10) {
                  scoreMax <- 10*scoreMax0/scoreMax
                }else{
                  scoreMax <- scoreMax0
                }
              }
              getHeight(lollipopData, 
                        ratio.yx, LINEW, GAP, cex, type,
                        scoreMax=scoreMax,
                        level="data")
            }
            maxHeight <- max(c(getMaxHeight(track@dat), getMaxHeight(track@dat2)), na.rm = TRUE)
            if(length(track@dat2)>0) maxHeight + .5
            plotLollipopData(track@dat, xlim, chr, style@yaxis@draw, gpar(),
                             ybase, side="top", main=style@yaxis@main,
                             baselineCol=style@color[1], maxHeight=maxHeight)
            if(length(track@dat2)>0) {
              plotLollipopData(track@dat2, xlim, chr, style@yaxis@draw, gpar(),
                               ybase, side="bottom", main=style@yaxis@main,
                               baselineCol=style@color[2], maxHeight=maxHeight)
            }
          }else{##interactionData
            ##plot yaxis
            drawYaxis(yscale, style@yaxis, curViewStyle)
            pushViewport(viewport(x=curViewStyle@margin[2], y=0, 
                                  height=1, 
                                  width=1-curViewStyle@margin[2]-curViewStyle@margin[4], 
                                  clip="on",
                                  just=c(0,0), 
                                  xscale=xscale, 
                                  yscale=yscale))
            ##grid.clip()
            ##for dat interaction: dat, dat2, pair.
            if(length(track@dat)==length(track@dat2)){
              plotInteractionDataTrack(track@dat, track@dat2, xlim, style@color, yscale=yscale, 
                                       breaks=style@breaks, NAcolor=style@NAcolor)
            }
          }
        }
    }else{##track@type=="transcript" or "gene"
      if(track@type %in% c("transcript", "gene")){
        pushViewport(viewport(x=curViewStyle@margin[2], y=0, 
                              height=1, 
                              width=1-curViewStyle@margin[2]-curViewStyle@margin[4], 
                              clip="on",
                              just=c(0,0), 
                              xscale=xscale))
        plotGeneModel(track, xlim, chr)
      }else{##others
        pushViewport(viewport(x=curViewStyle@margin[2], y=0, 
                              height=1, 
                              width=1-curViewStyle@margin[2]-curViewStyle@margin[4], 
                              clip="on",
                              just=c(0,0), 
                              xscale=xscale))
        plotGeneModel(track, xlim, chr) ### currently do the same thing as transcript.
      }
    }
    popViewport()## for data
    
    ## plot wavy line
    if(wavyLine){
        vgap <- convertWidth(unit(0.25, "lines"), 
                             unitTo = "npc",
                             valueOnly = TRUE)
        if(track@type=="data"){
          if(any(yscale==0)){
            y0 <- abs(0-min(yscale))/diff(yscale)
          }else{
            if(all(yscale>0)){
              y0 <- 0
            }else{
              if(all(yscale<0)){
                y0 <- 1
              }else{
                y0 <- abs(0-min(yscale))/diff(yscale)
              }
            }
          }
        }else{
            y0 <- .5
        }
        grid.text(label = "~", x = 1, y=y0, rot = 60)
        grid.text(label = "~", x = 1+vgap, y=y0, rot = 60)
    }
    
    popViewport() ## vp2
    popViewport() ## vp1
    
    ##return
    return(invisible(xy))
}