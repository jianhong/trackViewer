plotDataTrack <- function(.dat, chr, strand, scale, color){
    names(.dat) <- NULL
    mcols(.dat) <- mcols(.dat)[, "score"]
    colnames(mcols(.dat)) <- "score"
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
            grid.polygon(x2, y2, default.units="native", 
                         gp=gpar(col=NA, fill=color))
            grid.lines(x=scale, y=0, default.units="native",
                       gp=gpar(col=color))
            xt <- c(xt, x)
            yt <- c(yt, y)
        }
    }
    return(list(x=xt, y=yt))
}
plotTrack <- function(name, track, curViewStyle, curYpos,
                      yscale, height, xlim, chr, strand,
                      operator, wavyLine){
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
                      range(track@dat))
        } 
    }else{
        if(style@ylabpos=="upstream") style@ylabpos <- "left"
        if(style@ylabpos=="downstream") style@ylabpos <- "right"
        putYlab(curViewStyle, style, name, yHeightBottom, yHeightTop, height)
    }
    
    pushViewport(viewport(x=0, y=yHeightBottom, 
                          height=1 - yHeightTop - yHeightBottom,
                          width=1, 
                          just=c(0,0))) ## vp2
    xy <- list()
    if(track@type %in% c("data", "lollipopData")){
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
          xy1 <- plotDataTrack(track@dat, chr, strand, xlim, style@color[1])
          xy2 <- list(x=numeric(length=0L), y=numeric(length=0L))
          ##for dat2
          if(length(track@dat2)>0){
            if(is.null(operator)) track@dat2$score <- -1 * track@dat2$score ##convert to negtive value
            xy2 <- plotDataTrack(track@dat2, chr, strand, xlim, style@color[2])
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
          pushViewport(viewport(x=curViewStyle@margin[2], y=0, 
                                height=1, 
                                width=1-curViewStyle@margin[2]-curViewStyle@margin[4], 
                                clip="on",
                                just=c(0,0), 
                                xscale=xscale))
          ybase <- ifelse(length(track@dat2)>0, .5, 0)
          if(length(track@dat)>0) plotLollipopData(track@dat, xlim, chr, style@yaxis@draw, 
                                                   ybase, side="top", main=style@yaxis@main,
                                                   style@color[1])
          if(length(track@dat2)>0) {
            plotLollipopData(track@dat2, xlim, chr, style@yaxis@draw, 
                             ybase, side="bottom", main=style@yaxis@main,
                             baselineCol=style@color[2])
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
            y0 <- abs(0-min(yscale))/diff(yscale)
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