#' plot GInteractions
#' @description plot graph for GInteractions
#' @param gi an object of \link[InteractionSet:GInteractions-class]{GInteractions}
#' @param range the region to plot. an object of \link[GenomicRanges:GRanges-class]{GRanges}
#' @param ... Not used.
#' @importClassesFrom InteractionSet GInteractions
#' @importMethodsFrom InteractionSet regions subsetByOverlaps anchorIds 
#' @import GenomicRanges
#' @import S4Vectors
#' @importFrom BiocGenerics sort
#' @importClassesFrom graph graphNEL
#' @importFrom Rgraphviz layoutGraph renderGraph
#' @importFrom plotrix arctext
#' @importFrom stats plogis
#' @importFrom graphics curve lines par plot.default points segments strheight text
#' @export
#' @examples 
#' gi <- readRDS(system.file("extdata", "gi.rds", package="trackViewer"))
#' range <- GRanges("chr2", IRanges(1, 100000000))
#' plotGInteractions(gi, range)
plotGInteractions <- function(gi, range, ...){
  stopifnot(is(gi, "GInteractions"))
  if(!missing(range)){
    stopifnot(is(range, "GRanges"))
    gi <- subsetByOverlaps(gi, ranges = range)
  }
  gi <- sort(gi)
  reg <- regions(gi)
  ol <- findOverlaps(reg, drop.self=TRUE, drop.redundant=TRUE, minoverlap = 2)
  if(length(ol)>0){
    warning("There are overlaps in the input regions")
  }
  edgeL <- cbind(anchorIds(gi, type="first"), anchorIds(gi, type="second"))
  nodes <- unique(as.character(sort(as.numeric(edgeL))))
  edgeL <- c(split(edgeL[, 1], as.character(edgeL[, 2])), 
             split(edgeL[, 2], as.character(edgeL[, 1])))
  edgeL <- edgeL[nodes]
  gR <- new("graphNEL", nodes=nodes, edgeL=edgeL)
  gR <- layoutGraph(gR, layoutType = "neato")
  nodeX <- gR@renderInfo@nodes$nodeX
  nodeY <- gR@renderInfo@nodes$nodeY
  stopifnot(identical(names(nodeX), names(nodeY)))
  reg <- reg[as.numeric(names(nodeX))]
  nodesSize <- log(width(reg))
  nodesSize <- nodesSize* mean(gR@renderInfo@nodes$rWidth)/median(nodesSize)
  xlim <- range(nodeX)
  ylim <- range(nodeY)
  xlim <- c(xlim[1]-diff(xlim)/10, xlim[2]+diff(xlim)/10)
  ylim <- c(ylim[1]-diff(ylim)/10, ylim[2]+diff(ylim)/10)
  opar <- par(mar=rep(0, 4)+.1)
  on.exit(par(opar))
  plot.default(0, 0, type="n", asp = 1, xlab="", ylab="", xaxt="n", yaxt="n", 
               xaxs="i", yaxs="i", xlim=xlim, ylim=ylim, frame.plot=FALSE)
  init.angle <- 180*atan(diff(nodeY)/diff(nodeX))/pi + ifelse(diff(nodeX) > 0, 90, -90)
  init.angle <- c(init.angle, 0)
  nodeShape <- list()
  for(i in seq_along(nodeX)){
    nodeShape[[i]] <- archPlot(nodeX[i], nodeY[i], r=nodesSize[i], init.angle = init.angle[i])
  }
  for(i in seq.int(length(nodeShape)-1)){
    curveMaker(nodeShape[[i]]$x[100], nodeShape[[i]]$y[100], 
               nodeShape[[i+1]]$x[1], nodeShape[[i+1]]$y[1], 
               col="gray", lty=2, lwd=2)
  }
  edgeL.m <- data.frame(from=rep(as.numeric(names(edgeL)), lengths(edgeL)), 
                        to=unlist(edgeL))
  edgeL.m <- unique(edgeL.m)
  segments(nodeX[as.character(edgeL.m$from)], nodeY[as.character(edgeL.m$from)], 
           nodeX[as.character(edgeL.m$to)], nodeY[as.character(edgeL.m$to)],
           lwd=3, col = "gray80")
  for(i in seq_along(nodeX)){
    lines(nodeShape[[i]]$x, nodeShape[[i]]$y)
    arctext(x=paste(start(reg)[i], end(reg)[i], sep="-"), center=c(nodeX[i], nodeY[i]), 
            radius = nodesSize[i]+strheight("0-9"), clockwise = FALSE, middle=-pi/2)
    if(i==1){
      points(nodeShape[[i]]$x[1], nodeShape[[i]]$y[1], col="red", pch=16)
    }
    if(i==length(nodeX)){
      points(nodeShape[[i]]$x[100], nodeShape[[i]]$y[100], col="orange", pch=15)
    }
  }
  points(nodeX, nodeY, pch=16, col="white", cex=nodesSize/strheight("0"))
  text(nodeX, nodeY, labels = names(nodeX))
  return(invisible(gR))
}

archPlot <- function(x, y, r, angle=300, init.angle=0){
  gamma <- 2*pi * angle * (1:2) / 360 + init.angle * pi/180
  t2xy <- function(rx, t, x0, y0) list(x=rx*cos(t)+x0, y=rx*sin(t)+y0)
  P <- t2xy(r, seq.int(gamma[1], gamma[2], length.out = 100), x, y)
  #lines(P$x, P$y)
  #points(P$x[1], P$y[1])
  return(P)
}

curveMaker <- function(x1, y1, x2, y2, ...){
  x <- NULL
  if(x1<x2){
    curve( plogis( x, scale = exp(round(log10(mean(c(x1, x2)))-1)), location = (x1 + x2) /2 ) * (y2-y1) + y1, 
           x1, x2, add = TRUE, ...)
  }else{
    curve( plogis( x, scale = exp(round(log10(mean(c(x1, x2)))-1)), location = (x1 + x2) /2 ) * (y1-y2) + y2, 
           x2, x1, add = TRUE, ...)
  }
}
