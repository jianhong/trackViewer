#' plot GInteractions
#' @description plot graph for GInteractions
#' @param gi an object of \link[InteractionSet:GInteractions-class]{GInteractions}
#' @param range the region to plot. an object of \link[GenomicRanges:GRanges-class]{GRanges}
#' @param feature.gr the feature.gr to be added. an object of \link[GenomicRanges:GRanges-class]{GRanges}
#' @param ... Not used.
#' @importClassesFrom InteractionSet GInteractions
#' @importMethodsFrom InteractionSet regions anchorIds
#' @import GenomicRanges
#' @import S4Vectors
#' @importFrom BiocGenerics sort
#' @importClassesFrom graph graphNEL
#' @importFrom Rgraphviz layoutGraph renderGraph
#' @importFrom plotrix arctext
#' @importFrom stats plogis
#' @importFrom graphics curve lines par points segments strheight text arrows
#' @export
#' @examples
#' library(InteractionSet) 
#' gi <- readRDS(system.file("extdata", "gi.rds", package="trackViewer"))
#' range <- GRanges("chr2", IRanges(234500000, 235000000))
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' feature.gr <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' feature.gr <- subsetByOverlaps(feature.gr, regions(gi))
#' feature.gr$col <- sample(1:7, length(feature.gr), replace=TRUE)
#' feature.gr$type <- sample(c("promoter", "enhancer", "gene"), 
#'                          length(feature.gr), replace=TRUE, 
#'                          prob=c(0.1, 0.2, 0.7))
#' plotGInteractions(gi, range, feature.gr)
plotGInteractions <- function(gi, range, feature.gr, ...){
  stopifnot(is(gi, "GInteractions"))
  if(!missing(range)){
    stopifnot(is(range, "GRanges"))
    gi <- subsetByOverlaps(gi, ranges = range)
  }
  if(!missing(feature.gr)){
    stopifnot(is(feature.gr, "GRanges"))
  }else{
    feature.gr <- GRanges()
  }
  gi <- sort(gi)
  reg <- regions(gi)
  names(reg) <- seq_along(reg)
  feature.gr <- subsetByOverlaps(feature.gr, reg)
  if(length(feature.gr$col)==0) feature.gr$col <- rep("black", length(feature.gr))
  if(length(feature.gr$type)==0) feature.gr$type <- rep("gene", length(feature.gr))
  feature.gr <- promoters(feature.gr, upstream = 0, downstream = 1)
  ol <- findOverlaps(reg, drop.self=TRUE, drop.redundant=TRUE, minoverlap = 2)
  if(length(ol)>0){
    warning("There are overlaps in the input regions")
  }
  ol <- findOverlaps(feature.gr, reg)
  feature.gr <- split(feature.gr[queryHits(ol)], subjectHits(ol))
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
  reg <- reg[names(nodeX)]
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
  nodeShape2 <- list()
  spaces <- list()
  sH <- strheight("0")
  arrowLen <- convertUnit(stringHeight("0"), unitTo = "inches", valueOnly = TRUE)/3
  for(i in seq_along(nodeX)){
    nodeShape[[i]] <- archPlot(nodeX[i], nodeY[i], r=nodesSize[i], init.angle = init.angle[i])
    nodeShape2[[i]] <- archPlot(nodeX[i], nodeY[i], r=nodesSize[i]-sH, init.angle = init.angle[i])
    nodeShape3<- archPlot(nodeX[i], nodeY[i], r=nodesSize[i]-sH, init.angle = init.angle[i], angle = 360, length.out=120)
    dist <- (nodeShape3$x - nodeShape2[[i]]$x[1])^2 + (nodeShape3$y - nodeShape2[[i]]$y[1])^2
    dist <- which.min(dist)[1]
    nodeShape2[[i]] <- list(x=c(nodeShape3$x[dist:120], nodeShape3$x[seq.int(dist)]), 
                            y=c(nodeShape3$y[dist:120], nodeShape3$y[seq.int(dist)]))
    spaces[[i]] <- sH*120/(pi*2*(nodesSize[i]-sH))
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
    if(length(feature.gr[[names(nodeX)[i]]])>0){
      addGeneInNode(feature.gr[[names(nodeX)[i]]], nodeShape[[i]], nodeShape2[[i]], nodeX[i], nodeY[i], reg[names(nodeX)[i]], spaces[[i]], arrowLen)
    }
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

archPlot <- function(x, y, r, angle=300, init.angle=0, length.out=100){
  gamma <- 2*pi * angle * (1:2) / 360 + init.angle * pi/180
  t2xy <- function(rx, t, x0, y0) list(x=rx*cos(t)+x0, y=rx*sin(t)+y0)
  P <- t2xy(r, seq.int(gamma[1], gamma[2], length.out = length.out), x, y)
  #lines(P$x, P$y)
  #points(P$x[1], P$y[1])
  return(P)
}

curveMaker <- function(x1, y1, x2, y2, ...){
  x <- NULL
  if(abs(x1-x2)<10){
    segments(x1, y1, x2, y2, ...)
  }else{
    if(x1<x2){
      curve( plogis( x, scale = exp(round(log10(mean(c(x1, x2)))-1)), location = (x1 + x2) /2 ) * (y2-y1) + y1, 
             x1, x2, add = TRUE, ...)
    }else{
      curve( plogis( x, scale = exp(round(log10(mean(c(x1, x2)))-1)), location = (x1 + x2) /2 ) * (y1-y2) + y2, 
             x2, x1, add = TRUE, ...)
    }
  }
}

addGeneInNode <- function(gs, nShape, nShape2, x, y, gr, space, arrowLen){
  gs <- split(gs, gs$type)
  n <- length(nShape$x)
  gr.cut <- seq(start(gr), end(gr), length.out = n)
  mapply(function(gs.s, type){
    if(length(gs.s)>0){
      at <- as.numeric(cut(start(gs.s), breaks = gr.cut, labels = seq.int(n-1)))
      A.x <- nShape$x[at]
      A.y <- nShape$y[at]
      B.x <- nShape2$x[at]
      B.y <- nShape2$y[at]
      switch(type,
             "enhancer"={
               points(A.x, A.y, pch=11, col = gs.s$col)
             }, "promoter"={
               points(A.x, A.y, pch=13, col = gs.s$col)
             }, "gene"={
               segments(A.x, A.y, B.x, B.y, col = gs.s$col)
               at1 <- at + ifelse(as.character(strand(gs.s))=="+", space, -1*space)
               if(any(at1<1)) at1[at1<1] <- 120 + at1[at1<1]
               segments(B.x, B.y, nShape2$x[at1], nShape2$y[at1], col = gs.s$col)
               arrows(B.x, B.y, nShape2$x[at1], nShape2$y[at1], length = arrowLen, col = gs.s$col)
             },{
               segments(A.x, A.y, B.x, B.y, col = gs.s$col)
               at1 <- at + ifelse(as.character(strand(gs.s))=="+", space, -1*space)
               segments(B.x, B.y, nShape2$x[at1], nShape2$y[at1], col = gs.s$col)
               arrows(B.x, B.y, nShape2$x[at1], nShape2$y[at1], length = arrowLen, col = gs.s$col)
             })
    }
  }, gs, names(gs))
}