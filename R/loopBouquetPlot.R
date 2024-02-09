#' plot GInteractions
#' @description plot graph for GInteractions
#' @param gi An object of \link[InteractionSet:GInteractions-class]{GInteractions}
#' @param range The region to plot. an object of \link[GenomicRanges:GRanges-class]{GRanges}
#' @param feature.gr The annotation features to be added. An object of \link[GenomicRanges:GRanges-class]{GRanges}.
#' @param atacSig The ATAC-seq signals. An object of \link[GenomicRanges:GRanges-class]{GRanges} with scores or an object of \link{track}.
#' @param label_region Label the region node or not.
#' @param show_edges Plot the interaction edges or not.
#' @param show_cluster Plot the cluster background or not.
#' @param show_coor Plot ticks in the line to show the DNA compact tension.
#' @param reverseATACSig Plot the ATAC-seq signals in reverse values.
#' @param coor_tick_unit The bps for every ticks. Default is 1K.
#' @param coor_mark_interval The coordinates marker interval. Numeric(1). Set to 0
#' to turn it off. The default value 1e5 means show coordinates every 0.1M bp.
#' @param label_gene Show gene symbol or not.
#' @param lwd.backbone,lwd.gene,lwd.nodeCircle,lwd.edge,lwd.tension_line,lwd.maxAtacSig Line width for the 
#' linker, gene, interaction node circle, the dashed line of interaction edges, the tension line and the maximal reversed ATAC signal.
#' @param col.backbone,col.backbone_background,col.nodeCircle,col.edge,col.tension_line,col.coor Color
#' for the DNA chain, the compact DNA chain, the node circle, the linker, the tension line and the coordinates marker.
#' @param length.arrow Length of the edges of the arrow head (in inches).
#' @param safe_text_force The loops to avoid the text overlapping.
#' @param method Plot method. Could be 1 or 2.
#' @param ... Parameter will be passed to \link[igraph:layout_with_fr]{layout_with_fr}.
#' @importClassesFrom InteractionSet GInteractions
#' @importMethodsFrom InteractionSet regions anchorIds
#' @import GenomicRanges
#' @import S4Vectors
#' @importFrom BiocGenerics sort
#' @importFrom igraph graph_from_data_frame components layout_with_fr
#' norm_coords V
#' @importFrom plotrix arctext
#' @importFrom stats plogis
#' @importFrom graphics curve lines par points segments strheight text arrows
#' polygon
#' @importFrom scales rescale
#' @export
#' @examples
#' library(InteractionSet) 
#' gi <- readRDS(system.file("extdata", "gi.rds", package="trackViewer"))
#' range <- GRanges("chr2", IRanges(234500000, 235000000))
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' library(org.Hs.eg.db)
#' feature.gr <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' feature.gr <- subsetByOverlaps(feature.gr, range(regions(gi)))
#' symbols <- mget(feature.gr$gene_id, org.Hs.egSYMBOL, ifnotfound=NA)
#' feature.gr$label[lengths(symbols)==1] <- unlist(symbols[lengths(symbols)==1])
#' feature.gr$col <- sample(1:7, length(feature.gr), replace=TRUE)
#' feature.gr$type <- sample(c("cRE", "gene"), 
#'                          length(feature.gr), replace=TRUE, 
#'                          prob=c(0.1, 0.9))
#' loopBouquetPlot(gi, range, feature.gr)
loopBouquetPlot <- function(gi, range, feature.gr, atacSig,
                            label_region=FALSE, show_edges=TRUE, 
                            show_cluster=TRUE,
                            lwd.backbone = 2, col.backbone = 'gray',
                            lwd.maxAtacSig = 8, reverseATACSig = TRUE,
                            col.backbone_background = 'gray70',
                            lwd.gene = 2,
                            lwd.nodeCircle = 1, col.nodeCircle = '#DDDDDD25',
                            lwd.edge = 2, col.edge = "gray80",
                            coor_mark_interval = 1e5, col.coor = "black",
                            show_coor = TRUE,
                            coor_tick_unit = 1e3,
                            label_gene = TRUE,
                            col.tension_line = 'black',
                            lwd.tension_line = 1,
                            length.arrow = NULL,
                            safe_text_force = 3,
                            method = 1,
                            ...){
  stopifnot(is(gi, "GInteractions"))
  stopifnot('score' %in% colnames(mcols(gi)))
  stopifnot(is.numeric(coor_mark_interval))
  stopifnot(length(coor_mark_interval)==1)
  if(!missing(range)){
    stopifnot(is(range, "GRanges"))
    stopifnot('coor_tick_unit is too small.'=
                width(range(range))[1]/coor_tick_unit<1000)
    gi <- subsetByOverlaps(gi, ranges = range, ignore.strand=TRUE)
  }
  if(!missing(feature.gr)){
    stopifnot(is(feature.gr, "GRanges"))
  }else{
    feature.gr <- GRanges()
  }
  gi <- sort(gi)
  reg <- regions(gi)
  stopifnot('No interactions detected'=length(reg)>2)
  names(reg) <- seq_along(reg)
  if(!all(as.character(seqnames(reg))==as.character(seqnames(reg))[1])){
    warning('All interaction must within one chromosome.
            Interchromosomal interactions will be dropped.')
  }
  ol <- findOverlaps(reg, drop.self=TRUE, drop.redundant=TRUE, minoverlap = 2)
  if(length(ol)>0){
    warning("There are overlaps in the input regions")
  }
  if(length(feature.gr$col)==0) feature.gr$col <- rep("black", length(feature.gr))
  if(length(feature.gr$type)==0) feature.gr$type <- rep("gene", length(feature.gr))
  if(length(feature.gr$cex)==0) feature.gr$cex <- rep(1, length(feature.gr))
  stopifnot(length(feature.gr$type)==length(feature.gr))
  stopifnot(all(feature.gr$type %in% c('cRE', 'gene')))
  feature.gr[feature.gr$type=='cRE'] <- 
    promoters(feature.gr[feature.gr$type=='cRE'], upstream = 0, downstream = 1)
  
  nodes <- unique(as.character(sort(c(anchorIds(gi, type="first"),
                                      anchorIds(gi, type="second")))))
  ## add genomic coordinates to edge
  d0 <- distance(reg[nodes[-length(nodes)]],
                 reg[nodes[-1]])
  d0.ups.dws <- width(reg[nodes[-length(nodes)]]) + width(reg[nodes[-1]])
  d <- log10(d0 + 1) + 1
  # d_factor <- median(width(reg[nodes]))
  # d <- distance(reg[nodes[-length(nodes)]], reg[nodes[-1]])/d_factor + 1
  irq <- quantile(d, probs = c(.25, .75))
  edgeL_coor <- data.frame(
    from = nodes[-length(nodes)],
    to = nodes[-1],
    weight = max(d)/d
  )
  edgeL_link <- data.frame(
    from = as.character(anchorIds(gi, type="first")),
    to = as.character(anchorIds(gi, type="second")),
    weight = gi$score + 2*diff(irq) + irq[2])
  
  edgeL <- rbind(edgeL_coor, edgeL_link)
  m_w_reg <- min(width(reg[nodes]))
  nodes <- data.frame(names = nodes,
                      size = (width(reg[nodes]))/ifelse(m_w_reg==0, 1, m_w_reg))
  gL <- graph_from_data_frame(d = edgeL_link, directed=FALSE, vertices = nodes)
  cl <- igraph::components(gL, mode = "weak")
  sGnodes <- split(names(cl$membership), cl$membership)
  g <- graph_from_data_frame(d = edgeL, directed=FALSE, vertices = nodes)
  layout <- layout_with_fr(g, ...) ## only layout_with_fr and layout_with_kk work OK
  stopifnot('3 dim is not supported yet'=ncol(layout)==2)
  nodeXY <- norm_coords(layout, -1, 1, -1, 1)
  rownames(nodeXY) <- nodes$names
  colnames(nodeXY) <- c("X", "Y")
  vertex.factor <- 72
  vertex.size <- 1/vertex.factor * V(g)$size
  vertex.size[is.na(vertex.size)] <- 1/vertex.factor
  nodeXY <- fixXY(nodeXY, vertex.size, edgeL_link, lwd = lwd.backbone/300)
  maxv <- max(vertex.size)
  xlim <- range(nodeXY[, 1])
  ylim <- range(nodeXY[, 2])
  d_xlim <- diff(xlim)/5
  d_ylim <- diff(ylim)/5
  xlim <- c(xlim[1] - d_xlim - maxv, xlim[2] + d_xlim + maxv)
  ylim <- c(ylim[1] - d_ylim - maxv, ylim[2] + d_ylim + maxv)

  nodesSize <- nodes$size
  
  clusterCenter <- lapply(sGnodes, function(.ele){
    X <- nodeXY[.ele, "X"]
    Y <- nodeXY[.ele, "Y"]
    list(x=mean(X), y=mean(Y),
         r=sqrt((diff(range(X)/2) + median(vertex.size))^2 +
                  (diff(range(Y)/2) + median(vertex.size))^2),
         nodes=.ele)
  })
  
  opar <- par(mar=rep(0, 4)+.1)
  on.exit(par(opar))
  
  plotPoints <- list()

  # plot(0, 0, type="n", asp = 1, xlab="", ylab="", xaxt="n", yaxt="n",
  #      xaxs="i", yaxs="i", xlim=xlim, ylim=ylim, frame.plot=FALSE,
  #      new=TRUE)
  grid.newpage()
  vp <- viewport(default.units = 'native',
                 width=unit(min(1,diff(xlim)/diff(ylim)), "snpc"), # aspect ratio preserved
                 height=unit(min(1,diff(ylim)/diff(xlim)), "snpc"),
                 xscale=xlim, yscale=ylim)
  pushViewport(vp)
  on.exit(popViewport())
  if(!is.numeric(length.arrow)){
    sH <- stringWidth("0")
    arrowLen <- grid::convertUnit(grid::stringHeight("0"), unitTo = "inches", valueOnly = TRUE)/3
  }else{
    arrowLen <- length.arrow[1]
  }
  
  ## make sure the init.angle is facing to next node
  # init.angle <- 180*atan(diff(nodeY)/diff(nodeX))/pi +
  #   ifelse(diff(nodeX) >= 0, 90, -90)
  # init.angle <- c(init.angle, 0)
  ## make the init.angle is facing outside of the clusterCenter
  init.angle <- lapply(clusterCenter, function(.ele){
    dX <- nodeXY[.ele$nodes, "X"] - .ele$x
    dY <- nodeXY[.ele$nodes, "Y"] - .ele$y
    180*atan(dY/dX)/pi + ifelse(dX>=0, 90, -90)
  })
  init.angle <- c(unlist(unname(init.angle))[nodes$names], 0)
  nodeShape <- list() ## the curves connected the nodes
  for(i in seq.int(nrow(nodeXY))){
    nodeShape[[i]] <- archPlot(nodeXY[i, "X"], nodeXY[i, "Y"], r=nodesSize[i]/vertex.factor, init.angle = init.angle[i])
  }
  
  getNodesClusters <- function(...){
    which(vapply(sGnodes, function(.ele) any(... %in% .ele), logical(1L)))
  }
  last_a <- 1
  for(i in seq.int(length(nodeShape))){
    ## plot the curve connected the nodes
    if(i < length(nodeShape)){
      ## check the connection points
      if(i==1){ ## mark the start point
        reg0 <- reg[rownames(nodeXY)[i]]
        this_reg_gap <- GRanges(seqnames(reg0), IRanges(end=start(reg0),
                                                        width = 10*coor_tick_unit))
        idx <- c(i, i+1)
        ab <- getStartConnectionPoints(nodeShape[idx])
        idx <- ifelse(ab[1]==1, length(nodeShape[[1]]$x), 1)
        thisXY <- plotStartEndPoints(nodeShape[[i]], idx, 
                           xlim = xlim, ylim = ylim)
        plotPoints <- addPoints(plotPoints,
                               thisXY$x, thisXY$y, 
                               start(this_reg_gap),
                               end(this_reg_gap))
      }else{
        idx <- seq(i, min(i+2, length(nodeShape)))
        ab <- checkConnectionPoints(nodeShape[c(i, i+1)], last_a, d[i]>median(d))
      }
      ## The curve will be plot by the center of the subcluster
      ## or the center of the two nodes
      cn <- getNodesClusters(nodes$names[c(i, i+1)])
      if(length(cn)==1){
        cn <- clusterCenter[[cn]]
      }else{
        cn <- lapply(clusterCenter[cn], function(.ele)
          c(x=.ele$x, y=.ele$y, r=.ele$r))
        cn <- as.data.frame(do.call(rbind, cn))
      }
      reg0 <- reg[rownames(nodeXY)[i]]
      reg1 <- reg[rownames(nodeXY)[i+1]]
      if(ab[1]==1){
        plotPoints <- addPoints(plotPoints,
                                rev(nodeShape[[i]]$x),
                                rev(nodeShape[[i]]$y),
                                start(reg0),
                                end(reg0))
      }else{
        plotPoints <- addPoints(plotPoints,
                                nodeShape[[i]]$x,
                                nodeShape[[i]]$y,
                                start(reg0),
                                end(reg0))
      }
      thisXY <- curveMaker(nodeShape[[i]], 
                 nodeShape[[i+1]],
                 ab,
                 cn$x, cn$y, r=cn$r,
                 w=d0[i]/d0.ups.dws[i],
                 evaluation = floor(100*edgeL_coor$weight[i]),
                 method = method)
      plotPoints <- addPoints(plotPoints,
                              thisXY$x, thisXY$y, 
                              end(reg0),
                              start(reg1))
      last_a <- ifelse(ab[2]==1, length(nodeShape[[i]]$x), 1)
    }
    if(i==length(nodeShape)){ ## mark the end point
      reg0 <- reg[rownames(nodeXY)[i]]
      this_reg_gap <- GRanges(seqnames(reg0), IRanges(end(reg0),
                                                      width = 10*coor_tick_unit))
      if(last_a==1){
        plotPoints <- addPoints(plotPoints,
                                rev(nodeShape[[i]]$x),
                                rev(nodeShape[[i]]$y),
                                start(reg0),
                                end(reg0))
      }else{
        plotPoints <- addPoints(plotPoints,
                                nodeShape[[i]]$x,
                                nodeShape[[i]]$y,
                                start(reg0),
                                end(reg0))
      }
      thisXY <- plotStartEndPoints(nodeShape[[i]], last_a, 
                         xlim = xlim, ylim = ylim, 
                         start = nodeShape[[i]]$y[1])
      plotPoints <- addPoints(plotPoints,
                              thisXY$x, thisXY$y, 
                              start(this_reg_gap),
                              end(this_reg_gap))
    }
  }
  ## plot the bacground circle
  if(show_cluster){
    for(i in seq_along(clusterCenter)){
      grid.circle(clusterCenter[[i]]$x, clusterCenter[[i]]$y,
                  clusterCenter[[i]]$r,
                   default.units = "native",
              gp=gpar(fill = col.nodeCircle,
                      lty = 4, lwd=lwd.nodeCircle))
    }
  }
  ## plot the edge lines
  if(show_edges){
    grid.segments(nodeXY[as.character(edgeL_link$from), "X"],
                  nodeXY[as.character(edgeL_link$from), "Y"], 
                  nodeXY[as.character(edgeL_link$to), "X"],
                  nodeXY[as.character(edgeL_link$to), "Y"],
                  default.units = "native",
                  gp=gpar(lwd=lwd.edge, lty=2, col = col.edge))
  }
  ## plot the bouquet
  plotBouquet(pP=plotPoints,
              fgf=feature.gr,
              atacSig=atacSig,
              lwd.backbone=lwd.backbone,
              col.backbone=col.backbone,
              lwd.maxAtacSig = lwd.maxAtacSig,
              reverseATACSig = reverseATACSig,
              col.backbone_background=col.backbone_background,
              lwd.gene=lwd.gene,
              coor_mark_interval=coor_mark_interval,
              col.coor=col.coor,
              show_coor = show_coor,
              coor_tick_unit = coor_tick_unit,
              label_gene = label_gene,
              col.tension_line = col.tension_line,
              lwd.tension_line = lwd.tension_line,
              safe_text_force = safe_text_force,
              arrowLen = arrowLen,
              xlim=xlim, ylim=ylim)
  ## 
  if(label_region){
    ## plot the node dot
    grid.points(nodeXY[, 'X'],
                nodeXY[, 'Y'],
                pch=16,
                default.units = "native",
                gp=gpar(col="white",
                        cex=nodesSize/strheight("0")))
    ## label the node
    grid.text(nodeXY[, 'X'], nodeXY[, 'Y'],
              default.units = "native",
              label = rownames(nodeXY))
  }
  
  return(invisible(list(plotPoints=plotPoints, feature.gr=feature.gr,
                        xlim=xlim, ylim=ylim,
                        nodeXY=nodeXY, edgeL_link=edgeL_link,
                        clusterCenter=clusterCenter)))
}

addPoints <- function(pointscollection, x, y, s, e){
  pointscollection[[length(pointscollection)+1]] <- list(
    x=x, y=y, start=s, end=e)
  pointscollection
}

getStartConnectionPoints <- function(shapeX){
  ## return the nearest pairs
  ab <- c(1, length(shapeX[[1]]$x))
  p <- lapply(shapeX, function(.ele) data.frame(x=.ele$x[ab], y=.ele$y[ab]))
  pD <- apply(p[[1]], 1, function(p1) {
    apply(p[[2]], 1, function(p2) pointsDistance(p1, p2))
  }, simplify = FALSE)
  pD <- unlist(pD)
  cn <- list(c(1, 1), c(1, 2), c(2, 1), c(2, 2))
  cn <- cn[[which(pD==min(pD))[1]]]
  return(ab[cn]) ## return the connect pair (from to)
}
checkConnectionPoints <- function(shapeX, a, close = TRUE){
  ## check distance, which one will follow the distance of a, b points
  ## if the bps is long, select the far point,
  ## else select the close point
  ab0 <- c(1, length(shapeX[[1]]$x))
  ab <- c(a, 1)
  p1 <- c(shapeX[[1]]$x[1], shapeX[[1]]$y[1])
  p2 <- data.frame(x=shapeX[[2]]$x[ab0], y=shapeX[[2]]$y[ab0])
  pD <- apply(p2, 1, pointsDistance, p1=p1)
  i <- which(pD==min(pD))[1]
  if(close){
    ab[2] <- ab0[i]
  }else{
    ab[2] <- ab0[-i]
  }
  return(ab)
}

safeEndPoint <- function(x, lim, dx){
  if(x<lim[1]) x <- lim[1] + dx
  if(x>lim[2]) x <- lim[2] - dx
  x
}

plotStartEndPoints <- function(xy, idx, xlim, ylim, start, ...){
  # if it is the end point, the start will be the y of start point
  ## plot a horizontal line to the edge
  x <- xy$x[idx]
  y <- xy$y[idx]
  
  if(idx==1){
    idx <- c(1, 2)
  }else{
    idx <- c(idx-1, idx)
  }
  xs <- xy$x[idx]
  ys <- xy$y[idx]
  slope <- diff(ys)/diff(xs)
  b <- y - slope*x
  
  x0 <- ifelse(x>mean(xlim), xlim[2], xlim[1])
  y0 <- ifelse(y>mean(ylim), ylim[2], ylim[1])
  dx <- pointsDistance(c(x, y), c(x0, y))
  dy <- pointsDistance(c(x, y), c(x, y0))
  ddx <- diff(xlim)/100
  ddy <- diff(ylim)/100
  if(dx<=dy){
    # cross y axis
    x1 <- ifelse(x>mean(xlim), xlim[2]-ddx, xlim[1]+ddx)
    if(slope==0){
      y1 <- y0 <- y
      if(!missing(start)){
        if(start[1]==y){
          y1 <- y0 <- y + ddy
        }
      }
    }else{
      y0 <- slope * x0 + b
      y1 <- slope * x1 + b
    }
    y0 <- safeEndPoint(y0, ylim, 0)
    y1 <- safeEndPoint(y1, ylim, ddy)
  }else{
    # cross x axis
    y1 <- ifelse(y>mean(ylim), ylim[2]-ddy, ylim[1]+ddy)
    x0 <- (y0 - b)/slope
    x1 <- (y1 - b)/slope
    x0 <- safeEndPoint(x0, xlim, 0)
    x1 <- safeEndPoint(x1, xlim, ddx)
  }

  if(!missing(start)){
    return(list(x=c(x, x0), y=c(y, y0)))
  }else{
    return(list(x=c(x0, x), y=c(y0, y)))
  }
}

pointsDistance <- function(p1, p2){
  stopifnot(is.numeric(p1))
  stopifnot(is.numeric(p2))
  sqrt(sum((p1[c(1, 2)] - p2[c(1, 2)])^2))
}
overlapNodes <- function(from, to, lwd){
  pointsDistance(from, to) <= from[3] + to[3] + lwd
}
moveNodes <- function(from, to, lwd, reverse = FALSE){
  # from, to, data structure: c(x, y, r)
  lwd <- lwd/2
  center <- colMeans(rbind(from, to))[c(1, 2)]
  pd <- pointsDistance(from[c(1, 2)], center)
  if(pd!=0){
    from <- from[c(1, 2)] - (from[c(1, 2)] - center) * 
      (1 - (from[3] + lwd) / pd)
  }else{
    from <- from[c(1, 2)]
  }
  pd <- pointsDistance(to[c(1, 2)], center)
  if(pd!=0){
    to <- to[c(1, 2)] - (to[c(1, 2)] - center) * 
      (1 - (to[3] + lwd) / pd)
  }else{
    to <- to[c(1, 2)]
  }
  return(list(from = from, to = to))
}
reorderEdgeLinkBySize <- function(edgeL_link, nodeXYout){
  if(nrow(edgeL_link)==0) return(edgeL_link)
  size <- apply(edgeL_link, 1, function(.ele){
    sum(nodeXYout[.ele[c(1, 2)], 3])
  })
  edgeL_link[order(size, decreasing = TRUE), ]
}
fixXY <- function(nodeXY, vertex.size, edgeL_link, lwd = 10/300){
  ## rearrange nodes to make sure they are not overlap
  ## and make the interaction of nodes direct touch each other
  if(nrow(edgeL_link)==0){
    return(nodeXY)
  }
  nodeXYout <- cbind(nodeXY, vertex.size)
  edgeL_link <- reorderEdgeLinkBySize(edgeL_link, nodeXYout)
  ## touch each other
  for(i in seq.int(nrow(edgeL_link))){
    from <- nodeXYout[edgeL_link[i, 'from'], ]
    to <- nodeXYout[edgeL_link[i, 'to'], ]
    if(!overlapNodes(from, to, lwd=lwd)){
      newCoor <- moveNodes(from, to, lwd=lwd)
      nodeXYout[edgeL_link[i, 'from'], c(1, 2)] <- newCoor[["from"]]
      nodeXYout[edgeL_link[i, 'to'], c(1, 2)] <- newCoor[["to"]]
    }
  }
  ## Do not touch each other
  for(i in seq.int(nrow(edgeL_link))){
    from <- nodeXYout[edgeL_link[i, 'from'], ]
    to <- nodeXYout[edgeL_link[i, 'to'], ]
    if(overlapNodes(from, to, lwd=0)){
      newCoor <- moveNodes(from, to, lwd=lwd, reverse=TRUE)
      nodeXYout[edgeL_link[i, 'from'], c(1, 2)] <- newCoor[["from"]]
      nodeXYout[edgeL_link[i, 'to'], c(1, 2)] <- newCoor[["to"]]
    }
  }
  return(nodeXYout[, c(1, 2)])
}

# 2D Points P=[x,y] and R=[x,y] are arbitrary points on line, 
# Q=[x,y] is point for which we want to find reflection
# returns solution vector [x,y] 
# Q' = Q + 2(I - d^2*v*t(v))(P - Q)
# Where column vector v=R-P is parallel to line but is not unit (like n),
# P and R are two arbitrary line points; the v*t(v) is outer product;
# I is identity matrix; the d=1/sqrt(vx*vx + vy*vy) is inverse length of v.
# N is the tension strength.
mirrorP <- function(Q, P, R, N) {
  I <- matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2)
  n <- R - P
  d2 <- 1/sum(n^2)
  as.numeric(Q + 2*N * ( I - d2 * n %*% t(n) ) %*% ( P - Q ))
}

bezier0 <- function(t, ctrNodes){
  x <- 0
  y <- 0
  n <- nrow(ctrNodes) - 1
  for(i in seq.int(nrow(ctrNodes))){
    index <- i - 1
    if(!index){
      x <- x + ctrNodes[i, 1] * (( 1 - t ) ^ (n - index)) * (t^index) 
      y <- y + ctrNodes[i, 2] * (( 1 - t ) ^ (n - index)) * (t^index) 
    }else{
      x <- x + factorial(n) / factorial(index) / factorial(n - index) * ctrNodes[i, 1] * 
        (( 1 - t ) ^ (n - index)) * (t^index) 
      y <- y + factorial(n) / factorial(index) / factorial(n - index) * ctrNodes[i, 2] * 
        (( 1 - t ) ^ (n - index)) * (t^index) 
    }
  }
  return(c(x, y))
}

bezier <- function(coefs, center, r_coefs, r, evaluation=100, w, method=1){
  stopifnot('Dimention of coefs need to be 2'=length(dim(coefs))==2)
  stopifnot('Dimention of center need to be 2'=length(dim(center))==2)
  # mirror center
  center <- colMeans(center)
  # N is the strength of the bezier curve
  # strength is determined by the distance from the center of (p1, p2)
  # to the sum of perimeter of both nodes
  if(method==2){
    if(r==0){
      N <- min(w, 10)
    }else{
      # w is the ratio of distance to both node width.
      # distance of both nodes center points
      N <- pointsDistance(coefs[1, ], coefs[nrow(coefs), ])
      # distance of the center of nodes center points to sub-cluster center
      Nr <- sqrt((mean(coefs[, 1]) - center[1])^2 +
                   (mean(coefs[, 2]) - center[2])^2)
      # N/sum(2*pi*r_coefs) should be same as w
      N <- max(min(2*w, r/Nr, 2*N/sum(2*pi*r_coefs), na.rm = TRUE), 1.25)
    }
    N <- 2*N # 4 points in new algorithm need 2 times strengh.
  }else{
    if(r==0){
      N <- 2
    }else{
      Nr <- sqrt((mean(coefs[, 1]) - center[1])^2 +
                   (mean(coefs[, 2]) - center[2])^2)
      N <- r/Nr
    }
  }
  
  center_mirror <- mirrorP(center, coefs[1, ], coefs[nrow(coefs), ], N=N)
  half <- floor(nrow(coefs)/2)
  ctrNodes <- rbind(coefs[seq.int(half), ],
                    center_mirror,
                    coefs[seq.int(nrow(coefs))[-seq.int(half)], ])
  xy <- lapply((seq.int(101)-1)/100, bezier0, ctrNodes=ctrNodes)
  xy <- do.call(rbind, xy)
  list(x = xy[, 1], y = xy[, 2])
}

getMirrorPoint <- function(ns, a){
  n <- length(ns$x)
  b <- ifelse(a==1, 10, n - 10)
  c <- rotatePoint(ns$x[a], ns$y[a], ns$x[b], ns$y[b])
  mirrorP(c(ns$x[b], ns$y[b]), c(ns$x[a], ns$y[a]), c, N=1)
}

getPointInNodeCircle <- function(ns, a, c1, r){
  t2xy <- function(rx, t, x0, y0) list(x=rx*cos(t)+x0, y=rx*sin(t)+y0)
  t <- atan((ns$y[a]-c1[2])/(ns$x[a]-c1[1]))
  ap <- t2xy(r, t+ifelse(ns$x[a]-c1[1]>=0, 0, pi), c1[1], c1[2])
  c(ap$x[1], ap$y[1])
}

curveMaker <- function(ns1, ns2, ab, cx, cy, r, w, evaluation = 100,
                       method = 1,
                       ...){
  x1 <- ns1$x[ab[1]]
  y1 <- ns1$y[ab[1]]
  x2 <- ns2$x[ab[2]]
  y2 <- ns2$y[ab[2]]
  c1 <- c(cx[1], cy[1])
  if(length(r)==1){
    c2 <- c1
  }else{
    c2 <- c(cx[2], cy[2])
  }
  if(pointsDistance(c(x1, y1), c1) < r[1]){
    # find the extension point on the circle
    ep1 <- getPointInNodeCircle(ns1, ab[1], c1, r[1])
  }else{
    ep1 <- getMirrorPoint(ns1, ab[1])
  }
  r2 <- ifelse(length(r)==1, r[1] , r[2])
  if(pointsDistance(c(x2, y2), c2) < r2){
    ep2 <- getPointInNodeCircle(ns2, ab[2], c2, r2)
  }else{
    ep2 <- getMirrorPoint(ns2, ab[2])
  }
  # points(ep1[1], ep1[2])
  # lines(c(x1, ep1[1]), c(y1, ep1[2]))
  # points(ep2[1], ep2[2])
  # lines(c(x2, ep2[1]), c(y2, ep2[2]))
  if(method[1]==2){
    coefs <- c(x1, y1,
               ep1[1], ep1[2],
               ep2[1], ep2[2],
               x2, y2)
  }else{
    coefs <- c(x1, y1, 
               x2, y2)
  }
  coefs <- matrix(coefs, ncol = 2, byrow = TRUE)
  xy <- bezier(coefs = coefs,
               center = matrix(c(cx, cy), ncol = 2),
               r_coefs = c(pointsDistance(c(x1, y1), c1),
                           pointsDistance(c(x2, y2), c2)),
               r = ifelse(length(r)==1, r, 0), w = w,
               evaluation = evaluation,
               method = method[1])
  return(xy)
}

archPlot <- function(x, y, r, angle=300, init.angle=0, length.out=100, bck = FALSE){
  gamma <- 2*pi * angle * (1:2) / 360 + init.angle * pi/180
  t2xy <- function(rx, t, x0, y0) list(x=rx*cos(t)+x0, y=rx*sin(t)+y0)
  P <- t2xy(r, seq.int(gamma[1], gamma[2], length.out = length.out), x, y)
  # The first 10 points and last 10 points use the mirror points
  if(length.out>60 && !bck){
    rr <- 15
    N <- .6
    for(i in seq.int(rr)){
      xy <- mirrorP(Q=c(P$x[i], P$y[i]),
                    P=c(P$x[rr+1], P$y[rr+1]),
                    R=c(P$x[rr+2], P$y[rr+2]),
                    N=N)
      P$x[i] <- xy[1]
      P$y[i] <- xy[2]
    }
    ll <- length.out - rr
    for(i in seq.int(length.out)[-seq.int(ll)]){
      xy <- mirrorP(Q=c(P$x[i], P$y[i]),
                    P=c(P$x[ll-1], P$y[ll-1]),
                    R=c(P$x[ll-2], P$y[ll-2]),
                    N=N)
      P$x[i] <- xy[1]
      P$y[i] <- xy[2]
    }
  }
  return(P)
}


##TODO: put the gene symbols in proper location.
getCutPos <- function(x, breaks, n){
  x[x<=min(breaks)] <- min(breaks)+.1
  x[x>max(breaks)] <- max(breaks)
  at <- as.numeric(cut(x,
                       breaks = breaks,
                       labels = seq.int(n)))
  return(at)
}
getSrt <- function(x, y){
  ## length of x, y must be 2
  srt <- ifelse(diff(x)[1]==0, 0, 180/pi*atan(diff(y)[1]/diff(x)[1]))
  if(is.na(srt)[1]) srt <- 0
  return(srt)
}
rotatePoint <- function(x1, y1, x2, y2){
  c(-(y2-y1)+x1, x2-x1+y1)
}
prettyMark <- function(x, u){
  div <- findInterval(x, c(0, 1e3, 1e6, 1e9, 1e12))
  un <- findInterval(u, c(0, 1e3, 1e6, 1e9, 1e12))
  digit <- u/10^(3*(un-1))
  paste0(round(x/10^(3*(div-1)), ifelse(digit==1, 0, 1)), 
         c("","K","M","B","T")[div])
}

#' @importFrom grDevices col2rgb rgb
invertCol <- function (col){
  n <- names(col)
  col <- col2rgb(col,alpha=FALSE)
  int_col <- abs(255 - col)
  int_col <- apply(int_col, 2, function(.ele){
    rgb(.ele[1 ], .ele[2 ], .ele[3 ], maxColorValue = 255)})
  int_col <- unlist(int_col)
  names(int_col) <- n
  int_col
}


resampleDataByFun <- function(fromGR, targetGR, FUN=viewMeans, dropZERO = TRUE,
                              ...){
  ## resample the range
  strand(targetGR) <- "*"
  seqn <- as.character(seqnames(targetGR))[1]
  stopifnot('Only one seqlevel is supported.'=
              all(as.character(seqnames(targetGR))==
                    seqn))
  stopifnot('score' %in% colnames(mcols(fromGR)))
  strand(fromGR) <- "*"
  fromGR <- subsetByOverlaps(fromGR, targetGR)
  cvg <- coverage(fromGR, weight = fromGR$score)
  x <- cvg[[seqn]]
  .vw <- Views(x, start = ranges(targetGR))
  targetGR$score <- FUN(.vw, ...)
  targetGR$score[is.na(targetGR$score)] <- 0
  if(dropZERO) targetGR <- targetGR[targetGR$score!=0]
  orderedGR(targetGR)
}

breakPointByBin <- function(x, y, start, end, seqn){
  if(length(x)>2){
    n <- length(x)
  }else{
    n <- 11
    x <- seq(x[1], x[2], length.out=n)
    y <- seq(y[1], y[2], length.out=n)
  }
  rg <- seq(start, end, length.out=n)
  GRanges(seqn, IRanges(start=rg[-n], end = rg[-1]),
          x0=x[-n], x1=x[-1],
          y0=y[-n], y1=y[-1])
}

findOlPos <- function(query, subject){
  stopifnot(all(width(query)==1))
  ol <- findOverlaps(query, subject, type='within')
  ol <- split(subjectHits(ol), queryHits(ol))
  ol <- vapply(ol, function(.ol){
    .ol[ceiling(length(.ol)/2)]
  }, numeric(1L))
  ol[order(as.numeric(names(ol)))]
}

calTickPos <- function(feature.tick, curve_gr, arrowLen){
  ol <- findOlPos(feature.tick, curve_gr)
  a <- start(feature.tick)[as.numeric(names(ol))]
  ratio <- (a - start(curve_gr)[ol])/width(curve_gr[ol])
  x1 <- curve_gr$x0[ol] + (curve_gr$x1-curve_gr$x0)[ol]*ratio
  y1 <- curve_gr$y0[ol] + (curve_gr$y1-curve_gr$y0)[ol]*ratio
  srt <- atan((curve_gr$y1-curve_gr$y0)[ol]/(curve_gr$x1-curve_gr$x0)[ol]) + pi/2
  srt[is.na(srt)] <- 0
  tickLen <- arrowLen/3
  x2 <- x1 + tickLen*cos(srt)
  y2 <- y1 + tickLen*sin(srt)
  x3 <- x1 + 2*tickLen*cos(srt)
  y3 <- y1 + 2*tickLen*sin(srt)
  list(x1=x1, y1=y1, x2=x2, y2=y2,
       x3=x3, y3=y3,
       srt=srt,
       id=as.numeric(names(ol)),
       ol=ol)
}

calGenePos <- function(fgf, curve_gr, arrowLen){
  fgf <- subsetByOverlaps(fgf, curve_gr)
  fgf <- sort(fgf)
  s <- e <- fgf
  rg <- range(curve_gr)
  missing_start <- start(s)<start(rg)
  missing_end <- end(e)>end(rg)
  start(s)[missing_start] <- start(rg)+1
  width(s) <- 1
  end(e)[missing_end] <- end(rg)-1
  start(e) <- end(e)
  ol_s <- calTickPos(s, curve_gr, arrowLen)
  ol_e <- calTickPos(e, curve_gr, arrowLen)
  f <- ifelse(as.character(strand(fgf))=='-',
              ol_e$ol, ol_s$ol)
  t <- ifelse(as.character(strand(fgf))=='-',
              ol_s$ol, ol_e$ol)
  via_points <- mapply(seq, ol_s$ol, ol_e$ol, SIMPLIFY = FALSE)
  via_points <- lapply(via_points, function(.ele){
    .ele[-unique(c(1, length(.ele)))]
  })
  via_points <- lapply(via_points, function(.ele){
    curve_gr[.ele]
  })
  x1 <- mapply(function(p1, ps, pn){
    if(length(ps)>0){
      c(p1, ps$x0, ps$x1[length(ps)], pn)
    }else{
      c(p1, pn)
    }
  }, ol_s$x1, via_points, ol_e$x1, SIMPLIFY = FALSE)
  y1 <- mapply(function(p1, ps, pn){
    if(length(ps)>0){
      c(p1, ps$y0, ps$y1[length(ps)], pn)
    }else{
      c(p1, pn)
    }
  }, ol_s$y1, via_points, ol_e$y1, SIMPLIFY = FALSE)
  
  # start points
  neg_strand <- as.character(strand(fgf))=='-'
  x2 <- ifelse(neg_strand,
               ol_e$x1, ol_s$x1)
  y2 <- ifelse(neg_strand,
               ol_e$y1, ol_s$y1)
  x3 <- ifelse(neg_strand,
               ol_e$x3, ol_s$x3)
  y3 <- ifelse(neg_strand,
               ol_e$y3, ol_s$y3)
  getFirst2EleDiff <- function(.ele, isNeg){
    if(isNeg){
      .ele[length(.ele) - 1]-.ele[length(.ele)]
    }else{
      .ele[2] - .ele[1]
    }
  }
  x0diff <- mapply(getFirst2EleDiff, x1, neg_strand)
  y0diff <- mapply(getFirst2EleDiff, y1, neg_strand)
  srt <- atan(y0diff/x0diff) + ifelse(x0diff<0, pi, 0)
  srt[is.na(srt)] <- 0
  x4 <- x3 + 1.25*arrowLen*cos(srt)
  y4 <- y3 + 1.25*arrowLen*sin(srt)
  
  list(xs=x1, ys=y1,
       x1=x2, y1=y2,
       x2=x3, y2=y3,
       x3=x4, y3=y4,
       srt=srt,
       fgf=fgf,
       x0diff=x0diff,
       missing_start=ifelse(neg_strand,missing_end,missing_start))
}

textWidth <- function(xlim, ...){ 
  convertWidth(grobWidth(...),
               unitTo = 'npc', valueOnly = TRUE)*diff(xlim) }
textHeight <- function(ylim, ...){
  convertHeight(grobHeight(...),
               unitTo = 'npc', valueOnly = TRUE)*diff(ylim) }
## TODO fix the algorithm by the center of the cluster.
safeTextCoor <- function(textCoor, x, y, tg, srt, xlim, ylim, logic=TRUE, force=6){
  left <- textCoor$x-textCoor$w/4
  right <- textCoor$x+textCoor$w/4
  bottom <- textCoor$y-textCoor$h/4
  top <- textCoor$y+textCoor$h/4
  tg.w <- textWidth(xlim, tg)/4
  tg.h <- textHeight(ylim, tg)/4
  l <- any(x+tg.w>left & x-tg.w<right & y+tg.h>bottom & y-tg.h<top)
  if(is.na(l)){
    l <- FALSE
  }
  if(logic){
    return(!l)
  }
  srt[is.na(srt)] <- 0
  for(i in seq.int(force)){
    if(!l){
      return(c(x, y))
    }else{
      x <- x+tg.w*cos(srt+pi/2)*sample(c(0, 1), 1)
      y <- y+tg.h*sin(srt+pi/2)*sample(c(0, 1), 1)
      l <- any(x+tg.w>left & x-tg.w<right & y+tg.h>bottom & y-tg.h<top)
      if(is.na(l)) {
        l <- FALSE
      }
    }
  }
  return(c(x, y))
}

plotBouquet <- function(pP, fgf, atacSig,
                        lwd.backbone, col.backbone,
                        lwd.maxAtacSig,
                        reverseATACSig,
                        col.backbone_background,
                        lwd.gene,
                        coor_mark_interval=1e5,
                        col.coor='black',
                        show_coor = TRUE,
                        coor_tick_unit = 1e3,
                        label_gene = TRUE,
                        col.tension_line = 'black',
                        lwd.tension_line = 1,
                        safe_text_force = 6,
                        arrowLen,xlim,ylim){
  textCoor <- data.frame(x=numeric(), y=numeric(),
                         w=numeric(), h=numeric())
  seqn <- as.character(seqnames(fgf[1]))
  curve_gr <- lapply(pP, function(.ele){
    .ele$seqn <- seqn
    do.call(breakPointByBin, .ele)
    })
  curve_gr <- unlist(GRangesList(curve_gr))
  missing_atacSig <- FALSE
  if(missing(atacSig)){
    missing_atacSig <- TRUE
  }else{
    if(length(atacSig)==0){
      missing_atacSig <- TRUE
    }
  }
  if(!missing_atacSig){
    if(is(atacSig, 'track')){
      if(atacSig$format=='WIG'){
        atacSig <- parseWIG(trackScore=atacSig,
                            chrom=seqn,
                            from=start(range(curve_gr)),
                            to=end(range(curve_gr)))
      }
      atacSig <- atacSig$dat
    }
    stopifnot(is(atacSig, 'GRanges'))
    stopifnot('score' %in% colnames(mcols(atacSig)))
    atacSig <- resampleDataByFun(atacSig, curve_gr, dropZERO = FALSE,
                                 na.rm = TRUE)
    atacSigScoreRange <- quantile(log2(atacSig$score+1),
                                  probs = c(.1, .99))
    if(atacSigScoreRange[1]==atacSigScoreRange[2]){
      atacSigScoreRange <- range(log2(atacSig$score+1))
    }
    if(atacSigScoreRange[1]!=atacSigScoreRange[2]){
      atacSigBreaks <- c(-1,
                         seq(atacSigScoreRange[1],
                             atacSigScoreRange[2],
                             length.out = lwd.maxAtacSig-1),
                         max(log2(atacSig$score+1))+1)
      atacSiglabels <- seq_along(atacSigBreaks)[-length(atacSigBreaks)]
      if(reverseATACSig){
        atacSiglabels <- rev(atacSiglabels)
      }
      atacSig$lwd <- as.numeric(as.character(
        cut(log2(atacSig$score+1),
            breaks=atacSigBreaks,
            labels = atacSiglabels)))
      ## add atac signals
      grid.segments(curve_gr$x0, curve_gr$y0,
                    curve_gr$x1, curve_gr$y1,
                    default.units = "native",
                    gp=gpar(lwd=lwd.backbone+atacSig$lwd,
                            col=col.backbone_background,
                            lineend=1))
    }
  }else{
    grid.lines(c(curve_gr$x0,curve_gr$x1[length(curve_gr)]),
               c(curve_gr$y0,curve_gr$y1[length(curve_gr)]),
               default.units = "native",
               gp=gpar(lwd=lwd.maxAtacSig/2,
                       lty=1,
                       col=col.backbone_background))
  }
  ## add backbone
  grid.lines(c(curve_gr$x0,curve_gr$x1[length(curve_gr)]),
             c(curve_gr$y0,curve_gr$y1[length(curve_gr)]),
             default.units = "native",
             gp=gpar(lwd=lwd.backbone,
                     lty=1,
                     col=1))
  
  ## add genomic coordinates
  if(show_coor){
    r_tick <- range(curve_gr)
    end(r_tick) <- ceiling(end(r_tick)/coor_tick_unit)*coor_tick_unit
    start(r_tick) <- floor(start(r_tick)/coor_tick_unit)*coor_tick_unit
    strand(r_tick) <- "*"
    feature.tick <- GenomicRanges::slidingWindows(r_tick, width = 1, step=coor_tick_unit)[[1]]
    feature.tick$col <- col.tension_line
    tick.xy <- calTickPos(feature.tick, curve_gr, arrowLen)
    grid.segments(tick.xy$x1, tick.xy$y1,
                  tick.xy$x2, tick.xy$y2,
                  default.units = "native",
                  gp=gpar(col=col.tension_line,
                          lwd=lwd.tension_line))
    if(coor_mark_interval){
      feature.tick.mark <- feature.tick[tick.xy$id]
      mark <- start(feature.tick.mark)/coor_mark_interval
      keep <- which(mark==round(mark))
      if(length(keep)>0){
        grid.segments(tick.xy$x1[keep], tick.xy$y1[keep],
                 tick.xy$x3[keep], tick.xy$y3[keep],
                 default.units = "native",
                 gp=gpar(col=col.coor,
                         lwd=lwd.tension_line))
        for(k in keep){
          lab <- prettyMark(start(feature.tick.mark[k]),
                            coor_mark_interval)
          tg <- textGrob(label=lab,
                         x=tick.xy$x3[k], y=tick.xy$y3[k],
                         default.units = "native",
                         gp=gpar(col=col.coor),
                         just=c(.5, -.2),
                         rot=180*tick.xy$srt[k]/pi-90)
          if(safeTextCoor(textCoor,
                          x=tick.xy$x3[k],
                          y=tick.xy$y3[k],
                          tg=tg,
                          xlim=xlim,
                          ylim=ylim,
                          logic=TRUE,
                          force = safe_text_force)){
            grid.draw(tg)
            textCoor <- rbind(textCoor,
                              data.frame(x = tick.xy$x3[k],
                                         y = tick.xy$y3[k],
                                         w = textWidth(tg, xlim=xlim),
                                         h = textHeight(tg, ylim=ylim)))
          }
        }
      }
    }
  }
  ## add gene annotation
  genePos <- calGenePos(fgf, curve_gr, arrowLen)
  null <- mapply(function(x, y, col, lwd){
    grid.lines(x, y,
               default.units = "native",
               gp=gpar(col=col, lwd=lwd))
    }, x=genePos$xs, y=genePos$ys, col=genePos$fgf$col, lwd=lwd.gene)
  grid.segments(x0=genePos$x1, y0=genePos$y1,
                x1=genePos$x2, y1=genePos$y2,
                default.units = "native",
                gp=gpar(col=genePos$fgf$col))
  isGene <- genePos$fgf$type %in% 'gene' & !genePos$missing_start
  if(any(isGene)){
    grid.segments(x0=genePos$x2[isGene],
                  x1=genePos$x3[isGene],
                  y0=genePos$y2[isGene],
                  y1=genePos$y3[isGene],
                  default.units = "native",
                  gp=gpar(col=genePos$fgf$col[isGene],
                          fill=genePos$fgf$col[isGene]),
                  arrow=arrow(angle = 15,
                              type='closed',
                              length=unit(arrowLen, units = 'inch')))
  }
  isRE <- genePos$fgf$type %in% 'cRE'
  if(any(isRE)){
    grid.points(x = genePos$x1[isRE],
                y = genePos$y1[isRE],
                pch = 11,
                size = unit(0.25, "char"),
                default.units = "native",
                gp=gpar(col=genePos$fgf$col[isRE],
                        fill=genePos$fgf$col[isRE]))
  }
  for(k in seq_along(genePos$fgf)){
    if(is.na(genePos$x2[k])) next
    vadj <- -.2 #ifelse(genePos$x0diff[k]<0, 1.2, -.2)
    hadj <- 0
    srt <- 180*genePos$srt[k]/pi
    if(srt>90 && srt<270){
      srt <- srt - 180
      hadj <- 1
    }
    tg <- textGrob(label=genePos$fgf$label[k],
                    x=genePos$x2[k], y=genePos$y2[k],
                    default.units = "native",
                    just=c(hadj, vadj),
                    rot=srt,
                    gp=gpar(col=genePos$fgf$col[k],
                            cex=genePos$fgf$cex))
    lab.xy <- safeTextCoor(textCoor,
                           genePos$x2[k],
                           genePos$y2[k],
                           tg=tg,
                           xlim=xlim,
                           ylim=ylim,
                           srt=genePos$srt[k],
                           logic=FALSE,
                           force = safe_text_force)
    if(lab.xy[1]!=genePos$x2[k]){
      tg <- grid.text(label=genePos$fgf$label[k],
                      x=lab.xy[1], y=lab.xy[2],
                      default.units = "native",
                      just=c(hadj, vadj),
                      gp=gpar(col=genePos$fgf$col[k],
                              cex=genePos$fgf$cex))
      grid.segments(genePos$x3[k], genePos$y3[k],
               lab.xy[1],
               lab.xy[2],
               default.units = "native",
               gp=gpar(col = genePos$fgf$col[k], 
                       lwd = .5,
                       lty = 4))
    }else{
      grid.draw(tg)
    }
    textCoor <- rbind(textCoor,
                      data.frame(x = lab.xy[1],
                                 y = lab.xy[2],
                                 w = textWidth(tg, xlim=xlim),
                                 h = textHeight(tg, ylim=ylim)))
  }
  
}

