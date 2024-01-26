#' plot GInteractions
#' @description plot graph for GInteractions
#' @param gi an object of \link[InteractionSet:GInteractions-class]{GInteractions}
#' @param range the region to plot. an object of \link[GenomicRanges:GRanges-class]{GRanges}
#' @param feature.gr the feature.gr to be added. an object of \link[GenomicRanges:GRanges-class]{GRanges}
#' @param label_region label the region start-end or not.
#' @param show_tension plot ticks in the line to show the DNA compact tension.
#' @param tension_unit the bps for every ticks. Default is 1K.
#' @param label_gene show gene symbol or not.
#' @param lwd.backbone,lwd.gene,lwd.nodeCircle,lwd.edge line width for the 
#' linker, gene, interaction node circle and the dashed line of interaction edges.
#' @param ... Not used.
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
#' feature.gr <- subsetByOverlaps(feature.gr, regions(gi))
#' symbols <- mget(feature.gr$gene_id, org.Hs.egSYMBOL, ifnotfound=NA)
#' feature.gr$label <- vapply(symbols, function(.ele) .ele[1], character(1L))
#' feature.gr$col <- sample(1:7, length(feature.gr), replace=TRUE)
#' feature.gr$type <- sample(c("cRE", "gene"), 
#'                          length(feature.gr), replace=TRUE, 
#'                          prob=c(0.1, 0.9))
#' plotGInteractions(gi, range, feature.gr)
plotGInteractions <- function(gi, range, feature.gr,
                              label_region=FALSE,
                              lwd.backbone = 3,
                              lwd.gene = 2,
                              lwd.nodeCircle = 1,
                              lwd.edge = 2,
                              show_tension = TRUE,
                              tension_unit = 1e3,
                              label_gene = TRUE,
                              ...){
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
  if(!all(seqnames(reg)==seqnames(reg)[1])){
    warning('All interaction must within one chromosome.')
  }
  ol <- findOverlaps(reg, drop.self=TRUE, drop.redundant=TRUE, minoverlap = 2)
  if(length(ol)>0){
    warning("There are overlaps in the input regions")
  }
  if(length(feature.gr$col)==0) feature.gr$col <- rep("black", length(feature.gr))
  if(length(feature.gr$type)==0) feature.gr$type <- rep("gene", length(feature.gr))
  stopifnot(length(feature.gr$type)==length(feature.gr))
  stopifnot(all(feature.gr$type %in% c('cRE', 'gene')))
  feature.gr[feature.gr$type=='cRE'] <- 
    promoters(feature.gr[feature.gr$type=='cRE'], upstream = 0, downstream = 1)
  r_tick <- range(reg)
  tension_unit <- 1e3
  end(r_tick) <- ceiling(end(r_tick)/tension_unit)*tension_unit
  start(r_tick) <- floor(start(r_tick)/tension_unit)*tension_unit
  strand(r_tick) <- "*"
  feature.tick <- GenomicRanges::slidingWindows(r_tick, width = 1, step=tension_unit)[[1]]
  feature.tick$col <- 'black'
  feature.tick$type <- 'tick'
  feature.gr <- c(feature.gr, feature.tick)
  feature.gr_reg <- subsetByOverlaps(feature.gr, reg, ignore.strand = TRUE)
  ol <- findOverlaps(feature.gr_reg, reg)
  feature.gr_reg <- split(feature.gr_reg[queryHits(ol)], subjectHits(ol))
  
  nodes <- unique(as.character(sort(c(anchorIds(gi, type="first"),
                                      anchorIds(gi, type="second")))))
  ## add genomic coordinates to edge
  d <- log10(distance(reg[nodes[-length(nodes)]],
                      reg[nodes[-1]]) + 1) + 1
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
  nodes <- data.frame(names = nodes,
                      size = (width(reg[nodes]))/min(width(reg[nodes])))
  gL <- graph_from_data_frame(d = edgeL_link, directed=FALSE, vertices = nodes)
  cl <- igraph::components(gL, mode = "weak")
  sGnodes <- split(names(cl$membership), cl$membership)
  g <- graph_from_data_frame(d = edgeL, directed=FALSE, vertices = nodes)
  layout <- layout_with_fr(g) ## only layout_with_fr and layout_with_kk work OK
  # plot(g, layout = layout, main = 'Fruchterman-Reingold force-directed algorithm')
  # axis(1); axis(2)
  
  nodeXY <- norm_coords(layout, -1, 1, -1, 1)
  rownames(nodeXY) <- nodes$names
  colnames(nodeXY) <- c("X", "Y")
  vertex.factor <- 72
  vertex.size <- 1/vertex.factor * V(g)$size
  nodeXY <- fixXY(nodeXY, vertex.size, edgeL_link, lwd = lwd.backbone/300)
  maxv <- max(vertex.size)
  xlim <- range(nodeXY[, 1])
  ylim <- range(nodeXY[, 2])
  d_xlim <- diff(xlim)/5
  d_ylim <- diff(ylim)/5
  xlim <- c(xlim[1] - d_xlim - maxv, xlim[2] + d_xlim + maxv)
  ylim <- c(ylim[1] - d_ylim - maxv, ylim[2] + d_ylim + maxv)
  # nodeX <- nodeXY[, "X"]
  # nodeY <- nodeXY[, "Y"]
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

  plot(0, 0, type="n", asp = 1, xlab="", ylab="", xaxt="n", yaxt="n", 
       xaxs="i", yaxs="i", xlim=xlim, ylim=ylim, frame.plot=FALSE,
       new=TRUE)
  sH <- strheight("0")
  arrowLen <- grid::convertUnit(grid::stringHeight("0"), unitTo = "inches", valueOnly = TRUE)/3
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
  ## plot the bacground circle
  for(i in seq_along(clusterCenter)){
    polygon(archPlot(clusterCenter[[i]]$x, clusterCenter[[i]]$y,
                     clusterCenter[[i]]$r,
                     angle=360, init.angle=0, length.out=100,
                     bck = TRUE),
            col = '#DDDDDD25', lty = 4, lwd=lwd.nodeCircle)
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
                                                        width = 10*tension_unit))
        this_features <- subsetByOverlaps(feature.gr, this_reg_gap,
                                          ignore.strand = TRUE)
        this_features_tick_density <- 
          length(this_features[this_features$type=="tick"])
        this_features_tick_density <-
          rescale(
            this_features_tick_density,
            to = c(1, 8),
            from = c(0, max(
              length(subsetByOverlaps(feature.gr[feature.gr$type=='tick'],
                                      range))/length(nodeShape),
              this_features_tick_density)))
        idx <- c(i, i+1)
        ab <- getStartConnectionPoints(nodeShape[idx])
        idx <- ifelse(ab[1]==1, length(nodeShape[[1]]$x), 1)
        plotStartEndPoints(nodeShape[[i]], idx, 
                           xlim = xlim, ylim = ylim,
                           arrowLen = arrowLen,
                           col="red", pch=16,
                           feature.gr = this_features,
                           lwd_thick_weight = this_features_tick_density + 1,
                           gr = this_reg_gap,
                           lwd.backbone = lwd.backbone,
                           lwd.gene = lwd.gene,
                           label_gene = label_gene)
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
      this_features <- GRanges()
      this_reg_gap <- GRanges()
      this_features_tick_density <- 0
      reg0 <- reg[rownames(nodeXY)[i]]
      reg1 <- reg[rownames(nodeXY)[i+1]]
      if(start(reg1)>end(reg0)){
        this_reg_gap <- GRanges(seqnames(reg0), IRanges(end(reg0),
                                                          start(reg1)))
        this_features <- subsetByOverlaps(feature.gr, this_reg_gap,
                                          ignore.strand = TRUE)
        this_features_tick_density <- 
          length(this_features[this_features$type=="tick"])
        this_features_tick_density <-
          rescale(
            this_features_tick_density,
            to = c(1, 8),
            from = c(0, max(
              length(subsetByOverlaps(feature.gr[feature.gr$type=='tick'],
                            range))/length(nodeShape),
              this_features_tick_density)))
      }
      curveMaker(nodeShape[[i]], 
                 nodeShape[[i+1]],
                 ab,
                 cn$x, cn$y, r=cn$r,
                 w=2*d[i]/median(d),
                 feature.gr = this_features,
                 lwd_thick_weight = this_features_tick_density + 1,
                 gr = this_reg_gap,
                 arrowLen = arrowLen,
                 lwd.backbone = lwd.backbone,
                 lwd.gene = lwd.gene,
                 plot_tick = show_tension,
                 label_gene = label_gene)
      if(ab[1]!=1){
        nodeShape[[i]] <- lapply(nodeShape[[i]], rev)
      }
      last_a <- ifelse(ab[2]==1, length(nodeShape[[i]]$x), 1)
    }
    if(i==length(nodeShape)){ ## mark the end point
      reg0 <- reg[rownames(nodeXY)[i]]
      this_reg_gap <- GRanges(seqnames(reg0), IRanges(end(reg0),
                                                      width = 10*tension_unit))
      this_features <- subsetByOverlaps(feature.gr, this_reg_gap,
                                        ignore.strand = TRUE)
      this_features_tick_density <- 
        length(this_features[this_features$type=="tick"])
      this_features_tick_density <-
        rescale(
          this_features_tick_density,
          to = c(1, 8),
          from = c(0, max(
            length(subsetByOverlaps(feature.gr[feature.gr$type=='tick'],
                                    range))/length(nodeShape),
            this_features_tick_density)))
      plotStartEndPoints(nodeShape[[i]], last_a, 
                         xlim = xlim, ylim = ylim, 
                         start = nodeShape[[i]]$y[1],
                         arrowLen = arrowLen,
                         col="orange", pch=15,
                         feature.gr = this_features,
                         lwd_thick_weight = this_features_tick_density + 1,
                         gr = this_reg_gap,
                         lwd.backbone = lwd.backbone,
                         lwd.gene = lwd.gene,
                         label_gene = label_gene)
      if(last_a==1){
        nodeShape[[i]] <- lapply(nodeShape[[i]], rev)
      }
    }
  }
  ## plot the edge lines
  segments(nodeXY[as.character(edgeL_link$from), "X"], nodeXY[as.character(edgeL_link$from), "Y"], 
           nodeXY[as.character(edgeL_link$to), "X"], nodeXY[as.character(edgeL_link$to), "Y"],
           lwd=lwd.edge, lty=2, col = "gray80")
  ## 
  for(i in seq.int(nrow(nodeXY))){
    ## plot the shape of node
    lines(nodeShape[[i]]$x, nodeShape[[i]]$y,
          lwd=lwd.backbone, lty=1, col = 'darkblue')
    if(length(feature.gr_reg[[rownames(nodeXY)[i]]])>0){
      addGeneInNode(feature.gr_reg[[rownames(nodeXY)[i]]], nodeShape[[i]], reg[rownames(nodeXY)[i]], arrowLen=arrowLen, lwd.gene = lwd.gene,
                    plot_tick = show_tension,
                    label_gene = label_gene)
    }
    ## write the start and end position of the node
    if(label_region){
      ## plot the node dot
      points(nodeXY, pch=16, col="white", cex=nodesSize/strheight("0"), add=TRUE)
      ## label the node
      text(nodeXY, labels = rownames(nodeXY))
    }
  }
  
  return(invisible(g))
}

ABinSameSide <- function(A, B, c1, c2){
  return(
    ((c1[2] - c2[2])*(A[1] - c1[1]) +
       (c2[1] - c1[1])*(A[2] - c1[2])) *
      ((c1[2] - c2[2])*(B[1] - c1[1]) +
         (c2[1] - c1[1])*(B[2] - c1[2]))>0)
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
plotConnectionLines <- function(xy, lwd_thick_weight=1, lwd.backbone=3, ...){
  lines(xy$x[-c(1, length(xy$x))], xy$y[-c(1, length(xy$y))], 
        col="gray", lty=1, lwd=lwd.backbone+1+lwd_thick_weight)
  lines(xy,
        col="cyan", lty=1, lwd=lwd.backbone)
}

safeEndPoint <- function(x, lim, dx){
  if(x<lim[1]) x <- lim[1] + dx
  if(x>lim[2]) x <- lim[2] - dx
  x
}

plotStartEndPoints <- function(xy, idx, xlim, ylim, start, arrowLen, 
                               feature.gr,
                               lwd_thick_weight,
                               gr,
                               lwd.backbone,
                               lwd.gene, 
                               label_gene=TRUE,
                               ...){
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
          y1 <- y0 <- y + ddx
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

  plotConnectionLines(list(x=c(x0, x), 
                           y=c(y0, y)),
                      lwd_thick_weight=lwd_thick_weight,
                      lwd.backbone = lwd.backbone)
  if(length(feature.gr)>0){
    addGeneInNode(gs = feature.gr, nShape = list(x=c(x0, x), 
                                                 y=c(y0, y)), 
                  gr = gr,
                  arrowLen=arrowLen, lwd.gene = lwd.gene,
                  plot_tick = FALSE, label_gene = label_gene)
  }
  
  points(x0, y0, ...)
  if(!missing(start)){
    arrows(x1, y1, x0, y0, length = arrowLen)
  }else{
    arrows(x0, y0, x1, y1, length = arrowLen)
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
  size <- apply(edgeL_link, 1, function(.ele){
    sum(nodeXYout[.ele[c(1, 2)], 3])
  })
  edgeL_link[order(size, decreasing = TRUE), ]
}
fixXY <- function(nodeXY, vertex.size, edgeL_link, lwd = 10/300){
  ## rearrange nodes to make sure they are not overlap
  ## and make the interaction of nodes direct touch each other
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
mirrorP <- function(Q, P, R, N) {
  I <- matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2)
  n <- R - P
  d2 <- 1/sum(n^2)
  as.numeric(Q + 2*N * ( I - d2 * n %*% t(n) ) %*% ( P - Q ))
}

bezier <- function(coefs, center, r, evaluation=100, w){
  stopifnot('Dimention of coefs need to be 2'=length(dim(coefs))==2)
  stopifnot('Dimention of center need to be 2'=length(dim(center))==2)
  # mirror center
  center <- colMeans(center)
  # N is the strength of the bezier curve
  # strength is determined by the distance from the center of (p1, p2) to
  # (cx, cy)
  if(r==0){
    N <- w#2
  }else{
    N <- sqrt((mean(coefs[, 1]) - center[1])^2 +
                (mean(coefs[, 2]) - center[2])^2)
    N <- r/N
  }
  #n2 <- floor(nrow(coefs)/2)
  center_mirror <- mirrorP(center, coefs[1, ], coefs[nrow(coefs), ], N=N)
  x <- c(coefs[1, 1],
         center_mirror[1],
         coefs[nrow(coefs), 1])
  y <- c(coefs[1, 2],
         center_mirror[2],
         coefs[nrow(coefs), 2])
  n <- length(x)
  X <- Y <- single(evaluation)
  Z <- seq(0, 1, length = evaluation)
  X[1] <- x[1]
  X[evaluation] <- x[n]
  Y[1] <- y[1]
  Y[evaluation] <- y[n]
  for (i in 2:(evaluation - 1)) {
    z <- Z[i]
    xz <- yz <- 0
    const <- (1 - z)^(n - 1)
    for (j in 0:(n - 1)) {
      xz <- xz + const * x[j + 1]
      yz <- yz + const * y[j + 1]
      const <- const * (n - 1 - j)/(j + 1) * z/(1 - z)
    }
    X[i] <- xz
    Y[i] <- yz
  }
  list(x = as.numeric(X), y = as.numeric(Y))
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
                       feature.gr, gr, arrowLen, lwd_thick_weight = 1,
                       lwd.backbone = 3, lwd.gene = 2,
                       plot_tick = FALSE, label_gene = TRUE,
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
  # if(pointsDistance(c(x1, y1), c1) < r[1]){
  #   # find the extension point on the circle
  #   ep1 <- getPointInNodeCircle(ns1, ab[1], c1, r[1])
  # }else{
  #   ep1 <- getMirrorPoint(ns1, ab[1])
  # }
  # r2 <- ifelse(length(r)==1, r[1] , r[2])
  # if(pointsDistance(c(x2, y2), c2) < r2){
  #   ep2 <- getPointInNodeCircle(ns2, ab[2], c2, r2)
  # }else{
  #   ep2 <- getMirrorPoint(ns2, ab[2])
  # }
  # points(ep1[1], ep1[2])
  # lines(c(x1, ep1[1]), c(y1, ep1[2]))
  # points(ep2[1], ep2[2])
  # lines(c(x2, ep2[1]), c(y2, ep2[2]))
  # coefs <- c(x1, y1, 
  #            ep1[1], ep1[2],
  #            x2, y2,
  #            ep2[1], ep2[2])
  coefs <- c(x1, y1, 
             x2, y2)
  coefs <- matrix(coefs, ncol = 2, byrow = TRUE)
  xy <- bezier(coefs = coefs,
               center = matrix(c(cx, cy), ncol = 2),
               r = ifelse(length(r)==1, r, 0), w = w,
               evaluation = evaluation)
  plotConnectionLines(xy, lwd_thick_weight=lwd_thick_weight,
                      lwd.backbone = lwd.backbone, ...)
  if(length(feature.gr)>0){
    addGeneInNode(gs = feature.gr, nShape = xy, gr = gr,
                  arrowLen=arrowLen, lwd.gene = lwd.gene,
                  plot_tick = plot_tick,
                  label_gene = label_gene)
  }
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

# plogisMaker <- function(x1, y1, x2, y2, ...){
#   x <- NULL
#   if(abs(x1-x2)<10){
#     segments(x1, y1, x2, y2, ...)
#   }else{
#     if(x1<x2){
#       curve( plogis( x, scale = exp(round(log10(mean(c(x1, x2)))-1)), location = (x1 + x2) /2 ) * (y2-y1) + y1, 
#              x1, x2, add = TRUE, ...)
#     }else{
#       curve( plogis( x, scale = exp(round(log10(mean(c(x1, x2)))-1)), location = (x1 + x2) /2 ) * (y1-y2) + y2, 
#              x2, x1, add = TRUE, ...)
#     }
#   }
# }


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
}
rotatePoint <- function(x1, y1, x2, y2){
  c(-(y2-y1)+x1, x2-x1+y1)
}
addGeneInNode <- function(gs, nShape, gr, arrowLen, lwd.gene,
                          plot_tick = FALSE, label_gene = TRUE){
  gs <- split(gs, gs$type)
  n <- length(nShape$x)
  gr.cut <- seq(start(gr)-1, end(gr), length.out = n+1)
  mapply(function(gs.s, type){
    if(length(gs.s)>0){
      for(i in seq_along(gs.s)){
        at <- getCutPos(c(start(gs.s[i]), end(gs.s[i])),
                        breaks = gr.cut, n = n)
        if(at[1]==at[2]){
          if(at[1]==n){
            at[1] <- n -1
          }else{
            at[2] <- at[1] + 1
          }
        }
        A.x <- nShape$x[seq(at[1], at[2])]
        A.y <- nShape$y[seq(at[1], at[2])]
        switch(type,
               "tick" ={
                 if(plot_tick){
                   B <- rotatePoint(A.x[1], A.y[1], A.x[2], A.y[2])
                   segments(A.x[1], A.y[1], B[1], B[2], col = 'black')
                   # text(labels=paste0(round(start(gs.s[i])/1e6, digits = 1), 'M'),
                   #      x=A.x[1], y=A.y[1], adj=c(.5, 1.2), srt = getSrt(A.x, A.y))
                 }
               },
               "cRE"={
                 points(A.x[1], A.y[1], pch=11, col = gs.s$col[i])
               }, 
               "gene"={
                 lines(A.x, A.y, lwd=lwd.gene, col=gs.s$col[i])
                 pro <- promoters(gs.s[i], upstream = 0, downstream = 1)
                 pro <- subsetByOverlaps(pro, gr, type = 'within')
                 if(length(pro)>0){
                   pro <- promoters(gs.s[i], upstream = 0,
                                    downstream = ceiling(width(gr)/n))
                   at <- getCutPos(c(end(pro), start(pro)),
                                   breaks = gr.cut, n = n)
                   if(at[1]==at[2]){
                     at[1] <- at[2] + 1
                   }
                   if(as.character(strand(pro))=='-'){
                     at <- n - at
                   }
                   B.x <- nShape$x[at]
                   B.y <- nShape$y[at]
                   pD <- pointsDistance(c(B.x[1], B.y[1]), c(B.x[2], B.y[2]))
                   if(!is.na(pD)){
                     if(pD < arrowLen){
                       k <- ceiling(arrowLen/pD)
                       if(as.character(strand(pro))=='-'){
                         at[2] <- at[1] + k
                         if(at[2]>n) at[2] <- n
                       }else{
                         at[1] <- at[2] + k
                         if(at[1]>n) at[1] <- n
                       }
                       B.x <- nShape$x[at]
                       B.y <- nShape$y[at]
                     }
                   if(label_gene){
                     B1 <- rotatePoint(B.x[2], B.y[2], B.x[1], B.y[1])
                     segments(B.x[2], B.y[2], B1[1], B1[2], col = gs.s$col[i])
                     B2 <- rotatePoint(B1[1], B1[2], B.x[2], B.y[2])
                     arrows(B1[1], B1[2], B2[1], B2[2], length = arrowLen, col = gs.s$col[i], angle=15, code = 2)
                     text(labels=gs.s$label[i], B1[1], B1[2],
                          adj = c(ifelse(B1[1]>B2[1], 1, 0),
                                  ifelse(B.y[2]>B1[2], 1.5, -.5)),
                          srt = getSrt(B.x, B.y),
                          col = gs.s$col[i])
                   }else{
                     arrows(B.x[1], B.y[1], B.x[2], B.y[2], length = arrowLen, col = invertCol(gs.s$col[i]), angle=15, code = 1)
                   }
                   }
                 }
               })
      }
    }
  }, gs, names(gs))
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
