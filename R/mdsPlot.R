#' Plot genomic interactions by multi-dimensional scaling plot
#' @description This function will convert the interactions scores into a
#' distance matrix and then plot the matrix by multi-dimensional scaling plot.
#' @param gi An object of \link[InteractionSet:GInteractions-class]{GInteractions}
#' @param range The region to plot. an object of \link[GenomicRanges:GRanges-class]{GRanges}
#' @param feature.gr The annotation features to be added. An object of \link[GenomicRanges:GRanges-class]{GRanges}.
#' @param k The dimension of plot. 2: 2d, 3: 3d.
#' @param atacSig The ATAC-seq signals. An object of \link[GenomicRanges:GRanges-class]{GRanges} with scores or an object of \link{track}.
#' @param show_coor Plot ticks in the line to show the DNA compact tension.
#' @param reverseATACSig Plot the ATAC-seq signals in reverse values.
#' @param coor_tick_unit The bps for every ticks. Default is 1K.
#' @param coor_mark_interval The coordinates marker interval. Numeric(1). Set to 0
#' to turn it off. The default value 1e5 means show coordinates every 0.1M bp.
#' @param label_gene Show gene symbol or not.
#' @param lwd.backbone,lwd.gene,lwd.tension_line,lwd.maxAtacSig Line width for the 
#' linker, gene, interaction node circle, the dashed line of interaction edges, the tension line and the maximal reversed ATAC signal.
#' @param col.backbone,col.backbone_background,col.tension_line,col.coor Color
#' for the DNA chain, the compact DNA chain, the node circle, the linker, the tension line and the coordinates marker.
#' @param alpha.backbone_background Alpha channel for transparency of backbone background.
#' @param length.arrow Length of the edges of the arrow head (in inches).
#' @param safe_text_force The loops to avoid the text overlapping.
#' @param square A logical value that controls whether control points for the curve are created city-block fashion or obliquely. See \link[grid]{grid.curve}.
#' @param ... Parameter will be passed to \link[MASS]{isoMDS}.
#' @return Coordinates for 2d.
#' @importClassesFrom InteractionSet GInteractions
#' @importMethodsFrom InteractionSet regions anchorIds
#' @importFrom MASS isoMDS
#' @importFrom Matrix sparseMatrix
#' @importFrom stats as.dist splinefun
#' @import GenomicRanges
#' @import S4Vectors
#' @importFrom BiocGenerics sort
#' @export
#' @examples
#' library(InteractionSet) 
#' gi <- readRDS(system.file("extdata", "nij.chr6.51120000.53200000.gi.rds",
#'  package="trackViewer"))
#' range <- GRanges("chr6", IRanges(51120000, 53200000))
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
#' mdsPlot(gi, range, feature.gr)
mdsPlot <- function(gi, range, feature.gr, k=2,
                    atacSig,
                    lwd.backbone = 2, col.backbone = 'gray',
                    lwd.maxAtacSig = 8, reverseATACSig = TRUE,
                    col.backbone_background = 'gray70',
                    alpha.backbone_background = 0.5,
                    lwd.gene = 3,
                    coor_mark_interval = 5e5, col.coor = "black",
                    show_coor = TRUE,
                    coor_tick_unit = 5e4,
                    label_gene = TRUE,
                    col.tension_line = 'black',
                    lwd.tension_line = 1,
                    length.arrow = NULL,
                    safe_text_force = 3,
                    square = TRUE,
                    ...){
  gi <- checkGI(gi, fixedBin=TRUE)
  stopifnot(is.numeric(k))
  stopifnot(k==2 || k==3)
  stopifnot(is.numeric(coor_mark_interval))
  stopifnot(length(coor_mark_interval)==1)
  if(!missing(range)){
    stopifnot(is(range, "GRanges"))
    stopifnot('coor_tick_unit is too small.'=
                width(range(range))[1]/coor_tick_unit<1000)
    seqn <- as.character(seqnames(range)[1])
    gi <- subsetByOverlaps(gi, ranges = range, ignore.strand=TRUE)
    ol1 <- countOverlaps(first(gi), range)
    ol2 <- countOverlaps(second(gi), range)
    gi <- gi[ol1>0 & ol2>0]
  }else{
    seqn <- as.character(seqnames(first(gi))[1])
  }
  stopifnot('No interaction data available.'=length(gi)>0)
  x <- anchorIds(gi, type='first')
  y <- anchorIds(gi, type='second')
  dx <- diff(sort(unique(x)))
  dy <- diff(sort(unique(y)))
  message('Some region are missing from the input gi.'=
              all(dx==1) && all(dy==1) && identical(range(x), range(y)))
  ids <- sort(unique(c(x, y)))
  v <- mcols(gi)$score
  r <- regions(gi)[ids]
  x <- x - min(ids) + 1
  y <- y - min(ids) + 1
  l <- length(r)
  m <- sparseMatrix(x, y, x=v, dims = c(l, l),
                    dimnames = list(seq.int(l), seq.int(l)))
  m <- as.matrix(t(m))
  m <- cf2pd(m)
  m[upper.tri(m)] <- 0
  for(i in seq.int(l)){
    m[i, i] <- 0
  }
  d <- as.dist(m)
  mds <- isoMDS(d, k=k, ...)
  p <- r[as.numeric(rownames(mds$points))]
  mcols(p) <- mds$points
  colnames(mcols(p)) <- c('x', 'y', 'z')[seq.int(k)]
  view3dStructure(p=p, k=k,
                  feature.gr=feature.gr,
                  atacSig=atacSig,
                  lwd.backbone = lwd.backbone,
                  col.backbone = col.backbone,
                  lwd.maxAtacSig = lwd.maxAtacSig,
                  reverseATACSig = reverseATACSig,
                  col.backbone_background = col.backbone_background,
                  alpha.backbone_background = alpha.backbone_background,
                  lwd.gene = lwd.gene,
                  coor_mark_interval = coor_mark_interval,
                  col.coor = col.coor,
                  show_coor = show_coor,
                  coor_tick_unit = coor_tick_unit,
                  label_gene = label_gene,
                  col.tension_line = col.tension_line,
                  lwd.tension_line = lwd.tension_line,
                  length.arrow = length.arrow,
                  safe_text_force = safe_text_force,
                  square = square)
}

cf2pd <- function(m, alpha = -0.25, inf_dist){
  pd <- m ^ alpha
  if(missing(inf_dist)){
    inf_dist <- max(pd[!(is.na(pd)|is.infinite(pd))]) + 1
  }
  pd[is.na(pd) | is.infinite(pd)] <- inf_dist
  pd
}