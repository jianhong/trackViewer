#' Plot xyz data in 2d or 3d
#' @description
#' Plot xyz data with grid or rgl package.
#' @param p GRanges object with mcols x, y, or with z
#' @param k The dimension of plot. 2: 2d, 3: 3d.
#' @param feature.gr The annotation features to be added. An object of \link[GenomicRanges:GRanges-class]{GRanges}.
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
#' @param ... Not used.
#' @return Coordinates for 2d.
#' @export
#' @examples
#' p <- readRDS(system.file('extdata', '4DNFI1UEG1HD.chr21.FLAMINGO.res.rds',
#'     package='trackViewer'))
#' feature.gr <- readRDS(system.file('extdata', '4DNFI1UEG1HD.feature.gr.rds',
#'     package='trackViewer'))
#' view3dStructure(p, k=3, feature.gr=feature.gr, length.arrow=unit(0.000006, 'native'))
view3dStructure <- function(p, k=3, feature.gr,
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
  stopifnot(is(p, 'GRanges'))
  stopifnot(
    'Only work for data in a single chromsome.'=
      all(as.character(seqnames(p))==as.character(seqnames(p)[1])))
  seqn <- as.character(seqnames(p)[1])
  stopifnot(is.numeric(k))
  stopifnot(k==2 || k==3)
  stopifnot(all(c('x', 'y') %in% colnames(mcols(p))))
  if(k==3){
    stopifnot('z' %in% colnames(mcols(p)))
  }
  feature.gr <- parseFeature(feature.gr=feature.gr)
  xlim <- range(p$x)
  ylim <- range(p$y)
  d_xlim <- diff(xlim)/10
  d_ylim <- diff(ylim)/10
  xlim <- c(xlim[1] - d_xlim, xlim[2] + d_xlim)
  ylim <- c(ylim[1] - d_ylim, ylim[2] + d_ylim)
  if(k==2){
    ## plot 2D
    grid.newpage()
    vp <- viewport(default.units = 'native', xscale=xlim, yscale=ylim)
    pushViewport(vp)
    on.exit({
      popViewport()
    })
    ## fix the arrow length
    if(!is.numeric(length.arrow)){
      arrowLen <- grid::convertUnit(grid::stringHeight("0"),
                                    unitTo = "inch",
                                    valueOnly = FALSE)
    }else{
      arrowLen <- length.arrow[1]
      stopifnot(is(arrowLen, 'unit'))
    }
    rate <- max(abs(
      c(convertWidth(unit(diff(xlim), 'native'),
                     unitTo = "inch", valueOnly = TRUE),
        convertHeight(unit(diff(ylim), 'native'),
                      unitTo = "inch", valueOnly = TRUE))))/
      convertUnit(arrowLen, unitTo = "inches", valueOnly = TRUE)
    
    ncurve <- length(p) - 1
    x1 <- p$x[-length(p$x)]
    x2 <- p$x[-1]
    y1 <- p$y[-length(p$y)]
    y2 <- p$y[-1]
    if(square){
      cps <- calcSquareControlPoints(
        x1, y1, x2, y2,
        curvature = 1, angle = 90, ncp = 1)
      shape <- interleave(ncp = 1,
                          ncurve = ncurve,
                          val = 0.5,
                          sval = 1, 
                          eval = 1,
                          e = cps$end)
      ncp <- 2
    }else{
      cps <- calcControlPoints(
        x1, y1, x2, y2,
        curvature = 1, angle = 90, ncp = 1)
      ncp <- 1
      shape <- rep(0.5, ncurve)
    }
    idset <- seq.int(ncurve)
    splineGrob <- xsplineGrob(c(x1, cps$x, x2), 
                              c(y1, cps$y, y2),
                              id = c(idset, rep(idset,  each = ncp), idset),
                              default.units = "native", 
                              shape = c(rep(0, ncurve), shape, rep(0, ncurve)), 
                              arrow = NULL, open = TRUE, name = "xspline")
    coor <- grobCoords(splineGrob)
    coor <- lapply(coor, function(.ele){
      list(x = convertX(unit(.ele$x, 'inch'),
                        unitTo = 'native', valueOnly = TRUE),
           y = convertY(unit(.ele$y, 'inch'),
                        unitTo = 'native', valueOnly = TRUE))
    })
    midPoints <- ceiling((start(p) + end(p))/2)
    midPoints <- c(start(p)[1],
                   midPoints[-c(1, length(midPoints))],
                   end(p)[length(p)])
    pP <- mapply(coor, midPoints[-length(midPoints)]+1, midPoints[-1],
                 FUN=function(xys, startP, endP){
                   list(x=xys$x, y=xys$y, start=startP, end=endP)
                 }, SIMPLIFY = FALSE)
    objCoor <- plotBouquet(pP=pP,
                           fgf=feature.gr,
                           atacSig=atacSig,
                           lwd.backbone=lwd.backbone,
                           col.backbone=col.backbone,
                           lwd.maxAtacSig = lwd.maxAtacSig,
                           reverseATACSig = reverseATACSig,
                           col.backbone_background=col.backbone_background,
                           alpha.backbone_background=alpha.backbone_background,
                           lwd.gene=lwd.gene,
                           coor_mark_interval=coor_mark_interval,
                           col.coor=col.coor,
                           show_coor = show_coor,
                           coor_tick_unit = coor_tick_unit,
                           label_gene = label_gene,
                           col.tension_line = col.tension_line,
                           lwd.tension_line = lwd.tension_line,
                           safe_text_force = safe_text_force,
                           arrowLen = arrowLen, rate = rate,
                           xlim=xlim, ylim=ylim)
  }else{
    ## plot 3D
    if(!requireNamespace('rgl', quietly = TRUE)){
      stop('rgl package is required')
    }
    t <- seq_along(p)
    ## spline smooth for each bin with 100 points
    resolusion <- 100
    tt <- seq(1, length(p), len = resolusion*length(p)+1)
    sdata <- lapply(colnames(mcols(p)), function(j){
      splinefun(t, mcols(p)[, j, drop=TRUE])(tt)
    })
    sdata <- do.call(cbind, sdata)
    sdata <- data.frame(sdata)
    colnames(sdata) <- colnames(mcols(p))
    p <- tile(p, n=resolusion)
    p <- unlist(p)
    stopifnot('Unexpected happend.'=length(p)+1==nrow(sdata))
    sdata0 <- sdata[-nrow(sdata), , drop=FALSE]
    sdata1 <- sdata[-1, , drop=FALSE]
    colnames(sdata0) <- paste0(colnames(sdata), '0')
    colnames(sdata1) <- paste0(colnames(sdata), '1')
    mcols(p) <- cbind(sdata0, sdata1)
    rgl::open3d()
    
    ## fix the arrow length
    if(!is.numeric(length.arrow)){
      arrowLen <- grid::convertUnit(grid::stringHeight("0"),
                                    unitTo = "inch",
                                    valueOnly = FALSE)
    }else{
      arrowLen <- length.arrow[1]
      stopifnot(is(arrowLen, 'unit'))
    }
    rate <- 50
    
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
                              from=start(range(p)),
                              to=end(range(p)))
        }
        atacSig <- atacSig$dat
      }
      stopifnot(is(atacSig, 'GRanges'))
      stopifnot('score' %in% colnames(mcols(atacSig)))
      atacSig <- resampleDataByFun(atacSig, p, dropZERO = FALSE,
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
        ## add signal to -z axis
        px <- as.numeric(t(as.matrix(mcols(p)[, c('x0', 'x1')])))
        py <- as.numeric(t(as.matrix(mcols(p)[, c('y0', 'y1')])))
        pz <- as.numeric(t(as.matrix(mcols(p)[, c('z0', 'z1')]) -
                             data.frame(z0=atacSig$lwd/rate,
                                        z1=0)))
        # polygon3d is very slow
        rgl::lines3d(px, py, pz,
                     col=col.backbone_background,
                     alpha = alpha.backbone_background,
                     tag='atac_signal')
      }
    }
    
    ## add backbone
    rgl::lines3d(x = c(p$x0, p$x1[length(p)]),
                 y = c(p$y0, p$y1[length(p)]),
                 z = c(p$z0, p$z1[length(p)]),
                 lwd=lwd.backbone,
                 lty=1,
                 col=col.backbone,
                 tag='backbone')
    
    ## add genomic coordinates
    if(show_coor){
      r_tick <- range(p)
      end(r_tick) <- ceiling(end(r_tick)/coor_tick_unit)*coor_tick_unit
      start(r_tick) <- floor(start(r_tick)/coor_tick_unit)*coor_tick_unit
      strand(r_tick) <- "*"
      feature.tick <- GenomicRanges::slidingWindows(r_tick, width = 1, step=coor_tick_unit)[[1]]
      feature.tick$col <- col.tension_line
      tick.xy <- calTickPos(feature.tick, p, arrowLen, rate=rate, kd=k)
      null <- mapply(function(x0, y0, z0, x1, y1, z1){
        rgl::segments3d(c(x0, x1), c(y0, y1), c(z0, z1),
                        col=col.tension_line,
                        lwd=lwd.tension_line,
                        tag='tick_minor')
      }, x0=tick.xy$x1, y0=tick.xy$y1, z0=tick.xy$z1,
      x1=tick.xy$x2, y1=tick.xy$y2, z1=tick.xy$z2)
      
      if(coor_mark_interval){
        feature.tick.mark <- feature.tick[tick.xy$id]
        mark <- start(feature.tick.mark)/coor_mark_interval
        keep <- which(mark==round(mark) & !is.na(tick.xy$x3) & 
                        !is.na(tick.xy$y3) & !is.na(tick.xy$z3))
        if(length(keep)>0){
          null <- mapply(function(x0, y0, z0, x1, y1, z1){
            rgl::segments3d(c(x0, x1), c(y0, y1), c(z0, z1),
                            col=col.coor,
                            lwd=lwd.tension_line,
                            tag='tick_major')
          }, x0=tick.xy$x1[keep], y0=tick.xy$y1[keep], z0=tick.xy$z1[keep],
          x1=tick.xy$x3[keep], y1=tick.xy$y3[keep], z1=tick.xy$z3[keep])
          coor_text <- prettyMark(start(feature.tick.mark[keep]),
                                  coor_mark_interval)
          rgl::text3d(tick.xy$x3[keep], tick.xy$y3[keep], tick.xy$z3[keep],
                      texts = coor_text,
                      id = coor_text,
                      col=col.coor, tag='tick_labels',
                      pos = 3)
        }
      }
    }
    
    ## add gene annotation
    genePos <- calGenePos(feature.gr, p, arrowLen, rate=rate, kd=3)
    if(length(genePos)>0){
      null <- mapply(function(x, y, z, col, lwd, id){
        rgl::lines3d(x, y, z,
                     col=col,
                     lwd=lwd,
                     id=id,
                     tag='gene_body')
      }, x=genePos$xs, y=genePos$ys, z=genePos$zs, 
      col=genePos$fgf$col, lwd=lwd.gene, id=genePos$fgf$label)
      ## add a vertical line at the TSS.
      null <- mapply(function(x0, y0, z0, x1, y1, z1, col, lwd, id){
        rgl::segments3d(c(x0, x1), c(y0, y1), c(z0, z1),
                        col=col, lwd=lwd, id=id, tag='tss_labels')
      }, x0=genePos$x1, y0=genePos$y1, z0=genePos$z1,
      x1=genePos$x2, y1=genePos$y2, z1=genePos$z2,
      col=genePos$fgf$col, lwd=lwd.gene/2, 
      id=genePos$fgf$label)
      
      isGene <- genePos$fgf$type %in% 'gene' & !genePos$missing_start
      if(any(isGene)){
        ## add arrow
        null <- mapply(function(x0, y0, z0, x1, y1, z1, col, fill, id){
          rgl::arrow3d(c(x0, y0, z0),
                       c(tail(x1, n=1), tail(y1, n=1), tail(z1, n=1)),
                       type='rotation',
                       id=paste0(id, '_arrow'),
                       col=col, fill=fill,
                       tag='tss_labels')
        }, x0=genePos$x2[isGene], y0=genePos$y2[isGene], z0=genePos$z2[isGene],
        x1=genePos$x3[isGene], y1=genePos$y3[isGene], z1=genePos$z2[isGene], # keep same z
        col=genePos$fgf$col[isGene], fill=genePos$fgf$col[isGene],
        id=genePos$fgf$label[isGene])
      }
      
      notGene <- (!genePos$fgf$type %in% 'gene') & (!genePos$missing_start)
      if(any(notGene)){
        rgl::points3d(x = genePos$x2[notGene],
                      y = genePos$y2[notGene],
                      z = genePos$z2[notGene],
                      pch = genePos$fgf$pch[notGene],
                      size = genePos$fgf$size[notGene],
                      col=genePos$fgf$col[notGene],
                      fill=genePos$fgf$col[notGene],
                      tag='cRE')
      }
      if(label_gene){
        rgl::text3d(genePos$x2, genePos$y2, (genePos$z2 + genePos$z1)/2,
                    texts = genePos$fgf$label,
                    id = genePos$fgf$label,
                    col=genePos$fgf$col, tag='gene_labels',
                    pos = 4)
      }
    }
    rgl::rgl.bringtotop()
  }
}