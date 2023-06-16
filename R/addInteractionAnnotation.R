#' Add annotation markers to the figure at a given position
#' @description A function to add annotation markers for emphasizing interactions 
#' @param obj A \link[InteractionSet:GInteractions-class]{GInteractions} object,
#' \link[GenomicRanges:GRanges-class]{GRanges} object or numeric vector.
#' For numeric vector, the positive value will generate a line with slope 1 and
#' negative value will generate a line at the position with slope -1.
#' @param idx The layer number of track.
#' @param FUN Function for plot. Available functions are 
#' \link[grid:grid.polygon]{grid.polygon}, \link[grid:grid.lines]{grid.lines},
#' and \link[grid:grid.text]{grid.text} for GInteractions object;
#' \link[grid:grid.lines]{grid.lines},
#' and \link[grid:grid.text]{grid.text} for GRanges object;
#' FUN is not used for numeric vector.
#' @param panel Plot regions. Available values are "top", "bottom".
#' @param ... Parameters will be passed to FUN.
#' @return invisible viewport for plot region.
#' @import grid
#' @export
#' @seealso See Also as \code{\link{addGuideLine}}, \code{\link{addArrowMark}}
#' @examples 
#' library(trackViewer)
#' library(InteractionSet)
#' gi <- readRDS(system.file("extdata", "nij.chr6.51120000.53200000.gi.rds",
#'  package="trackViewer"))
#' tads <- GInteractions(
#' GRanges("chr6", 
#'         IRanges(c(51130001, 51130001, 51450001, 52210001), width = 20000)),
#' GRanges("chr6", 
#'         IRanges(c(51530001, 52170001, 52210001, 53210001), width = 20000)))
#' range <- GRanges("chr6", IRanges(51120000, 53200000))
#' tr <- gi2track(gi)
#' viewTracks(trackList(tr), 
#'            gr=range, autoOptimizeStyle = TRUE)
## add TAD information
#' addInteractionAnnotation(tads, "tr", grid.lines,
#'  gp=gpar(col = "#E69F00", lwd=3, lty=3))
#' 
addInteractionAnnotation <- function(obj, idx, FUN=grid.polygon,
                                     panel=c("top", "bottom"), ...){
  if(!inherits(idx, c("numeric", "integer", "character")))
    stop("idx is required as a integer indexes of tracks")
  stopifnot(inherits(obj, c("GInteractions", "GRanges", "numeric")))
  if(!all(panel %in% c("top", "bottom"))){
    stop("panel must be 'top' and/or 'bottom'")
  }
  if(!is.numeric(obj)){
    stopifnot(is.function(FUN))
    FUN.NAME <- deparse(substitute(FUN))
    if(is(obj, "GInteractions")){
      stopifnot("FUN must be one of grid.polygon, grid.lines, or grid.text"=
                  FUN.NAME %in% c("grid.polygon", "grid.lines", "grid.text"))
    }else{#GRanges
      stopifnot("FUN must be one of grid.lines, or grid.text"=
                  FUN.NAME %in% c("grid.lines", "grid.text"))
    }
  }else{
    message("FUN is not supported for numeric input.")
  }
  
  dots <- list(...)
  dots <- dots[!names(dots) %in% c("x", "y", "default.units")]
  lastTrackViewer <- getOption("LastTrackViewer")
  if(length(lastTrackViewer)<0) {
    stop("Can not locate the last trackViewer window.")
  }
  xscale <- sort(lastTrackViewer$xscale)
  yscales <- lastTrackViewer$yscales
  .dat2target <- lastTrackViewer$back2back["target2", , drop=TRUE]>0
  vp <- lastTrackViewer$vp
  if(is.character(idx)){
    if(!all(idx %in% names(lastTrackViewer$yHeights))){
      stop("Can not find the input idx from laste trackViewer window.")
    }
    idx <- which(names(lastTrackViewer$yHeights) %in% idx)
  }
  yHeightsLimits2 <- yHeights <- cumsum(lastTrackViewer$yHeights)
  yHeightsLimits1 <- c(0, yHeights)[seq_along(yHeights)]
  names(yHeightsLimits1) <- names(yHeights)
  ym <- (xscale[2]-xscale[1] + 1)/2
  plotInteractionAnno <- function(anchor1, anchor2, scale, FUN, dots, chromsome){
    xinr <- (inRange(start(anchor1), scale) &
               as.character(seqnames(anchor1)) %in% chromsome) |
      (inRange(end(anchor2), scale) &
         as.character(seqnames(anchor2)) %in% chromsome)
    anchor1 <- anchor1[xinr]
    anchor2 <- anchor2[xinr]
    xa <- (end(anchor1) + start(anchor2))/2
    xb <- (start(anchor1) + start(anchor2))/2
    xc <- (start(anchor1) + end(anchor2))/2
    xd <- (end(anchor1) + end(anchor2))/2
    ya <- (xa-end(anchor1)+1)/ym
    yb <- (xb-start(anchor1)+1)/ym
    yc <- (xc-start(anchor1)+1)/ym
    yd <- (xd-end(anchor1)+1)/ym
    if(FUN.NAME=="grid.polygon"){
      merged <- merge_shape(xa, xb, xc, xd, ya, yb, yc, yd)
      for(i in seq_along(merged)){
        args <- c(list(x=merged[[i]]$x,
                       y=merged[[i]]$y,
                       default.units="native"),
                  dots)
        do.call(FUN, args)
      }
    }else{
      for(i in seq_along(anchor1)){
        args <- switch (FUN.NAME,
                        "grid.lines" = c(list(x=c(start(anchor1)[i], xc[i], end(anchor2)[i]), 
                                              y=c(0, yc[i], 0), 
                                              default.units="native"),
                                         dots),
                        "grid.text" = c(list(x=xa[i], 
                                             y=yb[i], 
                                             default.units="native"),
                                        dots),
                        c(list(x=c(xa[i], xb[i], xc[i], xd[i]), 
                               y=c(ya[i], yb[i], yc[i], yd[i]), 
                               default.units="native"),
                          dots)
        )
        do.call(FUN, args)
      }
    }
  }
  plot1Danno <- function(posx, scale, dots){
    xinr <- inRange(abs(posx), scale)
    xc <- posx[xinr]
    for(i in seq_along(xc)){
      args <- c(list(x=c(abs(xc[i]), ifelse(xc[i]>0, scale[2], scale[1])), 
                     y=c(0, ifelse(xc[i]>0,
                                   (scale[2]-xc[i]+1)/ym,
                                   (abs(xc[i])-scale[1]+1)/ym)), 
                     default.units="native"),
                dots)
      do.call(grid.lines, args)
    }
    
  }
  plot2Danno <- function(gr, scale, FUN, dots, chromsome){
    xc <- (start(gr) + end(gr))/2
    yc <- (xc-start(gr)+1)/ym
    xinr <- (inRange(start(gr), scale) | inRange(end(gr), scale)) &
      as.character(seqnames(gr)) %in% chromsome
    for(i in seq_along(gr)){
      if(xinr[i]){
        args <- switch (FUN.NAME,
                        "grid.lines" = c(list(
                          x=c(start(gr)[i], xc[i], end(gr)[i]), 
                          y=c(0, yc[i], 0), 
                          default.units="native"),
                          dots),
                        "grid.text" = c(
                          list(x=xc[i], 
                               y=yc[i], 
                               default.units="native"),
                          dots),
                        c(list(x=xc[i], 
                               y=yc[i], 
                               default.units="native"),
                          dots)
        )
        do.call(FUN, args)
      }
    }
  }
  addAnno <- function(obj, xscale, FUN, dots, chromosome=lastTrackViewer$chromosome){
    if(is(obj, "GInteractions")){
      plotInteractionAnno(first(obj), second(obj), xscale, FUN, dots, chromosome)
    }else{
      if(is(obj, "numeric")){
        plot1Danno(obj, xscale, dots)
      }else{
        if(is(obj, "GRanges")){
          plot2Danno(obj, xscale, FUN, dots, chromosome)
        }
      }
    }
  }
  pushViewport(vp)
  for(i in idx){
    currentVP <- viewport(y=mean(c(yHeightsLimits1[i], yHeightsLimits2[i])),
                          height = yHeightsLimits2[i]-yHeightsLimits1[i])
    pushViewport(currentVP)
    pushViewport(viewport(x=0, y=lastTrackViewer$yHeightBottom[i], 
                          height=1-lastTrackViewer$yHeightBottom[i]-lastTrackViewer$yHeightTop[i], 
                          width=1, 
                          clip="on",
                          just=c(0,0), 
                          xscale=xscale, 
                          yscale=sort(yscales[[i]])))
    if(.dat2target[i]){## two interaction heatmap, back to back
      ## top triangle
      if("top" %in% panel){
        pushViewport(viewport(x=0, y=.5, 
                              height=.5, 
                              width=1, 
                              clip="on",
                              default.units = "npc",
                              just=c(0,0), 
                              xscale=xscale, 
                              yscale=sort(yscales[[i]])))
        addAnno(obj, xscale, FUN, dots, chromosome=lastTrackViewer$chromosome)
        popViewport()
      }
      if("bottom" %in% panel){
        ## bottom triangle
        pushViewport(viewport(x=0, y=0, 
                              height=.5, 
                              width=1, 
                              clip="on", 
                              default.units = "npc",
                              just=c(0,0), 
                              xscale=xscale, 
                              yscale=rev(sort(yscales[[i]]))))
        addAnno(obj, xscale, FUN, dots, chromosome=lastTrackViewer$chromosome)
        popViewport()
      }
    }else{
      addAnno(obj, xscale, FUN, dots, chromosome=lastTrackViewer$chromosome)
    }
    
    popViewport()
    popViewport()
  }
  popViewport()
  return(invisible(vp))
}

#' @importFrom utils combn
merge_shape <- function(xa, xb, xc, xd, ya, yb, yc, yd){
  rects <- lapply(seq_along(xa), function(i){
    return(list(x=c(xa[i], xb[i], xc[i], xd[i]),
                y=c(ya[i], yb[i], yc[i], yd[i])))
  })
  r_tree <- list()
  others <- seq_along(rects)
  if(length(others)<2){
    return(rects)
  }
  cbn <- combn(others, 2, simplify = FALSE)
  ol <- vapply(cbn, FUN=function(.ele){
    overlap_region(rects[.ele])
  }, FUN.VALUE = logical(1L))
  cbn <- cbn[ol]
  cbn <- reduce_cbn(cbn)
  others <- setdiff(others, unlist(cbn))
  r_tree <- lapply(cbn, function(.ele){
    filter_points(rects[.ele])
  })
  ## merge points
  r_tree <- mapply(r_tree, cbn, FUN=function(keep, ids){
    pts <- rects[ids]
    x <- unlist(mapply(pts, keep, FUN=function(pt, idx){
      pt$x[idx]
    }, SIMPLIFY = FALSE))
    y <- unlist(mapply(pts, keep, FUN=function(pt, idx){
      pt$y[idx]
    }, SIMPLIFY = FALSE))
    pts <- list(x=x, y=y)
    pts <- do.call(sort_pts, pts)
  }, SIMPLIFY = FALSE)
  r_tree <- c(r_tree, rects[others])
  return(r_tree)
}
overlap_region <- function(xy){
  rx1 <- range(xy[[1]]$x)
  ry1 <- range(xy[[1]]$y)
  rx2 <- range(xy[[2]]$x)
  ry2 <- range(xy[[2]]$y)
  return((inRange(rx1[1], rx2)||inRange(rx1[2], rx2)) &&
           (inRange(ry1[1], ry2)||inRange(ry1[2], ry2)))
}
in_rect <- function(pt, polygon4){
  # calculate the sum of the areas of all possible triangles
  # if the sum is greater than the area of the polygon, the point is outside
  # else if the area of any of the triangles is 0, the point is in the edge
  # else the point is inside the polygon
  cbm <- list(c(1, 4), c(3, 4), c(2, 3), c(1, 2))
  areas <- vapply(cbm, function(.ele){
    xs <- c(polygon4$x[.ele], pt$x)
    ys <- c(polygon4$y[.ele], pt$y)
    get_area_triangle(xs, ys)
  }, FUN.VALUE = numeric(1L))
  area_rect <- get_area_4edge_polygon(polygon4)
  if(sum(areas)>area_rect){
    return(1) ## outside
  }
  if(any(areas==0)){
    return(0) ## at on edge
  }
  return(-1) ## inside
}
get_area_triangle <- function(xs, ys){
  abs(xs[1]*(ys[2]-ys[3])+xs[2]*(ys[3]-ys[1])+xs[3]*(ys[1]-ys[2]))/2
}
get_area_4edge_polygon <- function(polygon){
  #abcd, area(abc) + area(adc)
  abc_x <- polygon$x[c(1, 2, 3)]
  abc_y <- polygon$y[c(1, 2, 3)]
  adc_x <- polygon$x[c(1, 3, 4)]
  adc_y <- polygon$y[c(1, 3, 4)]
  get_area_triangle(abc_x, abc_y) + get_area_triangle(adc_x, adc_y)
}
filter_points <- function(rts){
  # if inside of a rect, remove the points
  lapply(seq_along(rts), function(id){
    .ele <- rts[[id]]
    pts <- lapply(seq_along(.ele$x), function(i){
      list(x=.ele$x[i], y=.ele$y[i])
    })
    keep <- vapply(pts, FUN = function(pt){
      vs <- vapply(rts[-id], FUN = function(pg4){
        in_rect(pt, pg4)
      }, FUN.VALUE = numeric(1L))
      all(vs>0)
    }, FUN.VALUE = logical(1L))
  })
}
reduce_cbn <- function(cbn){
  new_cbn <- list()
  for(i in seq_along(cbn)){
    id <- in_cbn(cbn[[i]], new_cbn)
    if(length(id)>0){
      if(length(id)==1){
        new_cbn[[id[1]]] <- sort(unique(c(new_cbn[[id[1]]], cbn[[i]])))
      }else{
        new_cbn[[id[1]]] <- sort(unique(c(unlist(new_cbn[id]), cbn[[i]])))
        new_cbn <- new_cbn[-id[-1]]
      }
    }else{
      new_cbn <- c(new_cbn, cbn[i])
    }
  }
  return(new_cbn)
}
in_cbn <- function(x, y){
  id <- vapply(y, FUN=function(.ele){
    any(x %in% .ele)
  }, FUN.VALUE = logical(1L))
  which(id)
}
sort_pts <- function(x, y){
  ## sort points by anti-clockwise
  x1 <- x-mean(x)
  y1 <- y-mean(y)
  ord <- order(atan2(y1, x1), sqrt(x1^2 + y1^2))
  list(x=x[ord], y=y[ord])
}

