plotInteractionDataTrack <- function(.dat, .dat2, scale, color, yscale, breaks, NAcolor="white"){
  names(.dat) <- NULL
  ## if there is target metadata
  .dat1target <- NULL
  .dat2target <- NULL
  if(length(.dat$target)==length(.dat) && is(.dat$target, "GRanges")){
    .dat1target <- .dat$target
    if(length(.dat2)!=0){
      if(length(.dat2$target)==length(.dat2) && is(.dat2$target, "GRanges")){
        names(.dat2) <- NULL
        .dat2target <- .dat2$target
        mcols(.dat2) <- mcols(.dat2)[, "score"]
        colnames(mcols(.dat2)) <- "score"
      }
    }
  }else{
    .dat1target <- .dat2
  }
  mcols(.dat) <- mcols(.dat)[, "score"]
  colnames(mcols(.dat)) <- "score"
  if(missing(yscale)) yscale <- c(0, 1)
  if(length(.dat)<1){
    return()
  }
  ## plot rect at position
  ## x = (center1 + center2)/2
  ## y = unit((x-scale[1]+1)/((scale[2] - scale[1] + 1)/2), "npc")
  ## width = width(x1)
  ## height = uint(width(x2)/(scale[2] - scale[1] + 1), "npc")
  ## rot = 45 degree
  ## color = colorRampPalette(color)(100)
  if(length(breaks)<3){
    if(length(.dat2target)>0){
      rg <- range(c(.dat$score[!is.na(.dat$score)], 
                    .dat2$score[!is.na(.dat2$score)]))
    }else{
      rg <- range(.dat$score[!is.na(.dat$score)])
    }
    if(length(rg)!=2){
      return()
    }
    rg <- rg + c(-1, 1)*diff(rg)/1000
    breaks <- seq(min(0, rg[1]), rg[2], length.out = 101)
  }
  
  if(length(unique(color))==1){
    if(!tolower(color[1]) %in% c("white", "#ffffff", "#fff")){
      color <- c("white", color[1])
    }else{
      color <- c("black", color[1])
    }
  }
  if(length(color)==0){
    color <- c("white", "red")
  }
  crp <- colorRampPalette(color)(length(breaks)-1)
  inRange <- function(x, scale){
    x>=scale[1] & x<=scale[2]
  }
  ym <- (scale[2]-scale[1] + 1)/2
  getMC <- function(scores){
    mc <- cut(scores, breaks = breaks, labels = crp)
    mc <- as.character(mc)
    mc[is.na(mc)] <- NAcolor
    mc
  }
  
  plotInteractionHeatmap <- function(anchor1, anchor2){
    xa <- (end(anchor1) + start(anchor2))/2
    xb <- (start(anchor1) + start(anchor2))/2
    xc <- (start(anchor1) + end(anchor2))/2
    xd <- (end(anchor1) + end(anchor2))/2
    ya <- (xa-end(anchor1)+1)/ym
    yb <- (xb-start(anchor1)+1)/ym
    yc <- (xc-start(anchor1)+1)/ym
    yd <- (xd-end(anchor1)+1)/ym
    irx <- inRange(xa, scale) | inRange(xb, scale) | inRange(xc, scale) | inRange(xd, scale)
    iry <- inRange(ya, yscale) | inRange(yb, yscale) | inRange(yc, yscale) | inRange(yd, yscale)
    mc <- getMC(anchor1$score)
    for(i in seq_along(anchor1)){
      if(irx[i] && iry[i]){
        grid.polygon(x=c(xa[i], xb[i], xc[i], xd[i]), 
                     y=c(ya[i], yb[i], yc[i], yd[i]), 
                     default.units="native",
                     gp = gpar(fill=mc[i], col = NA))
      }
    }
  }
  if(length(.dat2target)){## two interaction heatmap, back to back
    ## top triangle
    pushViewport(viewport(x=0, y=.5, 
                          height=.5, 
                          width=1, 
                          clip="on",
                          default.units = "npc",
                          just=c(0,0), 
                          xscale=scale, 
                          yscale=yscale))
    plotInteractionHeatmap(.dat, .dat1target)
    popViewport()
    ## bottom triangle
    pushViewport(viewport(x=0, y=0, 
                          height=.5, 
                          width=1, 
                          clip="on", 
                          default.units = "npc",
                          just=c(0,0), 
                          xscale=scale, 
                          yscale=rev(yscale)))
    plotInteractionHeatmap(.dat2, .dat2target)
    popViewport()
  }else{
    plotInteractionHeatmap(.dat, .dat1target)
  }
  # legend
  vp <- viewport(x = 1 - convertWidth(unit(5, "char"), "npc", valueOnly = TRUE), 
                 y= 1 - convertHeight(unit(1, "char"), "npc", valueOnly = TRUE), 
                 width = convertWidth(unit(5, "char"), "npc", valueOnly = TRUE),
                 height = convertHeight(unit(.5, "char"), "npc", valueOnly = TRUE),
                 default.units = "npc", xscale = range(breaks))
  pushViewport(vp)
  grid.raster(t(crp), width=1, height=1)
  grid.xaxis(at=round(range(breaks)), gp=gpar(fontsize=8))
  popViewport()
}

