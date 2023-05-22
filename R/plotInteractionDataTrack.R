plotInteractionDataTrack <- function(.dat, .dat2, scale, color, yscale, breaks,
                                     NAcolor="white", style="heatmap"){
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
        #mcols(.dat2) <- mcols(.dat2)[, "score"]
        #colnames(mcols(.dat2)) <- "score"
      }
    }
  }else{
    .dat1target <- .dat2
  }
  #mcols(.dat) <- mcols(.dat)[, "score"]
  #colnames(mcols(.dat)) <- "score"
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
    if(rg[1]==rg[2] && rg[1]==0){
      rg[2] <- 1
    }
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
    cols <- rep(NA, length(anchor1))
    if(length(anchor1$border_color)==length(anchor1) && length(anchor1)>0){
      cols <- anchor1$border_color
    }
    for(i in seq_along(anchor1)){
      if(irx[i] && iry[i]){
        grid.polygon(x=c(xa[i], xb[i], xc[i], xd[i]), 
                     y=c(ya[i], yb[i], yc[i], yd[i]), 
                     default.units="native",
                     gp = gpar(fill=mc[i], col = cols[i]))
      }
    }
  }
  
  getBezierPoints <- function(x1, x2, y1, y2, h){
    convX <- function(x){
      (x-scale[1])/diff(scale)
    }
    reveX <- function(x){
      x*diff(scale)+scale[1]
    }
    convY <- function(y){
      h*(y-min(y, na.rm=TRUE))/diff(range(y, na.rm=TRUE))
    }
    x1 <- convX(x1)
    x2 <- convX(x2)
    y1 <- convX(y1)
    y2 <- convX(y2)
    bg1 <- bezierGrob(c(x2, x2, y1, y1),
                      c(0, h, h, 0),
                      default.units = "npc")
    trace1 <- bezierPoints(bg1)
    bg2 <- bezierGrob(c(x1, x1, y2, y2),
                      c(0, h, h, 0),
                      default.units = "npc")
    trace2 <- bezierPoints(bg2)
    return(list(x = reveX(c(convertX(trace2$x[1], unitTo = "npc", valueOnly = TRUE),
                            convertX(trace1$x, unitTo = "npc", valueOnly = TRUE),
                            convertX(rev(trace2$x), unitTo = "npc", valueOnly = TRUE))),
                y = convY(c(trace2$y[1], trace1$y, rev(trace2$y)))))
  }
  grid.link <- function(x1, x2, y1, y2, h, ...){
    pts <- getBezierPoints(x1, x2, y1, y2, h)
    grid.polygon(x=pts$x,
                 y=pts$y,
                 ...)
  }
  plotInteractioLink <- function(anchor1, anchor2){
    mc <- getMC(anchor1$score)
    hs <- sqrt(anchor1$score/
                max(anchor1$score[!is.infinite(anchor1$score)], 
                    na.rm = TRUE))
    irx <- inRange(start(anchor1), scale) | inRange(end(anchor1), scale) |
      inRange(start(anchor2), scale) | inRange(end(anchor2), scale)
    for(i in seq_along(anchor1)){
      if(irx[i]){
        grid.link(x1=start(anchor1)[i],
                  x2=end(anchor1)[i],
                  y1=start(anchor2)[i],
                  y2=end(anchor2)[i],
                  h = hs[i],
                  default.units="native",
                  gp = gpar(fill=mc[i], col = NA, alpha=hs[i]))
      }
    }
  }
  FUN = switch (tolower(style),
    'heatmap' = plotInteractionHeatmap,
    'link' = plotInteractioLink,
    plotInteractionHeatmap
  )
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
    FUN(.dat, .dat1target)
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
    FUN(.dat2, .dat2target)
    popViewport()
  }else{
    FUN(.dat, .dat1target)
  }
  # legend in y axis
  return(list(crp=crp, breaks=breaks))
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

