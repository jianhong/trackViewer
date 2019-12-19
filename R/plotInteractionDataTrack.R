plotInteractionDataTrack <- function(.dat, .dat2, scale, color, yscale, breaks, NAcolor="white"){
  names(.dat) <- NULL
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
    rg <- range(.dat$score[!is.na(.dat$score)])
    if(length(rg)!=2){
      return()
    }
    breaks <- seq(0, ceiling(rg[2]), length.out = 101)
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
  mc <- cut(.dat$score, breaks = breaks, labels = crp)
  mc <- as.character(mc)
  mc[is.na(mc)] <- NAcolor
  inRange <- function(x, scale){
    x>=scale[1] & x<=scale[2]
  }
  ym <- (scale[2]-scale[1] + 1)/2
  xa <- (end(.dat) + start(.dat2))/2
  xb <- (start(.dat) + start(.dat2))/2
  xc <- (start(.dat) + end(.dat2))/2
  xd <- (end(.dat) + end(.dat2))/2
  ya <- (xa-end(.dat)+1)/ym
  yb <- (xb-start(.dat)+1)/ym
  yc <- (xc-start(.dat)+1)/ym
  yd <- (xd-end(.dat)+1)/ym
  irx <- inRange(xa, scale) | inRange(xb, scale) | inRange(xc, scale) | inRange(xd, scale)
  iry <- inRange(ya, yscale) | inRange(yb, yscale) | inRange(yc, yscale) | inRange(yd, yscale)
  for(i in seq_along(.dat)){
    if(irx[i] && iry[i]){
      grid.polygon(x=c(xa[i], xb[i], xc[i], xd[i]), 
                   y=c(ya[i], yb[i], yc[i], yd[i]), 
                   default.units="native",
                   gp = gpar(fill=mc[i], col = NA))
    }
  }
  # legend
  vp <- viewport(x = 1 - convertWidth(unit(5, "char"), "npc", valueOnly = TRUE), 
                 y= 1 - convertHeight(unit(1, "char"), "npc", valueOnly = TRUE), 
                 width = convertWidth(unit(5, "char"), "npc", valueOnly = TRUE),
                 height = convertHeight(unit(.5, "char"), "npc", valueOnly = TRUE),
                 default.units = "npc", xscale = range(breaks))
  pushViewport(vp)
  grid.raster(t(crp), width=1, height=1)
  grid.xaxis(at=range(breaks), gp=gpar(fontsize=8))
  popViewport()
}

