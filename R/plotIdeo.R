#' plot ideogram
#' @description plot ideogram for one chromosome 
#' @param ideo output of \link{loadIdeogram}.
#' @param chrom A length 1 character vector of chromosome name.
#' @param gp parameters used for \link[grid]{grid.roundrect}.
#' @param ... parameters not used.
#' @import grid
#' @examples 
#' \dontrun{
#' ideo <- loadIdeogram("hg38")
#' library(grid)
#' grid.newpage()
#' plotIdeo(ideo)
#' }
#' 
#' 
plotIdeo <- function(ideo, chrom=seqlevels(ideo)[1], gp=gpar(fill=NA), ...){
  stopifnot(class(ideo)=="GRanges")
  stopifnot(length(chrom)==1)
  ideo <- ideo[seqnames(ideo) %in% chrom]
  strand(ideo) <- "*"
  ideo <- sort(ideo)
  ol <- findOverlaps(ideo, drop.self=TRUE, drop.redundant=TRUE, minoverlap = 2)
  if(length(ol)>0){
    stop("There is overlaps in ideo.")
  }
  ### gieStain #############################
  # #FFFFFF - gneg    - Giemsa negative bands
  # #999999 - gpos25  - Giemsa positive bands
  # #666666 - gpos50  - Giemsa positive bands
  # #333333 - gpos75  - Giemsa positive bands
  # #000000 - gpos100 - Giemsa positive bands
  # #660033 - acen    - centromeric regions
  # #660099 - gvar    - variable length heterochromatic regions
  # #6600cc - stalk   - tightly constricted regions on the short arms of
  #                     the acrocentric chromosomes
  colorSheme <- c(
    "gneg"    = "#FFFFFF",
    "acen"    = "#660033",
    "gvar"    = "#660099",
    "stalk"   = "#6600CC"
  )
  gposCols <- sapply(1:100, function(i){
    i <- as.hexmode(round(256-i*2.56, digits = 0))
    i <- toupper(as.character(i))
    if(nchar(i)==1) i <- paste0("0", i)
    return(paste0("#", i, i, i))
  })
  names(gposCols) <- paste0("gpos", 1:100)
  colorSheme <- c(gposCols, colorSheme)
  ### split the ideogram into p arm and q arm
  stopifnot(length(ideo$name)==length(ideo))
  ideoName <- tolower(as.character(ideo$name))
  ### first character should be p or q
  arms <- substr(ideoName, 1, 1)
  if(!all(arms %in% c("p", "q"))){
    stop("p arm or q arm only.")
  }
  ideo <- split(ideo, arms)
  ra <- lapply(ideo, function(.ele){## not use range to avoid warning
    c(min(start(.ele)) , max(end(.ele)))
  })
  if(all(c("p", 'q') %in% names(ra))){
    ## coordinates of p arm should be smaller than q arm
    if(ra$p[2]>ra$q[1]){
      stop("coordinates of p arm should be smaller than q arm")
    }
  }
  ## viewports with xscale=c(start.p.arm, end.q.arm)
  vp_pq <- viewport(xscale=range(unlist(ra)))
  pushViewport(vp_pq)
  GRanges2raster <- function(ranges, gieStain){
    wid <- width(ranges) - 1
    hcf <- function(x){
      if(length(x)>2){
        x <- x[-length(x)]
      }
      far <- range(nchar(as.character(x)))
      if(far[2]<3) return(1)
      fa <- 10^(far[2]:far[1]-2)
      fa <- fa[fa>1]
      for(i in fa){
        if(all(x %% i == 0)) return(i)
      }
      if(far[1]>2){
        return(10^(far[1]-2))
      }else{
        return(1)
      }
    }
    wid <- round(wid/hcf(wid), digits = 0)
    matrix(rep(colorSheme[as.character(gieStain)], wid), nrow=1)
  }
  
  rrpoints <- function(x) {
    left <- 0
    bottom <- 0
    right <- convertX(unit(1, "npc"), "inches", valueOnly=TRUE)
    top <- convertY(unit(1, "npc"), "inches", valueOnly=TRUE)
    rw <- convertWidth(x$r, "inches", valueOnly=TRUE)
    rh <- convertHeight(x$r, "inches", valueOnly=TRUE)
    r <- min(rw, rh)
    if(right>=top){
      rect <- c(x1=left+r, y1=bottom, x2=right-r, y2=top, r=r, horiz=1)
    }else{
      rect <- c(x1=left, y1=bottom+r, x2=right, y2=top-r, r=r, horiz=0)
    }
    return(rect)
  }
  
  inMask <- function(x, y, cx, cy, r){
    d <- sqrt((x-cx)^2+(y-cy)^2)
    return(d<=r)
  }
  imgpoints <- function(img, mask){
    img <- matrix(rep(img, 1000), nrow=1000, byrow = TRUE)
    x <- unit((1:ncol(img))/ncol(img), "npc")
    y <- unit((1000:1)/1000, "npc")
    xx <- convertX(x, "inch", valueOnly = TRUE)
    yy <- convertY(y, "inch", valueOnly = TRUE)
    if(mask["horiz"]==1){
      ## xx < mask["x1"] , leftCir
      for(i in which(xx < mask["x1"])){
        img[!inMask(xx[i], yy, cx=mask["x1"], cy=mask["r"], r=mask["r"]), i] <- NA
      }
      ## xx > mask["x2"], rightCir
      for(i in which(xx > mask["x2"])){
        img[!inMask(xx[i], yy, cx=mask["x2"], cy=mask["r"], r=mask["r"]), i] <- NA
      }
    }else{
      ## yy < mask["y1"], bottomCir
      for(i in which(yy < mask["y1"])){
        img[i, !inMask(xx, yy[i], cx=mask["r"], cy=mask["y1"], r=mask["r"])] <- NA
      }
      ## yy > mask["y2"], topCir
      for(i in which(yy > mask["y2"])){
        img[i, !inMask(xx, yy[i], cx=mask["r"], cy=mask["y2"], r=mask["r"])] <- NA
      }
    }
    
    img
  }
  # p arm
  if(length(ra$p)){
    vp_p <- viewport(x = ra$p[1], width=diff(ra$p), just="left", 
                     default.units="native", xscale=ra$p)
    pushViewport(vp_p)
    image <- GRanges2raster(ranges(ideo$p), ideo$p$gieStain)
    roundrect <- roundrectGrob(r=unit(.5, "snpc"))
    roundrect <- rrpoints(roundrect)
    image <- imgpoints(image, roundrect)
    grid.raster(image=image, height=1, width=1)
    grid.roundrect(r=unit(.5, "snpc"), gp=gp)
    popViewport()
  }
  # q arm
  if(length(ra$q)){
    vp_q <- viewport(x = ra$q[1], width=diff(ra$q), just="left", 
                     default.units="native", xscale=ra$q)
    pushViewport(vp_q)
    image <- GRanges2raster(ranges(ideo$q), ideo$q$gieStain)
    roundrect <- roundrectGrob(r=unit(.5, "snpc"))
    roundrect <- rrpoints(roundrect)
    image <- imgpoints(image, roundrect)
    grid.raster(image=image, height=1, width=1)
    grid.roundrect(r=unit(.5, "snpc"), gp=gp)
    popViewport()
  }
  popViewport()
}