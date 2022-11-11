plotScSeqTrack <- function(.dat, .dat2, scale, color, yscale, breaks, NAcolor="white"){
  if(length(names(.dat))<1){## no cell barcodes
    return()
  }
  if(length(.dat)<1){
    return()
  }
  ## remove unused mcols
  mcols(.dat) <- mcols(.dat)[, "score"] 
  colnames(mcols(.dat)) <- "score"
  ## split the data by cell barcode
  .dat <- split(.dat, names(.dat))
  ## resort the data by signals
  .dat <- .dat[order(vapply(.dat, FUN=function(.ele){
    mean(.ele$score, na.rm = TRUE)
  }, FUN.VALUE = 0.0), decreasing = TRUE)]
  id <- rep(seq_along(.dat), lengths(.dat))
  .dat <- unlist(GRangesList(.dat), use.names = FALSE)
  .dat$id <- id
  ## convert GRanges to matrix or plot it line by line?
  if(missing(yscale)) yscale <- c(0, 1)
  ## plot rect at position
  ## x = c(start, start, end, end)
  ## heightPerCell <- diff(yscale)/(numberOfCell+1)
  ## y = c(i-.5, i+.5, i+.5, i-.5)*heightPerCell # i is the cell order
  ## color = colorRampPalette(color)(100)
  if(length(breaks)<3){
    rg <- range(.dat$score[!is.na(.dat$score)])
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
  mc <- cut(.dat$score, breaks = breaks, labels = crp)
  mc <- as.character(mc)
  mc[is.na(mc)] <- NAcolor
  
  yh <- diff(yscale)/max(.dat$id)
  
  for(i in seq_along(.dat)){
      grid.polygon(x=c(start(.dat)[i], start(.dat)[i], 
                       end(.dat)[i], end(.dat)[i]), 
                   y=c(i-.5, i+.5, i+.5, i-.5)*yh, 
                   default.units="native",
                   gp = gpar(fill=mc[i], col = mc[i]))
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

