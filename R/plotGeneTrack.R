plotGeneTrack <- function(track, xscale, chr, yaxis.gp=gpar(), lollipop_style_switch_limit=10){
  if(length(track@dat$featureID)!=length(track@dat)){
    return(plotGeneModel(track, xscale, chr, yaxis.gp, lollipop_style_switch_limit))
  }
  
  if(length(track@dat2)>0){
    feature.height <- 
      if(is.list(track@dat2$feature.height)) track@dat2$feature.height[[1]] else track@dat2$feature.height[1]
    if(length(feature.height)==0) feature.height <- .5
    unit <- feature.height/4
    y <- feature.height/2
  }else{
    unit <- 0.25
    y <- 0.5
  }
  
  ## split transcripts by featureID, 
  ## get the ranges of each transcripts
  ## assign lines for all transcripts
  ## plot center line for each featureID
  ## plot exons by the feature type
  ## plot direction at TSS
  ## add gene name at TSS
  vp <- viewport(xscale = xscale, y=y, height = unit*2)
  pushViewport(vp)
  trs <- split(track@dat, track@dat$featureID)
  col <- track@style@color
  rgs <- unlist(GRangesList(sapply(trs, function(.ele)range(unname(.ele)))))
  if(length(trs)!=length(rgs)){
    stop("some genes are splited into multiple regions, eg multple strand.")
  }
  strand <- as.character(strand(rgs))
  dj <- disjoin(rgs, with.revmap=TRUE)
  lineId <- data.frame(id=seq_along(rgs), line=0)
  revmap <- dj$revmap
  for(i in seq_along(revmap)){
    revm <- sort(revmap[[i]])
    if(lineId[revm[1], "line"]==0){
      lineId[revm[1], "line"] <- 1
    }
    for(j in seq_along(revm)[-1]){
      lineId[revm[j], "line"] <- lineId[revm[j-1], "line"] + 1
    }
  }
  totalLines <- max(lineId$line)
  ## check height of plot region
  unity <- convertHeight(unit(1, "npc"), unitTo = "lines", valueOnly = TRUE)
  doLabels <- unity>=totalLines
  eachLineHeight <- 1/totalLines
  currLineBottom <- 1
  unitx <- convertWidth(unit(1, "npc"), unitTo = "lines", valueOnly = TRUE)
  arr_size <- min(.5, unitx/length(trs)) #determine arrow size by width/#events
  if(!any(c("fontsize", "cex") %in% names(track@style@ylabgp))){
    if(doLabels){
      track@style@ylabgp$cex <- optFontSize1(.2*y)
    }else{
      track@style@ylabgp$cex <- arr_size
    }
  }
  for(i in seq.int(totalLines)){
    vp <- viewport(y = currLineBottom - eachLineHeight/2, 
                   height = eachLineHeight, xscale = xscale)
    pushViewport(vp)
    if(doLabels){
      gene_y <- .75
      gene_h <- .5
    }else{
      gene_y <- .5
      gene_h <- 1
    }
    if(length(track@dat2)>0){
      ## add a baseline for lollipop plot
      grid.lines(x=xscale, 
                 y=c(gene_y, gene_y), 
                 gp = gpar(col=track@style@color),
                 default.units = "native")
    }
    
    trs.sub <- trs[which(lineId$line==i)]
    rgs.sub <- rgs[which(lineId$line==i)]
    str.sub <- strand[which(lineId$line==i)]
    oid <- order(start(rgs.sub))
    trs.sub <- trs.sub[oid]
    rgs.sub <- rgs.sub[oid]
    str.sub <- str.sub[oid]
    stringStopPos <- c(0, 0)
    for(j in seq_along(trs.sub)){
      curr_trs <- trs.sub[[j]]
      curr_rg <- rgs.sub[j]
      curr_str <- str.sub[j]
      str_neg <- curr_str=="-"
      this_color <- if(is.list(curr_trs$color)) curr_trs$color[[1]] else if(length(curr_trs$color)) curr_trs$color else col
      this_height <- 
        if(is.list(curr_trs$height)) curr_trs$height[[1]] else if(length(curr_trs$height)) curr_trs$height else gene_h/ifelse(curr_trs$feature %in% c("CDS", "exon"), 2, 4)
      hide_this_label <- FALSE
      if(length(curr_trs$hide_label)>0){
        hide_this_label <- unlist(curr_trs$hide_label)[1]
        if(!is.logical(hide_this_label)){
          hide_this_label <- FALSE
        }
      }
      ## plot center line
      grid.lines(x=c(start(curr_rg), end(curr_rg)), 
                 y=c(gene_y, gene_y), 
                 gp = gpar(col=col),
                 default.units = "native")
      ## plot exons by the feature type
      grid.rect(x=(start(curr_trs)+end(curr_trs))/2, y=gene_y, 
                width =width(curr_trs), 
                height=this_height, 
                gp=gpar(col=NA, fill=this_color), 
                default.units = "native")
      ## plot direction at TSS
      stringW <- convertWidth(stringWidth(names(curr_rg)), unitTo = "native", valueOnly = TRUE)
      pushViewport(viewport(x=ifelse(!str_neg, start(curr_rg), end(curr_rg)),
                            y=gene_y, 
                            width = unit(max(min(gene_h/2, 2*width(curr_rg)/abs(diff(xscale))),
                                             convertWidth(unit(8, "lines"),
                                                          unitTo = "npc",
                                                          valueOnly = TRUE)), "snpc"), 
                            height = unit(gene_h/2, "snpc"),
                            just = c(ifelse(!str_neg, 0, 1), 1),
                            default.units = "native"))
      if(!hide_this_label){
        if(!str_neg){
          grid.lines(x=unit(c(0, 0, 1), "npc"),
                     y=unit(c(.5, 0, 0), "npc"),
                     arrow = arrow(type="closed", angle = 15, length = unit(arr_size, "lines")),
                     gp=gpar(col=this_color, fill=this_color))
          ## add gene name at TSS
          if(doLabels){
            if(start(curr_rg) >= stringStopPos[1]){
              grid.text(label = names(curr_rg), x = 0, y = 0, hjust = 0, vjust=1.5,
                        gp = do.call(gpar, track@style@ylabgp))
              stringStopPos[1] <- start(curr_rg)+stringW
            }else{
              if(start(curr_rg) >= stringStopPos[2]){
                grid.text(label = names(curr_rg), x = 0, y = 1, hjust = 0, vjust=-1,
                          gp = do.call(gpar, track@style@ylabgp))
                stringStopPos[2] <- start(curr_rg)+stringW
              }
            }
          }else{#dynamic label
            if(i %in% c(1, totalLines)){
              if(start(curr_rg) >= stringStopPos[1] &&
                 end(curr_rg) >= stringStopPos[1] + stringW){
                if(i==1){#top
                  grid.text(label = names(curr_rg), x = 0, y = 1,
                            hjust = 0, vjust=-1.5,
                            gp = do.call(gpar, track@style@ylabgp))
                }else{
                  grid.text(label = names(curr_rg), x = 0, y = 0,
                            hjust = 0, vjust=1.5,
                            gp = do.call(gpar, track@style@ylabgp))
                }
              }
              stringStopPos[1] <- end(curr_rg)
            }
          }
        }else{
          grid.lines(x=unit(c(1, 1, 0), "npc"),
                     y=unit(c(.5, 0, 0), "npc"),
                     arrow = arrow(type="closed", angle = 15, length = unit(arr_size, "lines")),
                     gp=gpar(col=this_color, fill=this_color))
          ## add gene name at TSS
          if(doLabels){
            if(end(curr_rg)-stringW >= stringStopPos[1]){
              grid.text(label = names(curr_rg), x = 1, y = 0, hjust = 1, vjust=1.5,
                        gp = do.call(gpar, track@style@ylabgp))
              stringStopPos[1] <- end(curr_rg)
            }else{
              if(end(curr_rg)-stringW >= stringStopPos[2]){
                grid.text(label = names(curr_rg), x = 1, y = 1, hjust = 1, vjust=-1,
                          gp = do.call(gpar, track@style@ylabgp))
                stringStopPos[2] <- end(curr_rg)
              }
            }
          }else{#dynamic label
            if(i %in% c(1, totalLines)){
              if(start(curr_rg) >= stringStopPos[1] &&
                 end(curr_rg) >= stringStopPos[1] + stringW){
                if(i==1){#bottom
                  grid.text(label = names(curr_rg), x = 1, y = 1,
                            hjust = 1, vjust=-1.5,
                            gp = do.call(gpar, track@style@ylabgp))
                }else{
                  grid.text(label = names(curr_rg), x = 1, y = 0,
                            hjust = 1, vjust=1.5,
                            gp = do.call(gpar, track@style@ylabgp))
                }
              }
              stringStopPos[1] <- end(curr_rg)
            }
          }
        }
      }
      popViewport()
      
        
    }
    popViewport()
    currLineBottom <- currLineBottom - eachLineHeight
  }
  popViewport()
  
  if(length(track@dat2)>0){
    track@dat2 <- orderedGR(track@dat2)
    xscale.gr <- GRanges(seqnames = chr, ranges=IRanges(min(xscale), max(xscale)))
    track@dat2 <- subsetByOverlaps(track@dat2, xscale.gr, ignore.strand=TRUE)
    if(length(track@dat2)<1) return(invisible())
    width(track@dat2) <- 1
    TYPES <- c("circle", "pie", "pin", "pie.stack", "flag")
    type <- if(is.list(track@dat2$type)) track@dat2$type[[1]] else track@dat2$type[1]
    if(length(type)==0) type <- "circle"
    if(!type %in% TYPES) type <- "circle"
    if(type=="pin"){ ## read the pin shape file
      pinpath <- system.file("extdata", "map-pin-red.xml", package="trackViewer")
      pin <- readPicture(pinpath)
    }else{
      pin <- NULL
    }
    cex <- if(is.list(track@dat2$cex)) track@dat2$cex[[1]] else track@dat2$cex[1]
    if(length(cex)==0) cex <- 1
    dashline.col <- if(is.list(track@dat2$dashline.col)) track@dat2$dashline.col[[1]] else track@dat2$dashline.col[1]
    if(length(dashline.col)==0) dashline.col <- "gray80"
    jitter <- if(is.list(track@dat2$jitter)) track@dat2$jitter[[1]] else track@dat2$jitter[1]
    if(length(jitter)==0) jitter <- "node"
    if(!jitter %in% c("node", "label")) jitter <- "node"
    scoreMax0 <- scoreMax <- 
      if(length(track@dat2$score)>0) ceiling(max(c(track@dat2$score, 1), na.rm=TRUE)) else 1
    if(type=="pie.stack") scoreMax <- length(unique(track@dat2$stack.factor))
    if(!type %in% c("pie", "pie.stack")){
      if(scoreMax>lollipop_style_switch_limit) {
        track@dat2$score <- 10*track@dat2$score/scoreMax
        scoreMax <- 10*scoreMax0/scoreMax
      }else{
        scoreMax <- scoreMax0
      }
      scoreType <- 
        if(length(track@dat2$score)>0) all(floor(track@dat2$score)==track@dat2$score) else FALSE
    }else{
      scoreType <- FALSE
    }
    LINEW <- as.numeric(convertX(unit(1, "line"), "npc"))/2
    LINEH <- as.numeric(convertY(unit(1, "line"), "npc"))/2
    ## GAP the gaps between any elements
    GAP <- .2 * LINEH
    ratio.yx <- getYXratio()
    plotLollipops(track@dat2, feature.height=y+unit, bottomHeight=0, baseline=y, 
                  type=type, ranges=xscale.gr, yaxis=FALSE, yaxis.gp=yaxis.gp,
                  scoreMax=scoreMax, scoreMax0=scoreMax0, scoreType=scoreType, 
                  LINEW, cex, ratio.yx, GAP, pin, dashline.col,
                  side="top", jitter=jitter)
  }
  return(invisible())
}
