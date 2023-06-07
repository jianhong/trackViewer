#' plot the tracks
#' @description A function to plot the data for given range
#' @param trackList an object of \code{\link{trackList}}
#' @param chromosome chromosome
#' @param start start position
#' @param end end position
#' @param strand strand
#' @param gr an object of \code{\link[GenomicRanges:GRanges-class]{GRanges}}
#' @param ignore.strand ignore the strand or not when do filter. default TRUE
#' @param viewerStyle an object of \code{\link{trackViewerStyle}}
#' @param autoOptimizeStyle should use \code{\link{optimizeStyle}} to optimize style
#' @param newpage should be draw on a new page?
#' @param operator operator, could be +, -, *, /, ^, \%\%, and NA. "-" means dat - dat2, 
#' and so on. NA means do not apply any operator.
#' Note: if multiple operator is supplied, please make sure the length of operator keep
#' same as the length of trackList.
#' @param smooth logical(1) or numeric(). Plot smooth curve or not. If it is numeric, eg n,
#' mean of nearby n points will be used for plot.
#' If it is numeric, the second number will be the color. Default coloer is 2 (red).
#' @param lollipop_style_switch_limit The cutoff value for lollipop style for the 'circle' type.
#' If the max score is greater than this cutoff value, trackViewer will only plot one shape at
#' the highest score. Otherwise trackViewer will draw the shapes like `Tanghulu`.
#' @return An object of \code{\link[grid]{viewport}} for \code{\link{addGuideLine}}
#' @import GenomicRanges
#' @import grid
#' @importFrom grDevices as.raster col2rgb colorRampPalette palette
#' @importFrom scales rescale
#' @importFrom stats ksmooth
#' @importFrom IRanges viewMaxs Views
#' @export
#' @seealso See Also as \code{\link{addGuideLine}}, \code{\link{addArrowMark}}
#' @examples 
#' extdata <- system.file("extdata", package="trackViewer",
#'                        mustWork=TRUE)
#' files <- dir(extdata, "-.wig")
#' tracks <- lapply(paste(extdata, files, sep="/"), 
#'                  importScore, format="WIG")
#' tracks <- lapply(tracks, function(.ele) {strand(.ele@dat) <- "-"; .ele})
#' fox2 <- importScore(paste(extdata, "fox2.bed", sep="/"), format="BED")
#' dat <- coverageGR(fox2@dat)
#' fox2@dat <- dat[strand(dat)=="+"]
#' fox2@dat2 <- dat[strand(dat)=="-"]
#' gr <- GRanges("chr11", IRanges(122929275, 122930122), strand="-")
#' viewTracks(trackList(track=tracks, fox2=fox2), gr=gr, autoOptimizeStyle=TRUE)

viewTracks <- function(trackList, chromosome, start, end, strand, gr=GRanges(),
                       ignore.strand=TRUE,
                       viewerStyle=trackViewerStyle(), autoOptimizeStyle=FALSE,
                       newpage=TRUE, operator=NULL, smooth=FALSE,
                       lollipop_style_switch_limit=10){
  if(!is.null(operator)){
    if(!all(operator %in% c("+", "-", "*", "/", "^", "%%", NA))){
      stop('operator must be one of "+", "-", "*", "/", "^", "%%", NA')
    }
  }
  stopifnot(is.numeric(smooth)||is.logical(smooth))
  if(missing(trackList)){
    stop("trackList is required.")
  }
  if(!is(trackList, "trackList") && 
     !(is.list(trackList) && 
        all(vapply(trackList, function(.ele) is(.ele, "track"), 
                   FUN.VALUE = TRUE)))){
    stop("trackList must be an object of \"trackList\"
         (See ?trackList) or a list of track")
  }
  filterTracksFlag <- TRUE
  if(missing(gr)){
    if(missing(chromosome) || missing(start) || missing(end))
      stop("Please input the coordinate.")
    if(missing(strand) || !strand %in% c("+", "-"))
      strand <- "*"
    wavyLine <- FALSE
  }else{
    if(!is(gr, "GRanges")){
      stop("gr must be an object of GRanges.")
    }
    if(length(gr)==0){
      stop("the length of gr must greater than 0.")
    }
    if(length(gr$filterTracks)>0){
      filterTracksFlag <- as.logical(gr$filterTracks[1])
    }
    grs <- range(gr)
    chromosome <- as.character(GenomicRanges::seqnames(grs))[1]
    start <- GenomicRanges::start(grs)[1]
    end <- GenomicRanges::end(grs)[1]
    strand <- as.character(GenomicRanges::strand(grs))[1]
    if(ignore.strand) strand <- "*"
    if(length(gr)>1){ ## more than one region
      if(length(unique(as.character(seqnames(gr))))>1){
        stop("gr has multiple seqnames.")
      }
      if(viewerStyle@flip){
        gr <- sort(gr, decreasing = TRUE)
      }else{
        gr <- sort(gr, decreasing = FALSE)
      }
      if(length(gr$percentage)==0){
        gr$percentage <- 1/length(gr)
      }else{
        gr$percentage <- gr$percentage/sum(gr$percentage, na.rm = TRUE)
        gr$percentage[is.na(gr$percentage)] <- 0
      }
      gr$wavyLine <- TRUE
      gr$wavyLine[length(gr)] <- FALSE
      gr$filterTracks <- FALSE
      if(newpage) grid.newpage()
      
      if(autoOptimizeStyle){
        opt <- optimizeStyle(trackList, viewerStyle)
        trackList <- opt$tracks
        viewerStyle <- opt$style
      }
      trackList <- filterTracks(trackList, chromosome, start, end, strand) ## there is a bug if gr has different strand.
      for(i in seq_along(gr)){
        current.viewerStyle <- viewerStyle
        current.trackList <- trackList
        hgap <- convertWidth(unit(0.25, "lines"), 
                             unitTo = "npc", 
                             valueOnly = TRUE)
        if(i==1){
          setTrackViewerStyleParam(current.viewerStyle, 
                                   "margin", 
                                   c(current.viewerStyle@margin[1],
                                     current.viewerStyle@margin[2]/gr$percentage[i],
                                     current.viewerStyle@margin[3],
                                     hgap/gr$percentage[i]))
          for(j in seq_along(current.trackList)){
            tLjStyle <- current.trackList[[j]]@style
            if(tLjStyle@yaxis@draw && !tLjStyle@yaxis@main){
              setTrackYaxisParam(current.trackList[[j]],
                                 "draw",
                                 FALSE)
            }
            if(!tLjStyle@ylabpos %in% 
               c("upstream", "downstream") &&
               grepl("right", tLjStyle@ylabpos)){
              names(current.trackList)[[j]] <- ""
            }else{
              tLjStyle@yaxis@gp$cex <- tLjStyle@yaxis@gp$cex / gr$percentage[i]
              setTrackYaxisParam(current.trackList[[j]],
                                 "gp", tLjStyle@yaxis@gp)
            }
          }
        }else{
          if(i==length(gr)){
            setTrackViewerStyleParam(
              current.viewerStyle, 
              "margin", 
              c(current.viewerStyle@margin[1],
                hgap/gr$percentage[i],
                current.viewerStyle@margin[3],
                current.viewerStyle@margin[4]/gr$percentage[i])
            )
            for(j in seq_along(current.trackList)){
              tLjStyle <- current.trackList[[j]]@style
              if(tLjStyle@yaxis@draw && tLjStyle@yaxis@main){
                setTrackYaxisParam(current.trackList[[j]],
                                   "draw",
                                   FALSE)
              }
              if(!tLjStyle@ylabpos %in% 
                 c("upstream", "downstream") &&
                 grepl("left", tLjStyle@ylabpos)){
                names(current.trackList)[[j]] <- ""
              }else{
                tLjStyle@yaxis@gp$cex <- tLjStyle@yaxis@gp$cex / gr$percentage[i]
                setTrackYaxisParam(current.trackList[[j]],
                                   "gp", tLjStyle@yaxis@gp)
              }
            }
          }else{
            setTrackViewerStyleParam(
              current.viewerStyle, 
              "margin", 
              c(current.viewerStyle@margin[1],
                hgap/gr$percentage[i],
                current.viewerStyle@margin[3],
                hgap/gr$percentage[i])
            )
            for(j in seq_along(current.trackList)){
              setTrackYaxisParam(current.trackList[[j]],
                                 "draw",
                                 FALSE)
              if(!current.trackList[[j]]@style@ylabpos %in% 
                 c("upstream", "downstream")){
                names(current.trackList)[[j]] <- ""
              }else{
                tLjStyle@yaxis@gp$cex <- tLjStyle@yaxis@gp$cex / gr$percentage[i]
                setTrackYaxisParam(current.trackList[[j]],
                                   "gp", tLjStyle@yaxis@gp)
              }
            }
          }
        }
        pushViewport(vp=viewport(
          x = unit(ifelse(i==1, 
                          0, 
                          sum(gr$percentage[seq_len(i-1)]))+
                     gr$percentage[i]/2-hgap/2, 
                   "npc"),
          width = unit(gr$percentage[i]-hgap, "npc"),
          name = paste0("vp.track.v.", i)
        ))
        vp <- viewTracks(trackList = current.trackList,
                         gr = gr[i],
                         ignore.strand = ignore.strand,
                         viewerStyle = current.viewerStyle,
                         autoOptimizeStyle = FALSE,
                         newpage = FALSE, operator = operator)
        vp$name <- paste0("panel.", i)
        pushViewport(vp)
        upViewport()
        upViewport()
      }
      return(invisible(current.vpTree()))
    }
    wavyLine <- ifelse(length(gr$wavyLine)>0, as.logical(gr$wavyLine[1]),
                       FALSE)
  }
  if(ignore.strand) strand <- "*"
  if(end < start) stop("end should be greater than start.")
  
  if(!is(viewerStyle, "trackViewerStyle")){
    stop("viewerStyle must be an object of 'trackViewerStyle'.")
  }
  if(!is.logical(newpage)) stop("newpage should be a logical vector.")
  
  if(filterTracksFlag) {
    trackList <- filterTracks(trackList, chromosome, start, end, strand)
  }
  
  if(!is.null(operator)){
    ##change dat as operator(dat, dat2)
    ##dat2 no change
    operator <- rep(operator, length(trackList))[seq_along(trackList)]
    trackList <- mapply(trackList, operator, FUN=function(.ele, .operator){
      if(.ele@type=="data" && 
         length(.ele@dat2)>0 &&
         length(.ele@dat)>0 &&
         !is.na(.operator)){
        .ele@dat <- GRoperator(.ele@dat, 
                               .ele@dat2, 
                               col="score", 
                               operator=.operator)
        if(.operator !="+") .ele@dat2 <- GRanges()
      }else{
        warning("operator not supported for the inputs.")
      }
      .ele
    }, SIMPLIFY = FALSE)
  }else{
    operator <- rep(NA, length(trackList))
  }
  
  if(newpage) grid.newpage()
  
  if(autoOptimizeStyle){
    opt <- optimizeStyle(trackList, viewerStyle)
    trackList <- opt$tracks
    viewerStyle <- opt$style
  }
  
  margin <- viewerStyle@margin
  
  ##xscale, yscale
  xscale <- c(start, end)
  
  yscales <- getYlim(trackList, operator)
  yHeights <- getYheight(trackList)
  
  pushViewport(viewport(x=margin[2], y=margin[1], width=1-margin[2]-margin[4],
                        height=1-margin[1]-margin[3], just=c(0,0)))
  if(viewerStyle@xaxis){##draw x axis in bottom
    drawXaxis(xscale, viewerStyle)
  }
  popViewport()
  
  pushViewport(viewport(x=0, y=margin[1], width=1,
                        height=1-margin[1]-margin[3], just=c(0,0)))
  total <- length(trackList)
  #    if(interactive()) pb <- txtProgressBar(min=0, max=total, style=3)
  ht <- 0
  xy <- vector(mode="list", length=total)
  yHeightBottom <- yHeightTop <- rep(0.01, total)
  for(i in 1:total){
    xy[[i]] <- plotTrack(names(trackList)[i], trackList[[i]], 
                         viewerStyle, ht,
                         yscales[[i]], yHeights[i], xscale,
                         chromosome, strand, operator[i], wavyLine,
                         smooth=smooth,
                         lollipop_style_switch_limit=lollipop_style_switch_limit)
    ht <- ht + yHeights[i]
    if(length(trackList[[i]]@style@marginBottom)>0){
      yHeightBottom[i] <- trackList[[i]]@style@marginBottom
    }
    if(length(trackList[[i]]@style@marginTop)>0){
      yHeightTop[i] <- trackList[[i]]@style@marginTop
    }
    #        if(interactive()) setTxtProgressBar(pb, i)
  }
  #    if(interactive()) close(pb)
  popViewport()
  if(viewerStyle@flip) xscale <- rev(xscale)
  options(LastTrackViewer=list(vp=viewport(x=margin[2], y=margin[1], 
                                           height=1 - margin[1]- margin[3], 
                                           width=1 -margin[2] - margin[4],
                                           just=c(0,0), xscale=xscale, yscale=c(0,1)),
                               xy=xy,
                               yHeights=yHeights,
                               yscales=yscales,
                               xscale=xscale,
                               yHeightBottom=yHeightBottom,
                               yHeightTop=yHeightTop,
                               viewerStyle=viewerStyle))
  return(invisible(viewport(x=margin[2], y=margin[1], 
                            height=1 - margin[1]- margin[3], 
                            width=1 -margin[2] - margin[4],
                            just=c(0,0), xscale=xscale, yscale=c(0,1))))
  }