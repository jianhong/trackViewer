optFontSize <- function(axis, viewerStyle, height){
  stopifnot(length(axis)==1)
    if(axis=="x"){
        if(viewerStyle@xlas %in% c(0, 1)){
            fontHeight <- convertHeight(stringHeight("0123456789"), 
                                        unitTo="npc", valueOnly=TRUE)
            cex <- viewerStyle@margin[1]/6/fontHeight
        }else{
            fontWidth <- convertWidth(stringWidth("012345678"), 
                                      unitTo="npc", valueOnly=TRUE)
            cex <- viewerStyle@margin[1]/6/fontWidth
        }
    }else{
        if(axis=="y"){
            fontWidth <- convertWidth(stringWidth("0123"), 
                                      unitTo="npc", valueOnly=TRUE)
            cex <- viewerStyle@margin[2]/1.5/fontWidth
        }else{
            fontHeight <- convertHeight(stringHeight("0123456789"), 
                                        unitTo="npc", valueOnly=TRUE)
            s <- min(c(height, viewerStyle@margin[2]))
            cex <- s/fontHeight
        }
    }
    
    if(cex > 1) cex <- 1
    if(cex < .3) cex <- .3
    return(cex)
}
optFontSize1 <- function(height){
    fontHeight <- convertHeight(stringHeight("0123456789"), 
                                unitTo="npc", valueOnly=TRUE)
    cex <- height/3/fontHeight
    if(cex>2) cex <- log2(cex)
    cex
}

#' Optimize the style of plot
#' @description Automatic optimize the stlye of trackViewer
#' @param trackList An object of \code{\link{trackList}}
#' @param viewerStyle An object of \code{\link{trackViewerStyle}}
#' @param theme A character string. Could be "bw", "col" or "safe".
#' @return a list of a \code{\link{trackList}} and a \code{\link{trackViewerStyle}}
#' @seealso See Also as \code{\link{viewTracks}}
#' @export
#' @examples 
#' extdata <- system.file("extdata", package="trackViewer",
#'                        mustWork=TRUE)
#' files <- dir(extdata, ".wig")
#' tracks <- lapply(paste(extdata, files, sep="/"), 
#'                  importScore, format="WIG")
#' re <- optimizeStyle(trackList(tracks))
#' trackList <- re$tracks
#' viewerStyle <- re$style

optimizeStyle <- function(trackList, viewerStyle=trackViewerStyle(), theme=NULL){
    if(missing(trackList))
        stop("trackList must be an object of 'trackList'")
    if(!is(trackList, "trackList") && 
       !(is.list(trackList) && 
         all(vapply(trackList, function(.ele) is(.ele, "track"), 
                    FUN.VALUE = TRUE)))){
        stop("trackList must be an object of \"trackList\"
             (See ?trackList) or a list of track")
    }
    if(trackList[[length(trackList)]]@type %in% c("data", "interactionData") && 
       viewerStyle@margin[3] < .02)
        viewerStyle@margin[3] <- .02
    ##put x-axis?
    dataTracksIdx <- sapply(trackList, function(.ele) .ele@type=="data")
    dataTracksXscale <- sapply(trackList[dataTracksIdx], 
                               function(.ele) .ele@style@xscale@draw)
    if(sum(dataTracksIdx)<5 && !(any(dataTracksXscale))){
        ##put x-axis
        viewerStyle@xaxis <- TRUE
        ##increase space for x-axis
        if(viewerStyle@margin[1] < .05) viewerStyle@margin[1] <- .05
        ##adjust x fontsize
        viewerStyle@xgp <- c(viewerStyle@xgp, 
                             cex=optFontSize("x", viewerStyle))
    }
    ##put x-scale?
    if(!viewerStyle@xaxis){
        if(all(!(sapply(trackList, function(.ele) .ele@style@xscale@draw)))){
            for(i in 1:length(trackList)){
                if(trackList[[i]]@type=="data"){
                    trackList[[i]]@style@xscale@draw <- TRUE
                    break
                }
            }
        }
    }
    ##put y-axis?
    if(all(!(sapply(trackList, function(.ele) .ele@style@yaxis@label)))){
        for(i in 1:length(trackList)){
            if(trackList[[i]]@type %in% c("data", "interactionData")){
                trackList[[i]]@style@yaxis@label <- TRUE
                trackList[[i]]@style@yaxis@gp <- 
                    c(trackList[[i]]@style@yaxis@gp, 
                      cex=optFontSize("y", viewerStyle))
                trackList[[i]]@style@marginTop <- .1
            }
          if(trackList[[i]]@type %in% c("lollipopData")){
            trackList[[i]]@style@yaxis@draw <- FALSE
          }
          if(trackList[[i]]@type %in% c("interactionData")){
            trackList[[i]]@style@yaxis@main <- FALSE
            viewerStyle@margin[4] <- .05
          }
        }
    }
    ##about ylab, direction, fontsize
    viewerStyle@autolas <- TRUE
    if(!is.null(theme)){
      defaultGeneYlabPos <- "upstream"
      colors <- switch(theme,
                       "bw"  = rep("black",length(trackList)+1),
                       "col" = rep(palette(), length(trackList)),
                       "safe"= rep(c("#000000", "#D55E00", "#009E73",
                                     "#0072B2", "#56B4E9", "#CC79A7",
                                     "#E69F00", "#BEBEBE"),
                                   length(trackList)),
                       rep("black",length(trackList)+1))
      for(i in seq_along(trackList)){
        if(trackList[[i]]@type %in% c("data", "lollipopData", "interactionData")){
          trackList[[i]]@style@ylabpos="bottomleft"
          trackList[[i]]@style@marginBottom=.2
          trackList[[i]]@style@ylabgp=
            list(cex=optFontSize1(.4*trackList[[i]]@style@height), 
                 col=ifelse(theme=="col", colors[i+1], "black"))
        }else{
          trackList[[i]]@style@ylabpos=defaultGeneYlabPos
          if(trackList[[i]]@type=="gene"){
            trackList[[i]]@style@ylabgp=
              list(cex=optFontSize1(.4*trackList[[i]]@style@height), 
                   col=colors[i+1])
          }else{
            trackList[[i]]@style@ylabgp=
              list(cex=optFontSize1(trackList[[i]]@style@height), 
                   col=colors[i+1])
          }
        }
        trackList[[i]]@style@color=
          rep(colors[i+1], 2)
        trackList[[i]]@style@yaxis@main=FALSE
      }
      if(viewerStyle@margin[4] < .05) viewerStyle@margin[4] <- .05
    }
    return(list(tracks=trackList, style=viewerStyle))
    }