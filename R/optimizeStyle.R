optFontSize <- function(axis, viewerStyle, height){
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
}

#' Optimize the style of plot
#' @description Automatic optimize the stlye of trackViewer
#' @param trackList An object of \code{\link{trackList}}
#' @param viewerStyle An object of \code{\link{trackViewerStyle}}
#' @param theme A character string. Could be "bw" or "col".
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
       !((is.list(trackList) && all(sapply(trackList, class)=="track")))){
        stop("trackList must be an object of \"trackList\"
             (See ?trackList) or a list of track")
    }
    if(trackList[[length(trackList)]]@type=="data" && viewerStyle@margin[3] < .02)
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
            if(trackList[[i]]@type=="data"){
                trackList[[i]]@style@yaxis@label <- TRUE
                trackList[[i]]@style@yaxis@gp <- 
                    c(trackList[[i]]@style@yaxis@gp, 
                      cex=optFontSize("y", viewerStyle))
                trackList[[i]]@style@marginTop <- .1
            }
          if(trackList[[i]]@type=="lollipopData"){
            trackList[[i]]@style@yaxis@draw <- FALSE
          }
        }
    }
    ##about ylab, direction, fontsize
    viewerStyle@autolas <- TRUE
    if(!is.null(theme)){
        if(theme=="bw"){
            for(i in 1:length(trackList)){
                if(trackList[[i]]@type %in% c("data", "lollipopData")){
                    trackList[[i]]@style@ylabpos="bottomleft"
                    trackList[[i]]@style@marginBottom=.2
                    trackList[[i]]@style@ylabgp=
                        list(cex=optFontSize1(.4*trackList[[i]]@style@height), 
                             col="black")
                }else{
                    trackList[[i]]@style@ylabpos="upstream"
                    trackList[[i]]@style@ylabgp=
                        list(cex=optFontSize1(trackList[[i]]@style@height), 
                             col="black")
                }
                trackList[[i]]@style@color=c("black", "black")
                trackList[[i]]@style@yaxis@main=FALSE
            }
            if(viewerStyle@margin[4] < .05) viewerStyle@margin[4] <- .05
        }
        if(theme=="col"){
            for(i in 1:length(trackList)){
                if(trackList[[i]]@type %in% c("data", "lollipopData")){
                    trackList[[i]]@style@ylabpos="bottomleft"
                    trackList[[i]]@style@marginBottom=.2
                    trackList[[i]]@style@ylabgp=
                        list(cex=optFontSize1(.4*trackList[[i]]@style@height), 
                             col=rep(palette(), ceiling(i/7))[i+1])
                }else{
                    trackList[[i]]@style@ylabpos="upstream"
                    trackList[[i]]@style@ylabgp=
                        list(cex=optFontSize1(trackList[[i]]@style@height), 
                             col=rep(palette(), ceiling(i/7))[i+1])
                }
                trackList[[i]]@style@color=
                    rep(rep(palette(), ceiling(i/7))[i+1], 2)
                trackList[[i]]@style@yaxis@main=FALSE
            }
            if(viewerStyle@margin[4] < .05) viewerStyle@margin[4] <- .05
        }
    }
    return(list(tracks=trackList, style=viewerStyle))
    }