#' Add guide lines to the tracks
#' @description A function to add lines for emphasizing the positions 
#' @param guideLine The genomic coordinates to draw the lines
#' @param col A vector for the line color
#' @param lty A vector for the line type
#' @param lwd A vector for the line width
#' @param vp A Grid viewport object. It must be output of \code{\link{viewTracks}}
#' @import grid
#' @export
#' @return NULL
#' @seealso See Also as \code{\link{getCurTrackViewport}}, \code{\link{addArrowMark}}, 
#' \code{\link{viewTracks}}
#' @examples
#' vp <- getCurTrackViewport(trackViewerStyle(), 10000, 10200)
#' addGuideLine(c(10010, 10025, 10150), vp=vp)

addGuideLine <- function(guideLine, col="gray", lty="dashed", lwd=1, vp=NULL){
  if(missing(guideLine) | 
     !(class(guideLine) %in% c("numeric", "integer")) |
     length(guideLine) < 1)
    stop("guideLine is required as a numeric vector of coordinates of genome")
  len <- length(guideLine)
  trimLen <- function(obj, len){
    if(length(obj)<len) 
      obj <- rep(obj, len)[1:len]
    obj
  }
  vpmultiple <- FALSE
  if(length(vp)>0){
    stopifnot(is(vp, "viewport"))
    if(is(vp, "vpTree")){
      ## check the range
      vpmultiple <- TRUE
    }
  }
  selectVP <- function(x, tree){
    xscales <- sapply(tree$children, function(.ele){
      seekViewport(names(.ele$children))
      current.viewport()$xscale
    }, simplify = FALSE)
    xscales <- do.call(cbind, xscales)
    i <- which(x>=xscales[1, ] & x<=xscales[2, ])
    if(length(i)<1) {
      message(x, " out of the range.")
      return(NULL)
    }
    seekViewport(paste0("panel.", i))
    current.viewport()
  }
  col <- trimLen(col, len)
  lty <- trimLen(lty, len)
  lwd <- trimLen(lwd, len)
  for(i in seq_along(guideLine)){
    if(vpmultiple){
      currentVP <- selectVP(guideLine[i], vp)
      grid.lines(x=guideLine[i], y=c(0, 1), 
                 gp=gpar(col=col[i], lty=lty[i], lwd=lwd[i]), 
                 default.units="native")
    }else{
      currentVP <- vp
      grid.lines(x=guideLine[i], y=c(0, 1), 
                 gp=gpar(col=col[i], lty=lty[i], lwd=lwd[i]), 
                 default.units="native", vp=currentVP)
    }
  }
  return(invisible())
}
