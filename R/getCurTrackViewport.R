#' Get current track viewport
#' @description Get current track viewport for addGuideLine
#' @param curViewerStyle an object of \code{\link{trackViewerStyle}}
#' @param start start position of current track
#' @param end end position of current track
#' @return an object of \code{\link[grid]{viewport}}
#' @import grid
#' @export
#' @seealso See Also as \code{\link{addGuideLine}}
#' @examples 
#' vp <- getCurTrackViewport(trackViewerStyle(), 10000, 10200)
#' addGuideLine(c(10010, 10025, 10150), vp=vp)

getCurTrackViewport <- function(curViewerStyle, start, end){
  if(!is(curViewerStyle, "trackViewerStyle"))
    stop("curViewerStyle must be an object of trackViewerStyle")
  margin <- curViewerStyle@margin
  xscale <- c(start, end)
  if(curViewerStyle@flip) xscale <- rev(xscale)
  return(viewport(x=margin[2], y=margin[1], 
                  height=1 - margin[1]- margin[3], 
                  width=1 -margin[2] - margin[4],
                  just=c(0,0), 
                  xscale=xscale, 
                  yscale=c(0,1)))
}
