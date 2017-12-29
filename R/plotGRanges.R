#' plot GRanges data
#' @description A function to plot GRanges data for given range
#' @param \dots one or more objects of \code{\link[GenomicRanges]{GRanges}}
#' @param range an object of \code{\link[GenomicRanges]{GRanges}}
#' @param viewerStyle an object of \code{\link{trackViewerStyle}}
#' @param autoOptimizeStyle should use \code{\link{optimizeStyle}} to optimize style
#' @param newpage should be draw on a new page?
#' @return An object of \code{\link[grid]{viewport}} for \code{\link{addGuideLine}}
#' @import GenomicRanges
#' @export
#' @seealso See Also as \code{\link{addGuideLine}}, \code{\link{addArrowMark}}
#' @examples 
#' gr1 <- GRanges("chr1", IRanges(1:50, 51:100))
#' gr2 <- GRanges("chr1", IRanges(seq(from=10, to=80, by=5),
#'                                seq(from=20, to=90, by=5)))
#' vp <- plotGRanges(gr1, gr2, range=GRanges("chr1", IRanges(1, 100)))
#' addGuideLine(guideLine=c(5, 10, 50, 90), col=2:5, vp=vp)
#' 
#' gr <- GRanges("chr1", IRanges(c(1, 11, 21, 31), width=9), 
#'               score=c(5, 10, 5, 1))
#' plotGRanges(gr, range=GRanges("chr1", IRanges(1, 50)))

plotGRanges <- function(..., range=GRanges(), 
                        viewerStyle=trackViewerStyle(),
                        autoOptimizeStyle=FALSE,
                        newpage=TRUE){
    GRList <- list(...)
    n <- length(GRList)
    if(n==1 && inherits(GRList[[1]], c("list", "GRangesList"))){
        GRList <- GRList[[1]]
        names <- names(GRList)
        if(is.null(names)) names <- paste("track", 1:length(GRList), sep="")
    }else{
        dots <- substitute(list(...))[-1]
        names <- make.names(unlist(sapply(dots, deparse)))
    }
    sta <- sapply(GRList, inherits, what="GRanges")
    if(any(!sta)) stop("All inputs must be objests of GRanges")
    if(missing(range)){
        stop("range is required!")
    }
    if(length(range)==0){
        stop("range is not given!")
    }
    cov <- lapply(GRList, coverageGR)
    tracks <- lapply(cov, function(.ele){
        new("track", dat=orderedGR(.ele), type="data", format="BED")
    })
    names(tracks) <- names
    viewTracks(tracks, gr=range, 
               viewerStyle=trackViewerStyle(),
               autoOptimizeStyle=FALSE,
               newpage=TRUE)
}