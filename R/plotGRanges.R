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