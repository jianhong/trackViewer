viewTracks <- function(trackList, chromosome, start, end, strand, gr=GRanges(),
                       viewerStyle=trackViewerStyle(), autoOptimizeStyle=FALSE,
                       newpage=TRUE){
    if(missing(trackList)){
        stop("trackList is required.")
    }
    if(class(trackList)!="trackList" && 
           !((is.list(trackList) && all(sapply(trackList, class)=="track")))){
        stop("trackList must be an object of \"trackList\"
             (See ?trackList) or a list of track")
    }
    if(missing(gr)){
        if(missing(chromosome) || missing(start) || missing(end))
            stop("Please input the coordinate.")
        if(missing(strand) || !strand %in% c("+", "-"))
            strand <- "*"
    }else{
        if(class(gr)!="GRanges" || length(gr)!=1){
            stop("gr must be an object of GRanges with one element.")
        }
        chromosome <- as.character(GenomicRanges::seqnames(gr))[1]
        start <- GenomicRanges::start(gr)[1]
        end <- GenomicRanges::end(gr)[1]
        strand <- as.character(GenomicRanges::strand(gr))[1]
    }
    
    if(class(viewerStyle)!="trackViewerStyle"){
        stop("viewerStyle must be an object of 'trackViewerStyle'.")
    }
    if(!is.logical(newpage)) stop("newpage should be a logical vector.")
    
    trackList <- filterTracks(trackList, chromosome, start, end, strand)
    
    if(newpage) grid.newpage()
    
    if(autoOptimizeStyle){
        opt <- optimizeStyle(trackList, viewerStyle)
        trackList <- opt$tracks
        viewerStyle <- opt$style
    }
    
    margin <- viewerStyle@margin
    
    ##xscale, yscale
    xscale <- c(start, end)
    
    yscales <- getYlim(trackList)
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
    for(i in 1:total){
        plotTrack(names(trackList)[i], trackList[[i]], 
                  viewerStyle, ht,
                  yscales[[i]], yHeights[i], xscale,
                  chromosome, strand)
        ht <- ht + yHeights[i]
#        if(interactive()) setTxtProgressBar(pb, i)
    }
#    if(interactive()) close(pb)
    popViewport()
    
    return(invisible(viewport(x=margin[2], y=margin[1], 
                              height=1 - margin[1]- margin[3], 
                              width=1 -margin[2] - margin[4],
                              just=c(0,0), xscale=xscale, yscale=c(0,1))))
}