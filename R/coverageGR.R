###coverage
coverageGR <- function(gr){##todo, check GAlignmentPairs
    if(!class(gr) %in% c("GRanges", "GAlignments", "GAlignmentPairs"))
        stop("input is not an object of RGanges, GAlignments or GAlignmentPairs.")
    gr <- split(gr, strand(gr))
    gr <- mapply(function(.ele, strand){
        coverage <- coverage(.ele)
        re <- mapply(function(cov, seqname, strand){
            len <- length(cov)
            if(len>0){
                width <- cov@lengths
                start <- cumsum(width)
                start <- c(1, start[-length(start)]+1)
                value <- cov@values
            }else{
                start <- 1
                width <- 1
                value <- 0
            }
            idx <- value > 0
            if(all(!idx)) {
                GRanges(score=numeric(0))
            }else{
                GRanges(seqnames=seqname, 
                        ranges=IRanges(start=start[idx],
                                       width=width[idx]),
                        strand=strand,
                        score=value[idx])
            }
        }, coverage, names(coverage), strand)
        unlist(GRangesList(re))
    }, gr, names(gr))
    gr <- unlist(GRangesList(gr))
    names(gr) <- NULL
    return(orderedGR(gr))
}