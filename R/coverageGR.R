#' calculate coverage
#' @description calculate coverage for \code{\link[GenomicRanges:GRanges-class]{GRanges}},
#'  \code{\link[GenomicAlignments:GAlignments-class]{GAlignments}} or 
#'  \code{\link[GenomicAlignments:GAlignmentPairs-class]{GAlignmentPairs}}
#' @param gr an object of RGanges, GAlignments or GAlignmentPairs
#' @return an object of GRanges
#' @seealso See Also as \code{\link[IRanges:coverage-methods]{coverage}}, 
#' \code{\link[GenomicRanges]{coverage-methods}}
#' @export
#' @import GenomicAlignments
#' @import GenomicRanges
#' @import IRanges
#' @examples 
#' bed <- system.file("extdata", "fox2.bed", package="trackViewer",
#'                    mustWork=TRUE)
#' fox2 <- importScore(bed)
#' fox2$dat <- coverageGR(fox2$dat)
#' 
coverageGR <- function(gr){##todo, check GAlignmentPairs
    if(!inherits(gr, c("GRanges", "GAlignments", "GAlignmentPairs")))
        stop("input is not an object of RGanges, GAlignments or GAlignmentPairs.")
    flag <- TRUE
    if(is(gr, "GRanges")){
        if(is.null(score(gr))) gr$score <- rep(1, length(gr))
        if(length(unique(score(gr)))!=1) flag <- FALSE
    }
    if(flag){
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
    }else{
        disj <- disjoin(gr)
        if(length(disj)>length(gr)){
            ol <- findOverlaps(disj, gr)
            score <- sapply(split(subjectHits(ol), queryHits(ol)), 
                            function(.ids){
                                sum(score(gr)[.ids])
                            })
            disj$score <- rep(0, length(disj))
            disj$score[as.numeric(names(score))] <- score
            ##merge continuous region
            block <- rle(score(disj))
            if(any(block$lengths!=1)){
                block$values <- 1:length(block$values)
                block <- inverse.rle(block)
                dup <- block %in% block[duplicated(block)]
                disj.uni <- disj[!dup]
                disj.dup <- disj[dup]
                block <- rle(score(disj.dup))
                block$values <- 1:length(block$values)
                block <- inverse.rle(block)
                disj.dup <- split(disj.dup, block)
                
                disj.rd <- sapply(disj.dup, reduce)
                disj.rd <- unlist(GRangesList(disj.rd))
                disj.rd$score <- sapply(disj.dup, function(.ele) score(.ele)[1])
                disj <- c(disj.uni, disj.rd)
            }
            gr <- unname(disj)
        }
    }
    return(orderedGR(gr))
}