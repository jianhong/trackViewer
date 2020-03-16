#' Reading data from a BAM file
#' @description Read a \code{\link{track}} object from a BAM file
#' @param file The path to the BAM file to read.
#' @param file2 The path to the second BAM file to read.
#' @param ranges An object of \code{\link[GenomicRanges:GRanges-class]{GRanges}} to indicate
#' the range to be imported
#' @param pairs logical object to indicate the BAM is paired or not. See
#' \code{\link[GenomicAlignments]{readGAlignments}}
#' @return a \code{\link{track}} object 
#' @import GenomicAlignments
#' @export
#' @seealso See Also as \code{\link{importScore}}, \code{\link{track}}, 
#' \code{\link{viewTracks}}
#' @examples 
#' bamfile <- system.file("extdata", "ex1.bam", package="Rsamtools",
#' mustWork=TRUE)
#' dat <- importBam(file=bamfile, ranges=GRanges("seq1", IRanges(1, 50), strand="+"))

importBam <- function(file, file2, ranges=GRanges(), pairs=FALSE){
    if(!is(ranges, "GRanges"))
        stop("ranges must be an object of GRanges.")
    param <- ScanBamParam(which=ranges)
    stopifnot(length(pairs)==1)
    if(pairs){
        bam <- readGAlignmentPairs(file, param=param)
    }else{
        bam <- readGAlignments(file, param=param)
    }
    gr <- coverageGR(bam)
    if(!missing(file2)){
        if(pairs){
            bam2 <- readGAlignmentPairs(file2, param=param)
        }else{
            bam2 <- readGAlignments(file2, param=param)
        }
        gr2 <- coverageGR(bam2)
        return(new("track", dat=orderedGR(gr), dat2=orderedGR(gr2),
                   type="data", format="BAM"))
    }
    return(new("track", dat=orderedGR(gr), type="data", format="BAM"))
}
