importBam <- function(file, file2, ranges=GRanges(), pairs=FALSE){
    if(class(ranges)!="GRanges")
        stop("ranges must be an object of GRanges.")
    param <- ScanBamParam(which=ranges)
    if(pairs){
        bam <- readGAlignmentPairsFromBam(file, param=param)
    }else{
        bam <- readGAlignmentsFromBam(file, param=param)
    }
    gr <- coverageGR(bam)
    if(!missing(file2)){
        if(pairs){
            bam2 <- readGAlignmentPairsFromBam(file2, param=param)
        }else{
            bam2 <- readGAlignmentsFromBam(file2, param=param)
        }
        gr2 <- coverageGR(bam2)
        return(new("track", dat=orderedGR(gr), dat2=orderedGR(gr2),
                   type="data", format="BAM"))
    }
    return(new("track", dat=orderedGR(gr), type="data", format="BAM"))
}