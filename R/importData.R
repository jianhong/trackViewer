#' Reading data from a BED or WIG file to RleList
#' @description Read a \code{\link{track}} object from a BED, bedGraph, 
#' WIG or BigWig file to RleList
#' @param files The path to the files to read.
#' @param format The format of import file. Could be BAM, BED, bedGraph, WIG or BigWig
#' @param ranges An object of \code{\link[GenomicRanges:GRanges-class]{GRanges}} to indicate
#' the range to be imported
#' @return a list of \code{\link[IRanges:AtomicList]{RleList}}. 
#' @import Rsamtools
#' @importFrom tools file_ext
#' @export
#' @examples 
#' #import a BED file
#' bedfile <- system.file("tests", "test.bed", package="rtracklayer",
#'                        mustWork=TRUE)
#' dat <- importData(files=bedfile, format="BED",
#'                   ranges=GRanges("chr7", IRanges(127471197, 127474697)))
#' 
#' ##import a WIG file
#' wigfile <- system.file("tests", "step.wig", package = "rtracklayer",
#'                        mustWork=TRUE)
#' dat <- importData(files=wigfile, format="WIG", 
#'                   ranges=GRanges("chr19", 
#'                                  IRanges(59104701, 59110920)))
#' 
#' ##import a BigWig file
#' if(.Platform$OS.type!="windows"){
#'   ##this is because we are using rtracklayer::import
#'   bwfile <- system.file("tests", "test.bw", package = "rtracklayer",
#'                         mustWork=TRUE)
#'   dat <- importData(files=bwfile, format="BigWig", 
#'                     ranges=GRanges("chr19", IRanges(1500, 2700)))
#' }

importData <- 
    function(files, 
             format=NA,
             ranges=GRanges()){
    stopifnot(class(ranges)=="GRanges")
    ranges <- reduce(ranges)
    formats <- c("BAM", "BED", "bedGraph", "WIG", "BigWig")
    if(is.na(format[1])){
        format <- sapply(files, file_ext) ## need tools package
    }
    format <- rep(format, length(files))[1:length(files)]
    format <- formats[match(tolower(format), tolower(formats))]
    if(any(is.na(format)))
        stop("Can not determine the file format.",
             "Format must be BAM, BED, bedGraph, WIG, or BigWig")
    covBAM <- NULL
    covBW <- NULL
    covBED <- NULL
    covBedGraph <- NULL
    covWIG <- NULL
    bams <- files[format=="BAM"]
    BigWigs <- files[format=="BigWig"]
    BEDs <- files[format=="BED"]
    bedGraphs <- files[format=="bedGraph"]
    WIGs <- files[format=="WIG"]
    if(length(bams)>0){
        ## import bams, need Rsamtools, need GenomicAlignments
        param <- ScanBamParam(which=ranges)
        bam <- lapply(bams, readGAlignments, param=param)
        names(bam) <- bams
        covBAM <- sapply(bam, coverage)
    }
    if(length(BigWigs)>0){
        ## import by rtracklayer
        ### if format == "BigWig"
        ## rtracklayer can't read bigWig files on a Windows computer. 
        ## Type in  ?`BigWigFile-class` to get the help.
        if (.Platform$OS.type != "windows") {
            covBW <- sapply(BigWigs, import, 
                            format="BigWig", 
                            which=ranges, 
                            as="RleList")
        }else{
            stop("Can not import the bigWig files on a Windows computer.")
        }
    }
    if(length(BEDs)>0){
        covBED <- sapply(BEDs, importScore, 
                         format="BED",
                         ranges=ranges,
                         ignore.strand=TRUE)
        covBED <- sapply(covBED, function(.ele){
            coverage(.ele@dat, weight=score(.ele@dat))
        })
    }
    if(length(bedGraphs)>0){
        covBedGraph <- sapply(bedGraphs, importScore, 
                              format="bedGraph",
                              ranges=ranges,
                              ignore.strand=TRUE)
        covBedGraph <- sapply(covBedGraph, function(.ele){
            coverage(.ele@dat, weight=score(.ele@dat))
        })
    }
    if(length(WIGs)>0){
        ranges <- reduce(ranges, ignore.strand=TRUE)
        covWIG <- sapply(WIGs, importScore, 
                              format="WIG",
                              ranges=ranges,
                              ignore.strand=TRUE)
        covWIG <- sapply(covWIG, function(.ele){
            .tracks <- lapply(seq_len(length(ranges)), function(id){
                .gr <- ranges[id]
                .out <- parseWIG(trackScore=.ele,
                                 chrom=as.character(seqnames(.gr)),
                                 from=start(.gr),
                                 to=end(.gr))
                .out@dat
            })
            .tracks <- unlist(GRangesList(.tracks))
            coverage(.tracks, weight=score(.tracks))
        })
    }
    covs <- c(covBAM, covBW, covBED, covBedGraph, covWIG)
    covs[files]
}