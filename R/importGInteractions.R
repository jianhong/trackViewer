#' Reading data from a ginteractions, hic, cool, or validPairs file
#' @description Read a \code{\link{track}} object from a ginteractions, hic, mcool, or validPairs file
#' @param file The path to the file to read.
#' @param format The format of import file. Could be ginteractions, hic, cool or validPairs
#' @param ranges An object of \code{\link[GenomicRanges:GRanges-class]{GRanges}} to indicate
#' the range to be imported
#' @param ignore.strand ignore the strand or not when do filter. default TRUE
#' @param out output format. Default is track. Possible values: track, GInteractions.
#' @param resolution Resolutions for the interaction data.
#' @param unit BP (base pair) or FRAG (fragment) (.hic file only). 
#' @param normalization Type of normalization, NONE, VC, VC_SORT or KR for .hic and
#' NONE, balanced for .cool.
#' @param ... NOT used.
#' @return a \code{\link{track}} object 
#' @import GenomicRanges
#' @import IRanges
#' @import Rcpp
#' @importFrom InteractionSet GInteractions
#' @importFrom rhdf5 H5Fopen h5ls h5read h5closeAll
#' @useDynLib trackViewer
#' @export
#' @seealso See Also as \code{\link{importBam}}, \code{\link{track}}, 
#' \code{\link{viewTracks}}
#' @examples 
#' #import a ginteractions file
#' #gi <- system.file("extdata", "test.ginteractions.tsv", package="trackViewer",
#' #                       mustWork=TRUE)
#' #dat <- importGInteractions(file=gi, format="ginteractions",
#' #                   ranges=GRanges("chr7", IRanges(127471197, 127474697)))
#' 
#' ##import a hic file
#' if(.Platform$OS.type!="windows"){
#' hic <- system.file("extdata", "test_chr22.hic", package = "trackViewer",
#'                        mustWork=TRUE)
#' dat <- importGInteractions(file=hic, format="hic", 
#'                            ranges=GRanges("22", IRanges(1500000, 100000000)))
#' }
#' 
#' ##import a cool file
#' cool <- system.file("extdata", "test.mcool", package = "trackViewer",
#'                      mustWork=TRUE)
#' dat <- importGInteractions(file=cool, format="cool",
#'                            resolution = 2,
#'                            ranges=GRanges("chr1", IRanges(10, 28)))
#' 
#' ##import a validPairs file
#' #validPairs <- system.file("extdata", "test.validPairs", package = "trackViewer",
#' #                       mustWork=TRUE)
#' #dat <- importGInteractions(file=validPairs, format="validPairs")
#' 

importGInteractions <- function(file, 
                        format=c("ginteractions", "hic", "cool", "validPairs"),
                        ranges=GRanges(), ignore.strand=TRUE,
                        out=c("track", "GInteractions"),
                        resolution = 100000,
                        unit=c("BP", "FRAG"),
                        normalization =c("NONE", "VC", "VC_SORT", "KR", "balanced"),
                        ...){
    if(missing(file))
        stop("file is required.")
    ##    on.exit(closeAllConnections())
    format <- match.arg(format)
    out <- match.arg(out)
    unit <- match.arg(unit)
    normalization <- match.arg(normalization)
    
    #    res <- GRanges(score=numeric(0))
    if(!is(ranges, "GRanges")) stop("ranges must be an object of GRanges.")
    gr <- orderedGR(ranges)
    readginteractions <- function(file){
        localFUN = function(buf){
            buf <- do.call(rbind, strsplit(buf, "\\s"))
            dat <- GRanges(buf[, 1], IRanges(as.numeric(buf[, 2]), 
                                             as.numeric(buf[, 3])))
            dat2 <- GRanges(buf[, 4], IRanges(as.numeric(buf[, 5]), 
                                              as.numeric(buf[, 6])))
            gi <- GInteractions(anchor1=dat, anchor2 = dat2, 
                                score = as.numeric(buf[, 7]))
            subsetByOverlaps(gi, gr)
        }
        readFile(file, localFUN);
    }
    readhic <- function(file){
        if(length(gr)>2){
            stop("length of gr > 2 is not supported for .hic file.")
        }
        if(length(gr)==0){
            return(GInteractions())
        }
        if(length(gr)==1){
            seq1 <- seq2 <- as.character(seqnames(gr))
            start1 <- start2 <- start(gr)
            end1 <- end2 <- end(gr)
        }
        if(length(gr)==2){
            seq1 <- as.character(seqnames(gr)[1])
            start1 <- start(gr)[1]
            end1 <- end(gr)[1]
            seq2 <- as.character(seqnames(gr)[2])
            start2 <- start(gr)[2]
            end2 <- end(gr)[2]
        }
        df <- getContactRecords(
            hicfilename = file,
            qname1 = seq1,
            start1 = start1,
            end1 = end1,
            qname2 = seq2,
            start2 = start2,
            end2 = end2,
            binSize = resolution,
            normalization = normalization,
            unit = unit)
        if(length(df)==0){
            return(GInteractions())
        }
        anchor1 <- GRanges(df$seq1, IRanges(df$start1+1, df$end1))
        anchor2 <- GRanges(df$seq2, IRanges(df$start2+1, df$end2))
        gi <- GInteractions(anchor1=anchor1, anchor2=anchor2, 
                            score=df$score)
    }
    readcool <- function(file){
        gi <- cooler_pixels(file, resolution, gr)
        if(normalization=="NONE"){
            gi$score <- gi$count
        }else{
            gi$score <- gi$balanced
        }
        gi$count <- NULL
        gi$balanced <- NULL
        gi
    }
    readvalidPairs <- function(file){
        localFUN = function(buf){
            buf <- do.call(rbind, strsplit(buf, "\\s"))
            dat <- GRanges(buf[, 2], IRanges(as.numeric(buf[, 3]), 
                                             width = 1),
                           strand = buf[, 4])
            dat2 <- GRanges(buf[, 5], IRanges(as.numeric(buf[, 6]), 
                                              width = 1),
                            strand = buf[, ])
            gi <- GInteractions(anchor1=dat, anchor2 = dat2, 
                                score = as.numeric(buf[, 7]))
            subsetByOverlaps(gi, gr)
        }
        readFile(file, localFUN);
    }
    readFile <- function(file, FUN){
        s <- file.info(file)$size
        if(s<100000000){
            buf <- readChar(file, s, useBytes=TRUE)
            buf <- strsplit(buf, "\n", fixed=TRUE, useBytes=TRUE)[[1]]
            res <- FUN(buf)
        }else{
            message("file is too huge. Please consider to subset the data before import.")
            res <- NULL
            con <- file(file, open="r")
            while(length(buf <- readLines(con, n=1000000, warn=FALSE))>0){
                buf <- FUN(buf)
                if(length(res)<1) {
                    res <- buf
                }else{
                    suppressWarnings(res <- c(res, buf))
                }
            }
            close(con)
        }
        res <- unique(res)
        return(res)
    }
    readFiles <- function(file, format){
        FUN <- get(paste("read", format, sep=""))
        res <- FUN(file)
        return(res)
    }
    res <- readFiles(file, format)
    if(any(duplicated(res))){
        res_uniq <- res[!duplicated(res)]
        res_dup <- res[duplicated(res)]
        ol <- findOverlaps(res_uniq, res_dup, type = "equal")
        qol <- split(res_dup[subjectHits(ol)]$score, queryHits(ol))
        qol <- vapply(qol, FUN=sum, FUN.VALUE = 0)
        res_uniq[as.numeric(names(qol))]$score <- 
            res_uniq[as.numeric(names(qol))]$score + qol
        res <- res_uniq
    }
    if(out=="track"){
        return(gi2track(res))
    }else{
        return(res)
    }
}
