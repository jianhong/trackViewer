#' Reading data from a ginteractions, hic, cool, or validPairs file
#' @description Read a \code{\link{track}} object from a ginteractions, hic, cool, or validPairs file
#' @param file The path to the file to read.
#' @param format The format of import file. Could be ginteractions, hic, cool or validPairs
#' @param ranges An object of \code{\link[GenomicRanges:GRanges-class]{GRanges}} to indicate
#' the range to be imported
#' @param ignore.strand ignore the strand or not when do filter. default TRUE
#' @return a \code{\link{track}} object 
#' @import GenomicRanges
#' @import IRanges
#' @importFrom InteractionSet GInteractions
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
#' #hic <- system.file("extdata", "test.hic", package = "trackViewer",
#' #                       mustWork=TRUE)
#' #dat <- importGInteractions(file=hic, format="hic")
#' 
#' ##import a cool file
#' #cool <- system.file("extdata", "test.cool", package = "trackViewer",
#' #                       mustWork=TRUE)
#' #dat <- importGInteractions(file=cool, format="cool")
#' 
#' ##import a validPairs file
#' #validPairs <- system.file("extdata", "test.validPairs", package = "trackViewer",
#' #                       mustWork=TRUE)
#' #dat <- importGInteractions(file=validPairs, format="validPairs")
#' 

importGInteractions <- function(file, 
                        format=c("ginteractions", "hic", "cool", "validPairs"),
                        ranges=GRanges(), ignore.strand=TRUE){
    if(missing(file))
        stop("file is required.")
    ##    on.exit(closeAllConnections())
    format <- match.arg(format)
    #    res <- GRanges(score=numeric(0))
    if(!is(ranges, "GRanges")) stop("ranges must be an object of GRanges.")
    gr <- orderedGR(ranges)
    readginteractions <- function(buf){
        buf <- do.call(rbind, strsplit(buf, "\\s"))
        dat <- GRanges(buf[, 1], IRanges(as.numeric(buf[, 2]), 
                                         as.numeric(buf[, 3])))
        dat2 <- GRanges(buf[, 4], IRanges(as.numeric(buf[, 5]), 
                                          as.numeric(buf[, 6])))
        gi <- GInteractions(anchor1=dat, anchor2 = dat2, 
                            score = as.numeric(buf[, 7]))
        subsetByOverlaps(gi, gr)
    }
    readhic <- function(buf){
        stop("not support yet.")
    }
    readcool <- function(buf){
        stop("not support yet.")
    }
    readvalidPairs <- function(buf){
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
    readFile <- function(file, format, FUN){
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
        res <- readFile(file, format, FUN)
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
    return(gi2track(res))
}
