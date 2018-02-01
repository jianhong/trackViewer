#' convert WIG format track to BED format track
#' @description convert WIG format track to BED format track for a given range
#' @param trackScore an object of track with WIG format
#' @param chrom sequence name of the chromosome
#' @param from start coordinate
#' @param to end coordinate
#' @return an object of \code{\link{track}}
#' @import GenomicRanges
#' @export
#' @examples 
#' extdata <- system.file("extdata", package="trackViewer", mustWork=TRUE)
#' repA <- importScore(file.path(extdata, "cpsf160.repA_-.wig"),
#'                     file.path(extdata, "cpsf160.repA_+.wig"),
#'                     format="WIG")
#' strand(repA$dat) <- "-"
#' strand(repA$dat2) <- "+"
#' parseWIG(repA, chrom="chr11", from=122929275, to=122930122)

parseWIG <- function(trackScore, chrom, from, to){
    if(!is(trackScore, "track"))
        stop("trackScore must be an object of track")
    if(trackScore@format!="WIG")
        stop("format must be WIG")
    
    parser <- function(data, chrom, from, to){
        data <- data[seqnames(data)==chrom]
        if(length(data)<1) return(GRanges(score=numeric(0)))
        res <- mapply(function(chr, start, data, span, step, 
                               structure, from, to, strand){
            if(length(data)>1000){
                dataList <- split(data, 
                                  rep(1:ceiling(length(data)/1000), each=1000)[length(data)])
            }else{
                dataList <- list(data)
            }
            gr <- GRanges(score=numeric(0))
            for(data in dataList){
                data <- do.call(rbind, strsplit(data, "\\s+"))
                if(structure=="variableStep"){
                    pos1 <- as.numeric(data[, 1])
                    if(is.na(span)) span <- 1
                    pos2 <- pos1 + as.numeric(span) - 1
                    score <- as.numeric(data[, 2])
                }else{##structure=="fixedStep"
                    pos1 <- start + (1:nrow(data) - 1) * step
                    if(!is.na(span)){
                        pos2 <- pos1 + span
                    }else{
                        pos2 <- start + (1:nrow(data)) * step - 1
                    }
                    score <- as.numeric(data[, 1])
                }
                idx <- pos2>=from & pos1<=to
                if(sum(idx)>0){
                    gr <- c(gr, GRanges(seqnames=chr, 
                                        ranges=IRanges(start=pos1[idx],
                                                       end=pos2[idx]),
                                        strand=strand,
                                        score=score[idx]))
                }
            }
            gr
        }, as.character(seqnames(data)), start(data), data$score, 
                      data$span, data$step, data$structure, from, to, 
                      as.character(strand(data)))
        res <- unlist(GRangesList(res))
        res <- split(res, res$score)
        res <- lapply(res, reduce)
        res <- mapply(function(.ele, .sc) {
            .ele$score <- .sc
            .ele}, res, as.numeric(names(res)))
        res <- unlist(GRangesList(res))
        names(res) <- NULL
        if(length(res)<1) res <- GRanges(score=numeric(0))
        return(orderedGR(res))
    }
    trackScore@dat <- parser(trackScore@dat, chrom, from, to)
    gc(reset=TRUE)
    trackScore@dat2 <- parser(trackScore@dat2, chrom, from, to)
    gc(reset=TRUE)
    trackScore@format <- "BED"
    return(trackScore)
}