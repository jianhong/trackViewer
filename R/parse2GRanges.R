#' parse text into GRanges
#' @description parse text like "chr13:99,443,451-99,848,821:-" into GRanges
#' @param text character vector like "chr13:99,443,451-99,848,821:-" or 
#' "chr13:99,443,451-99,848,821"
#' @return an object of \link[GenomicRanges:GRanges-class]{GRanges}
#' @import GenomicRanges
#' @export
#' @examples 
#' parse2GRanges("chr13:99,443,451-99,848,821:-")
#' 

parse2GRanges <- function(text){
    ## chr13:99,443,451-99,848,821:-
    x <- gsub(",", "", text)
    x <- do.call(rbind, strsplit(x, ":"))
    if(ncol(x)==1){
        stop("Can not handle the text.")
    }
    coor <- do.call(rbind, strsplit(x[,2], "\\-"))
    if(ncol(coor)!=2){
        stop("Can not handle the text.")
    }
    if(ncol(x)==2){
        x <- cbind(x, "*")
    }
    GRanges(x[,1], 
            IRanges(as.numeric(coor[,1]), as.numeric(coor[,2])),
            strand=x[,3])
}