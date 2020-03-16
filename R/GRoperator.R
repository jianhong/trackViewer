GRbin <- function(A, B, col, ignore.strand){
    mcols(A) <- DataFrame(source="A", id=seq_along(A), value=mcols(A)[, col])
    mcols(B) <- DataFrame(source="B", id=seq_along(B), value=mcols(B)[, col])
    D <- c(A, B)
    C <- disjoin(D, with.revmap=TRUE, 
                 ignore.strand=ignore.strand)
    revmap <- data.frame(oid=unlist(C$revmap), 
                         nid=rep(seq_along(C), lengths(C$revmap)))
    revmap$source <- D$source[revmap$oid]
    C$revmap <- NULL
    setMeta <- function(group="A"){
        X <- C
        mcols(X) <- DataFrame(source="", id=0, value=0)
        id <- revmap[revmap$source==group, ]
        mcols(X)[id$nid, ] <- mcols(D)[id$oid, ]
        X$id[X$id==0] <- NA
        X
    }
    A1 <- setMeta("A")
    B1 <- setMeta("B")
    GRangesList(A=A1, B=B1)
}

#' GRanges operator
#' @description GRanges operations (add, aubtract, multiply, divide)
#' @param A an object of GRanges
#' @param B an object of GRanges
#' @param col colname of A and B to be calculated
#' @param operator operator, "+" means A + B, and so on. 
#' User-defined function also could be used.
#' @param ignore.strand When set to TRUE, the strand information is 
#' ignored in the overlap calculations.
#' @return an object of GRanges
#' @import GenomicRanges
#' @export
#' @examples 
#' gr2 <- GRanges(seqnames=c("chr1", "chr1"),
#' ranges=IRanges(c(7,13), width=3),
#' strand=c("-", "-"), score=3:4)
#' gr3 <- GRanges(seqnames=c("chr1", "chr1"),
#'                ranges=IRanges(c(1, 4), c(3, 9)),
#'                strand=c("-", "-"), score=c(6L, 2L))
#' GRoperator(gr2, gr3, col="score", operator="+")
#' GRoperator(gr2, gr3, col="score", operator="-")
#' GRoperator(gr2, gr3, col="score", operator="*")
#' GRoperator(gr2, gr3, col="score", operator="/")
#' GRoperator(gr2, gr3, col="score", operator=mean)

GRoperator <- function(A, B, col="score", 
                       operator=c("+", "-", "*", "/", "^", "%%"),
                       ignore.strand=TRUE){
    if(!is(A, "GRanges") || !is(B, "GRanges")){
        stop("A and B must be objects of GRanges")
    }
    if(any(countOverlaps(A)>1) || any(countOverlaps(B)>1)){
        stop("A or B has overlaps.")
    }
    if(length(A)<1) return(B)
    if(length(B)<1) return(A)
    argA <- deparse(substitute(A))
    argB <- deparse(substitute(B))
    if(!is.function(operator)){
        operator <- match.arg(operator)
        operator <- .Primitive(operator)
    }
    
    if(!(col %in% colnames(mcols(A)))){
        stop("col is not in metadata of A")
    }
    if(!(col %in% colnames(mcols(B)))){
        stop("col is not in metadata of B")
    }
    C <- GRbin(A, B, col, ignore.strand=ignore.strand)
    A <- C$A
    B <- C$B
    out <- A
    mcols(out) <- NULL
    out$value <- tryCatch(operator(A$value, B$value),
                          error=function(e) mapply(operator, A$value, B$value))
    out$A <- A$source=="A"
    out$A_id <- A$id
    out$B <- B$source=="B"
    out$B_id <- B$id
    colnames(mcols(out)) <- c(col, 
                              argA, paste0(argA, "_id"),
                              argB, paste0(argB, "_id"))
    out
}

