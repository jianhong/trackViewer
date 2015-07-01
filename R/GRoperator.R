GRbin <- function(A, B, col){
    A1 <- B1 <- C <- disjoin(c(A, B), ignore.strand=FALSE)
    ol1 <- findOverlaps(C, A, maxgap=0L, minoverlap=1L,
                        type="any", ignore.strand=FALSE)
    ol2 <- findOverlaps(C, B, maxgap=0L, minoverlap=1L,
                        type="any", ignore.strand=FALSE)
    A1$value <- B1$value <- 0
    A1$value[queryHits(ol1)] <- mcols(A)[subjectHits(ol1), col]
    B1$value[queryHits(ol2)] <- mcols(B)[subjectHits(ol2), col]
    GRangesList(A=A1, B=B1)
}
GRoperator <- function(A, B, col="score", operator=c("+", "-", "*", "/", "^", "%%")){
    if(class(A)!="GRanges" || class(B)!="GRanges"){
        stop("A and B must be objects of GRanges")
    }
    operator <- match.arg(operator)
    if(!(col %in% colnames(mcols(A)))){
        stop("col is not in metadata of A")
    }
    if(!(col %in% colnames(mcols(B)))){
        stop("col is not in metadata of B")
    }
    C <- GRbin(A, B, col)
    A <- C$A
    B <- C$B
    out <- A
    operator <- .Primitive(operator)
    out$value <- operator(A$value, B$value)
    colnames(mcols(out)) <- col
    out
}

