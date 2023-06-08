#' GInteractions operator
#' @description GInteractions operations (add, aubtract, multiply, divide)
#' @param gi_list a list of GInteractions objects
#' @param col colname of metadata to be calculated
#' @param operator operator, "+" means A + B, and so on. 
#' User-defined function also could be used.
#' @return an object of GInteractions
#' @export
#' @examples
#' library(InteractionSet)
#' gr2 <- GRanges(seqnames=c("chr1", "chr1"),
#'                ranges=IRanges(c(7,13), width=3))
#' gr3 <- GRanges(seqnames=c("chr1", "chr1"),
#'                ranges=IRanges(c(1, 4), c(3, 9)))
#' gi <- GInteractions(gr2, gr3, score=c(1, 2))
#' gi2 <- GInteractions(gr2, gr3, score=c(3, 4))
#' GIoperator(list(gi, gi2), col="score", operator="+")
#' GIoperator(list(gi, gi2), col="score", operator="-")

GIoperator <- function(gi_list, col="score",
                       operator=c("+", "-", "*", "/")){
  null <- lapply(gi_list, function(.ele){
    stopifnot(is(.ele, "GInteractions"))
    stopifnot("col in not in metadata of input" =
                col %in% colnames(mcols(.ele)))
  })
  if(!is.function(operator)){
    operator <- match.arg(operator)
    operator <- .Primitive(operator)
  }
  gi_add <- unique(Reduce(c, gi_list))
  mcols(gi_add)[, col] <- 0
  gi_add <- sort(gi_add)
  ol <- findOverlaps(gi_add, gi_list[[1]], type="equal")
  mcols(gi_add)[queryHits(ol), col] <- mcols(gi_list[[1]])[subjectHits(ol), col]
  for(i in seq_along(gi_list)[-1]){
    ol <- findOverlaps(gi_add, gi_list[[i]], type="equal")
    mcols(gi_add)[queryHits(ol), col] <- 
      do.call(operator,
              args=list(mcols(gi_add)[queryHits(ol), col], 
                        mcols(gi_list[[i]])[subjectHits(ol), col]))
  }
  gi_add
}
