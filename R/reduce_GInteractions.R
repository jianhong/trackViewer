#' Reduce method for 'GInteractions'
#' @description
#' Reduce returns an object of the same type as x containing reduced ranges for
#' each distinct (seqname, strand) pairing.
#' @importMethodsFrom GenomicRanges reduce range
#' @importFrom InteractionSet anchorIds regions GInteractions
#' @importFrom S4Vectors first second
#' @exportMethod reduce
#' @aliases reduce,GInteractions
#' @param x GInteractions object.
#' @param min.gapwidth Ranges separated by a gap of at least min.gapwidth
#'  positions are not merged.
#' @param ignore.strand	TRUE or FALSE. Whether the strand of the input ranges
#' should be ignored or not.
#' @param ... Not used.
#' @examples
#' \dontrun{
#' library(InteractionSet) 
#' gi <- readRDS(system.file("extdata", "gi.rds", package="trackViewer"))
#' reduce(head(gi, n=20))
#' }
#' 
setMethod("reduce", "GInteractions",
          function(x, min.gapwidth=1L,
                   ignore.strand=TRUE, ...){
            anchor1 <- first(x)
            anchor2 <- second(x)
            an1 <- reduce(anchor1, with.revmap=TRUE,
                          min.gapwidth=min.gapwidth,
                          ignore.strand=ignore.strand)
            an2 <- reduce(anchor2, with.revmap=TRUE,
                          min.gapwidth=min.gapwidth,
                          ignore.strand=ignore.strand)
            il1 <- an1$revmap[lengths(an1$revmap)>1]
            il2 <- an2$revmap[lengths(an2$revmap)>1]
            tobemerged <- mapply(il1[rep(seq_along(il1), each=length(il2))],
                                 il2[rep(seq_along(il2), length(il1))],
                                 FUN=function(a, b){
              intersect(a, b)
            }, SIMPLIFY = FALSE)
            tobemerged <- tobemerged[lengths(tobemerged)>1]
            anchor1Ids <- anchorIds(x, type="first")
            anchor2Ids <- anchorIds(x, type="second")
            reg <- regions(x)
            if(ignore.strand){
              strand(reg) <- '*'
            }
            tobemerged <- lapply(tobemerged, function(.ele){
              c(range(reg[anchor1Ids[.ele]]),
                range(reg[anchor2Ids[.ele]]))
            })
            tobemerged <- unlist(GRangesList(tobemerged))
            tobemerged <- unique(tobemerged)
            tobemerged <- c(tobemerged, reg)
            tobemerged <- reduce(tobemerged, min.gapwidth=0L)
            tobemerged <- sort(tobemerged)
            ol1 <- findOverlaps(anchor1, tobemerged, type = 'within',
                                minoverlap = 1L)
            ol2 <- findOverlaps(anchor2, tobemerged, type = 'within',
                                minoverlap = 1L)
            anchor1[queryHits(ol1)] <- tobemerged[subjectHits(ol1)]
            anchor2[queryHits(ol2)] <- tobemerged[subjectHits(ol2)]
            newGI <- GInteractions(anchor1, anchor2, score=mcols(x)$score)
            newGI <- unique(newGI)
            newGI[order(seqnames(first(newGI)),
                        start(first(newGI)),
                        seqnames(second(newGI)),
                        start(second(newGI)),
                        decreasing = FALSE)]
          }
)