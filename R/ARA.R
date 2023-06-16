#' Aggregate Region Analysis
#' @description Extract the interaction signal means from given coordinates. 
#' @param gr A `GRanges` object. The center of the object will be used for
#'  alignment for all the given regions.
#' @param upstream,downstream numeric(1L). Upstream and downstream from the 
#' center of given `gr` input will be used to extract the signals.
#' @param resolution numeric(1L). The resolution will be passed to
#'  \link{importGInteractions} function.
#' @param ... The parameters used by \link{importGInteractions} function.
#' Please note that the ranges resolution and out parameter should not be
#' involved.
#' @return A \link[InteractionSet:GInteractions-class]{GInteractions} object
#' with scores which represent the mean values of the interactions.
#' @export
#' @importFrom InteractionSet regions `regions<-`
#' @examples
#' hic <- system.file("extdata", "test_chr22.hic", package = "trackViewer",
#'                    mustWork=TRUE)
#' gr <- GRanges("22", c(seq(20000001, 50000001, by=10000000), width=1))
#' gi <- ARA(gr, file=hic, format="hic")
#' rg <- GRanges("22", IRanges(1, 400000))
#' op <- optimizeStyle(trackList(gi2track(gi)))
#' heatmap <- op$tracks
#' sty <- op$style
#' setTrackViewerStyleParam(sty, "xat", c(1, 200000, 400000))
#' setTrackViewerStyleParam(sty, "xlabel",c("-20K", "center", "20K"))
#' viewTracks(heatmap, viewerStyle=sty, gr=rg)

ARA <- function(gr, upstream=200000, downstream=upstream, resolution=10000,
                ...){
  stopifnot(is(gr, "GRanges"))
  stopifnot(is(upstream, "numeric"))
  stopifnot(is(downstream, "numeric"))
  stopifnot(is(resolution, "numeric"))
  resolution <- resolution[1]
  start(gr) <- floor((start(gr)+end(gr))/2/resolution)*resolution+1
  width(gr) <- resolution
  gr <- promoters(gr, upstream = upstream, downstream = downstream)
  gr <- gr[start(gr)>0]
  gr <- trim(gr)
  gr <- gr[width(gr)==upstream + downstream]
  gi <- lapply(seq_along(gr), function(.ele){
    importGInteractions(ranges=gr[.ele], resolution = resolution,
                        out = "GInteractions", ...)
  })
  seqn <- as.character(seqlevels(gr))[1]
  gi_shift <- mapply(gi, start(gr), FUN = function(.ele, start){
    .ele <- shift(.ele, -1*start)
    seqlevels(regions(.ele)) <- seqn
    .ele
  })
  gi_accum <- GIoperator(gi_shift, operator = "+")
  gi_accum$score <- gi_accum$score/length(gi_accum)
  gi_accum
}