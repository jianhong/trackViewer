#' plot tracks based on gene name
#' @description given a gene name, plot the tracks.
#' @param symbol Gene symbol
#' @param upstream upstream from anchor
#' @param downstream downstream from anchor
#' @param anchor TSS, or gene
#' @param filenames files used to generate tracks
#' @param format file format used to generate tracks
#' @param txdb txdb will be used to extract the genes
#' @param org org package name
#' @param plot plot the tracks or not.
#' @return an invisible list of a \code{\link{trackList}},
#'  a \code{\link{trackViewerStyle}} and a \code{\link[GenomicRanges]{GRanges}}
#' @import GenomicRanges
#' @import GenomicFeatures
#' @import Rsamtools
#' @importFrom BiocGenerics basename
#' @export
#' @examples 
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' library(org.Hs.eg.db)
#' extdata <- system.file("extdata", package="trackViewer", mustWork=TRUE)
#' filename = file.path(extdata, "fox2.bed")
#' optSty <- viewGene("HSPA8", filenames=filename, format="BED", 
#'                    txdb=TxDb.Hsapiens.UCSC.hg19.knownGene, 
#'                    org="org.Hs.eg.db")
#'


viewGene <- function(symbol, filenames, format, txdb, org,
                     upstream=1000, downstream=1000, 
                     anchor=c("gene", "TSS"),
                     plot=FALSE){
  anchor <- match.arg(anchor)
  stopifnot(is(txdb, "TxDb"))
  gr <- getLocation(symbol, txdb, org)
  if(length(gr)>1){
    stop("Can not handle multiple positions.")
  }
  if(anchor=="TSS"){
    gr <- promoters(gr, upstream = upstream, downstream = downstream)
  }else{
    if(as.character(strand(gr))=="-"){
      tmp <- upstream
      upstream <- downstream
      downstream <- upstream
    }
    start(gr) <- start(gr) - upstream
    end(gr) <- end(gr) + downstream
  }
  trs <- geneModelFromTxdb(txdb, get(org), gr=gr)
  if(format=="BAM"){
    data <- sapply(filenames, function(filename){
      importBam(filename, ranges = gr, 
                pairs = testPairedEndBam(filename))
    })
  }else{
    data <- sapply(filenames, importScore, format = format, 
                   ranges=gr, ignore.strand=TRUE)
  }
  names(data) <- basename(filenames)
  tL <- trackList(trs, data, heightDist=c(length(trs), length(data)))
  optSty <- optimizeStyle(trackList = tL)
  trackList <- optSty$tracks
  viewerStyle <- optSty$style
  optSty$range <- gr
  if(plot) viewTracks(trackList, gr=gr, viewerStyle = viewerStyle)
  return(invisible(optSty))
}
