#' load ideogram from UCSC
#' @description Download ideogram table from UCSC 
#' @param genome Assembly name  assigned by UCSC, such as hg38, mm10.
#' @param chrom A character vector of chromosome names, or NULL.
#' @param ranges A \link[IRanges:IRanges-class]{Ranges} object with the intervals.
#' @param ... Additional arguments to pass to the 
#' \link[GenomicRanges:GRanges-class]{GRanges} constructor.
#' @importFrom rtracklayer import getTable ucscTableQuery GRangesForUCSCGenome browserSession
#' @importFrom GenomeInfoDb "genome<-"
#' @export
#' @return A \link[GenomicRanges:GRanges-class]{GRanges} object.
#' @seealso See Also as \code{\link{ideogramPlot}}
#' @examples 
#' \dontrun{
#' head(loadIdeogram("hg38", chrom = "chr1"))
#' }
#' 
#' 
loadIdeogram <- function(genome, chrom = NULL, ranges = NULL, ...){
  range <- GRangesForUCSCGenome(genome, chrom, ranges, ...)
  session <- browserSession()
  genome(session) <- genome
  ideo <- getTable(ucscTableQuery(session, 
                          table="cytoBandIdeo", 
                          range = range))
  ideo <- ideo[ideo$name!="", , drop=FALSE]
  ideo$chrom <- as.character(ideo$chrom)
  ideo <- ideo[ideo$chrom %in% seqlevels(range), , drop=FALSE]
  gr <- GenomicRanges::trim(GRanges(ideo))
  suppressWarnings(seqinfo(gr) <- seqinfo(range)[seqlevels(gr)])
  gr
}

