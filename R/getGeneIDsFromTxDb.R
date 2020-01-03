#' get gene ids by genomic location
#' @description retrieve gene ids from txdb object by genomic location.
#' @param gr GRanges object.
#' @param txdb An object of \code{\link[GenomicFeatures:TxDb-class]{TxDb}}. 
#' @return A character vector of gene ids
#' @importFrom AnnotationDbi select
#' @importFrom GenomicFeatures genes
#' @import GenomicRanges
#' @import IRanges
#' @export
#' @examples 
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' gr <- parse2GRanges("chr11:122,830,799-123,116,707")
#' ids <- getGeneIDsFromTxDb(gr, TxDb.Hsapiens.UCSC.hg19.knownGene)

getGeneIDsFromTxDb <- function(gr, txdb){
  stopifnot(is(gr, "GRanges"))
  stopifnot(length(gr)>0)
  stopifnot(is(txdb, "TxDb"))
  if(length(gr)>1){
    warning("The length of gr is greater than 1. Only first genomic location will be used.")
    gr <- gr[1]
  }
  genes <- genes(txdb, columns="gene_id")
  genes <- subsetByOverlaps(genes, gr)
  return(genes$gene_id)
}