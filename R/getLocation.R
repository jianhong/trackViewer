#' get genomic location by gene symbol
#' @description given a gene name, get the genomic coordinates.
#' @param symbol Gene symbol
#' @param txdb txdb will be used to extract the genes
#' @param org org package name
#' @import GenomicRanges
#' @importFrom AnnotationDbi mget
#' @importFrom GenomeInfoDb "seqlevelsStyle<-" seqlevels "seqlevels<-" seqinfo "seqinfo<-"
#' @importFrom GenomicFeatures genes
#' @import IRanges
#' @export
#' @examples 
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' library(org.Hs.eg.db)
#' getLocation("HSPA8", TxDb.Hsapiens.UCSC.hg19.knownGene, "org.Hs.eg.db")

getLocation <- function(symbol, txdb, org){
  stopifnot(is(symbol, "character"))
  stopifnot(is(txdb, "TxDb"))
  stopifnot(is.character(org))
  stopifnot(length(org)==1)
  if(!grepl("^org.*.db", org)){
    stop("Org package name is required.")
  }
  dbPrefix <- sub(".db", "", org)
  eg <- mget(symbol, get(paste0(dbPrefix, "SYMBOL2EG")), ifnotfound = NA)
  eg <- lapply(eg, `[`, i=1)
  eg <- unlist(eg, use.names = TRUE)
  eg <- eg[!is.na(eg)]
  if(length(eg)<1){
    return(NULL)
  }
  genes <- genes(txdb)
  gr <- genes[eg]
  seqlevelsStyle(gr) <- "UCSC"
  seqlevels(gr) <- seqlevels(gr)[seqlevels(gr) %in% as.character(seqnames(gr))]
  if(length(seqinfo(gr))>0){
    seqinfo(gr) <- seqinfo(gr)[seqlevels(gr)]
  }
  gr
}