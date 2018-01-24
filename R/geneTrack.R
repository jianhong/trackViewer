#' track from TxDb
#' @description Generate a track object from TxDb by given gene ids 
#' @param ids Gene IDs. A vector of character. It should be keys in txdb.
#' @param txdb An object of \code{\link[GenomicFeatures:TxDb-class]{TxDb}}
#' @param type Output type of track, "gene" or "transcript".
#' @return An object of \link{track}
#' @importFrom AnnotationDbi select
#' @import GenomicRanges
#' @import IRanges
#' @export
#' @examples 
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' geneTrack(c("3312", "3313"), TxDb.Hsapiens.UCSC.hg19.knownGene)
#' 
geneTrack <- function(ids, txdb, type=c("gene", "transcript")){
  stopifnot(is(txdb, "TxDb"))
  stopifnot(is(ids, "character"))
  type <- match.arg(type)
  if(type == "gene"){
    suppressMessages({
      exons <- select(txdb, ids, 
                      c("GENEID", "EXONCHROM", "EXONEND", "EXONSTART", "EXONSTRAND"),
                      "GENEID")
    })
    exons <- split(exons, exons$GENEID)
    exons <- mapply(function(.ele, .name){
      .ele <- GRanges(seqnames = .ele$EXONCHROM, 
                      ranges = IRanges(.ele$EXONSTART, .ele$EXONEND),
                      strand = .ele$EXONSTRAND, 
                      feature = rep("exon", nrow(.ele)))
      new("track", dat=.ele, type="gene", 
          name=.name,
          style=new("trackStyle", color="lightblue"))
    }, exons, names(exons))
    return(exons)
  }
  ## type == "transcript"
  suppressMessages({
    exons <- select(txdb, ids, 
                  c("GENEID", "EXONCHROM", "EXONSTRAND", "TXNAME", "EXONSTART", "EXONEND"),
                  "GENEID")
  })
  exons <- split(exons, exons$GENEID)
  trs <- mapply(function(exon, .name){
    exon <- GRanges(seqnames = exon$EXONCHROM,
                    ranges = IRanges(exon$EXONSTART, exon$EXONEND),
                    strand = exon$EXONSTRAND, TXNAME=exon$TXNAME)
    txs <- split(exon, exon$TXNAME)
    txs <- lapply(txs, function(.ele){
      suppressMessages({
        cds <- select(txdb, .ele$TXNAME[1], 
                      c("TXNAME", "CDSCHROM", "CDSSTRAND", "CDSSTART", "CDSEND"),
                      "TXNAME")
      })
      cds <- GRanges(seqnames = cds$CDSCHROM,
                     ranges = IRanges(cds$CDSSTART, cds$CDSEND),
                     strand = cds$CDSSTRAND)
      cds$TXNAME <- .ele$TXNAME[1]
      .ele <- c(.ele, cds)
      .ele.disjoin <- disjoin(.ele, with.revmap = TRUE)
      .ele.disjoin$feature <- ifelse(lengths(.ele.disjoin$revmap)==1, "ncRNA", "CDS")
      .ele.disjoin <- sort(.ele.disjoin, 
                           decreasing = as.character(strand(.ele.disjoin)[1])=="-")
      cds <- range(which(.ele.disjoin$feature=="CDS"))
      if(cds[1]>1) .ele.disjoin$feature[seq.int(cds[1]-1)] <- "utr5"
      if(cds[2]<length(.ele.disjoin)) .ele.disjoin$feature[seq(cds[2]+1, length(.ele.disjoin))] <- "utr3"
      .ele.disjoin$id <- "unknown"
      .ele.disjoin$exon <- paste(.ele$TXNAME[1], sapply(.ele.disjoin$revmap, `[`, i=1), sep="_")
      .ele.disjoin$transcript <- .ele$TXNAME[1]
      .ele.disjoin$gene <- .name
      .ele.disjoin$revmap <- NULL
      new("track", dat=.ele.disjoin, type="transcript", 
          name=.ele$TXNAME[1],
          style=new("trackStyle", color="lightblue"))
    })
  }, exons, names(exons))
  unlist(trs)
}