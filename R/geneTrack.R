#' track from TxDb
#' @description Generate a track object from TxDb by given gene ids 
#' @param ids Gene IDs. A vector of character. It should be keys in txdb.
#' @param txdb An object of \code{\link[GenomicFeatures:TxDb-class]{TxDb}}.
#' @param type Output type of track, "gene" or "transcript".
#' @param symbols symbol of genes. 
#' @param asList Output a list of tracks or not. Default TRUE.
#' @return An object of \link{track}
#' @importFrom AnnotationDbi select
#' @import GenomicRanges
#' @import IRanges
#' @export
#' @examples 
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' library(org.Hs.eg.db)
#' ids <- c("3312", "390259", "341056", "79827")
#' symbols <- mget(ids, org.Hs.egSYMBOL)
#' geneTrack(ids, TxDb.Hsapiens.UCSC.hg19.knownGene, symbols)
#' 
geneTrack <- function(ids, txdb, symbols, type=c("gene", "transcript"), asList=TRUE){
  stopifnot(is(txdb, "TxDb"))
  stopifnot(is(ids, "character"))
  stopifnot(is(asList, "logical"))
  if(missing(symbols)) symbols <- ids
  if(is.null(names(symbols))){
    names(symbols) <- ids
  }
  symbols <- symbols[ids]
  stopifnot(length(ids)==length(symbols))
  symbols <- sapply(symbols, `[`, i=1)
  type <- match.arg(type)
  getFeatureType <- function(.ele, .cds, .id, .txname, .gene){
    .ele <- disjoin(.ele)
    .cds <- disjoin(.cds)
    .ele <- c(.ele, .cds)
    .ele.disjoin <- disjoin(.ele, with.revmap = TRUE)
    .ele.disjoin$feature <- ifelse(lengths(.ele.disjoin$revmap)==1, "ncRNA", "CDS")
    .ele.disjoin <- sort(.ele.disjoin, 
                         decreasing = as.character(strand(.ele.disjoin)[1])=="-")
    cds <- range(which(.ele.disjoin$feature=="CDS"))
    if(cds[1]>1) .ele.disjoin$feature[seq.int(cds[1]-1)] <- "utr5"
    if(cds[2]<length(.ele.disjoin)) .ele.disjoin$feature[seq(cds[2]+1, length(.ele.disjoin))] <- "utr3"
    .ele.disjoin$id <- .id
    .ele.disjoin$exon <- paste(.id, sapply(.ele.disjoin$revmap, `[`, i=1), sep="_")
    .ele.disjoin$transcript <- .txname
    .ele.disjoin$gene <- .gene
    .ele.disjoin$revmap <- NULL
    .ele.disjoin
  }
  if(type == "gene"){
    suppressMessages({
      exons <- select(txdb, ids, 
                      c("GENEID", "EXONCHROM", "EXONEND", "TXNAME", "EXONSTART", "EXONSTRAND"),
                      "GENEID")
      cds <- select(txdb, unique(exons$TXNAME), 
                    c("TXNAME", "CDSCHROM", "CDSSTRAND", "CDSSTART", "CDSEND"),
                    "TXNAME")
      cds <- merge(unique(exons[, c("GENEID", "TXNAME")]), cds, all.x=TRUE)
    })
    exons <- split(exons, exons$GENEID)
    cds <- split(cds, cds$GENEID)
    exons <- mapply(function(.ele, .cds, .name, .id){
      .ele <- GRanges(seqnames = .ele$EXONCHROM, 
                      ranges = IRanges(.ele$EXONSTART, .ele$EXONEND),
                      strand = .ele$EXONSTRAND)
      .cds <- .cds[!is.na(.cds$CDSCHROM), , drop=FALSE]
      if(nrow(.cds)>0){
        .cds <- GRanges(seqnames = .cds$CDSCHROM,
                        ranges = IRanges(.cds$CDSSTART, .cds$CDSEND),
                        strand = .cds$CDSSTRAND)
        .ele <- getFeatureType(.ele, .cds, .id, "unknown", .name)
      }else{
        .ele$feature <- "ncRNA"
      }
      new("track", dat=.ele, type="gene", 
          name=.name,
          style=new("trackStyle", color="lightblue"))
    }, exons, cds[names(exons)], symbols[names(exons)], names(exons))
    if(!asList){
      exons <- lapply(exons, function(.ele){
        .e <- .ele@dat
        .e$featureID <- .ele@name
        .e
      })
      exons <- unlist(GRangesList(exons))
      exons <- new("track", dat=exons, type="gene",
                   name="genes",
                   style=new("trackStyle", color="black"))
    }
    return(exons)
  }
  ## type == "transcript"
  suppressMessages({
    exons <- select(txdb, ids, 
                  c("GENEID", "EXONCHROM", "EXONSTRAND", "TXNAME", "EXONSTART", "EXONEND"),
                  "GENEID")
  })
  exons <- split(exons, exons$GENEID)
  trs <- mapply(function(exon, .id, .name){
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
      cds <- cds[!is.na(cds$CDSCHROM), , drop=FALSE]
      if(nrow(cds)>0){
        cds <- GRanges(seqnames = cds$CDSCHROM,
                       ranges = IRanges(cds$CDSSTART, cds$CDSEND),
                       strand = cds$CDSSTRAND)
        cds$TXNAME <- .ele$TXNAME[1]
        .ele <- getFeatureType(.ele, cds, .ele$TXNAME[1], .ele$TXNAME[1], .name)
      }else{
        .ele$feature <- "ncRNA"
        .ele$id <- .id
        .ele$transcript <- .ele$TXNAME
        .ele$gene <- .name
      }
      
      new("track", dat=.ele, type="transcript", 
          name=.ele$transcript[1],
          style=new("trackStyle", color="lightblue"))
    })
  }, exons, names(exons), symbols[names(exons)])
  names(trs) <- symbols[names(exons)]
  trs <- unlist(trs)
  
  if(!asList){
    trs <- lapply(trs, function(.ele){
      .e <- .ele@dat
      .e$featureID <- .ele@name
      .e
    })
    trs <- unlist(GRangesList(trs))
    trs <- new("track", dat=trs, type="transcript",
                 name="transcripts",
                 style=new("trackStyle", color="black"))
  }
  trs
}