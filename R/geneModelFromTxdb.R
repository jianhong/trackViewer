#' Prepare gene model from an object of TxDb
#' @description Generate an object of \code{\link{track}} for 
#' \code{\link{viewTracks}} by given parameters.
#' @param txdb An object of \code{\link[GenomicFeatures:TxDb-class]{TxDb}}
#' @param orgDb An object of "OrgDb"
#' @param gr An object of GRanges.
#' @param chrom chromosome name, must be a seqname of txdb
#' @param start start position
#' @param end end position
#' @param strand strand
#' @param txdump output of as.list(txdb), a list of data frames that can be used 
#' to make the db again with no loss of information.
#' @import GenomicRanges
#' @importFrom GenomicFeatures transcripts genes exonsBy cdsBy fiveUTRsByTranscript threeUTRsByTranscript
#' @importFrom txdbmaker makeTxDb
#' @importFrom Seqinfo seqnames
#' @importFrom BiocGenerics strand start end
#' @importFrom AnnotationDbi mapIds columns
#' @return Generate a list of \code{\link{track}} from a TxDb object.
#' @export
#' @seealso See Also as \code{\link{importScore}}, \code{\link{importBam}}, 
#' \code{\link{viewTracks}}
#' @examples 
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' library(org.Hs.eg.db)
#' trs <- geneModelFromTxdb(TxDb.Hsapiens.UCSC.hg19.knownGene,
#'                          org.Hs.eg.db,
#'                          chrom="chr20", 
#'                          start=22560000, 
#'                          end=22565000, 
#'                          strand="-")
geneModelFromTxdb <- function(txdb, orgDb, gr, 
                              chrom, start, end, 
                              strand=c("*", "+", "-"), 
                              txdump=NULL){
    if(missing(txdb))
        stop("txdb is required!")
    stopifnot(is(txdb, "TxDb"))
    if(!missing(orgDb)) stopifnot(is(orgDb, "OrgDb"))
    strand <- match.arg(strand)
    if(missing(gr)){
        if(missing(chrom)||missing(start)||missing(end)){
            stop("chromosome location is required.")
        }
        gr <- GRanges(chrom, IRanges(start, end), strand=strand)
    }
    if(length(gr)>1){
        warning("only the first element of gr is used")
        gr <- gr[1]
    }
    stopifnot(is(gr, "GRanges"))
    chrom <- as.character(seqnames(gr))
    strand <- as.character(strand(gr))
    start <- BiocGenerics::start(gr)
    end <- BiocGenerics::end(gr)
    if(strand=="*"){
        genes <- suppressMessages(transcripts(txdb, columns="exon_id",
                                              filter=list(tx_chrom=chrom)))
    }else{
        genes <- suppressMessages(transcripts(txdb, columns="exon_id", 
                                              filter=list(tx_chrom=chrom, 
                                                          tx_strand=strand)))
    }
    ignore.strand <- strand=="*"
    ol <- findOverlaps(gr, genes, ignore.strand=ignore.strand)
    if(length(ol)>0){
        exons <- getGeneModel(txdb, chrom, start, end, strand)
        if(!missing(orgDb)){
            if("SYMBOL" %in% columns(orgDb)){
              suppressMessages(symbol <- tryCatch(mapIds(x=orgDb, keys=exons$gene, 
                                                         column="SYMBOL", keytype="ENTREZID",
                                                         multiVals="first"), 
                                                  error=function(e) NA))
                if(!is.na(symbol[1]) && length(symbol)==length(exons)){
                    exons$symbol <- symbol
                }
            }
        }
        trs <- split(exons, as.character(exons$transcript))
        idx <- sapply(trs, function(.ele){
            .r <- range(.ele)
            if(end(.r)>=start && start(.r)<=end) return(TRUE)
            return(FALSE)
        })
        trs <- trs[idx]
        trs <- lapply(trs, function(.ele){
            new("track", dat=.ele, type="transcript", 
                name=as.character(.ele$symbol)[1],
                style=new("trackStyle", color="lightblue"))
        })
        return(trs)
    }else{
        stop("No transcripts in the given range.")
    }
}

getGeneModel <- function(txdb, chrom, start, end, strand){
  txs <- transcripts(txdb, columns=c("tx_id", "tx_name", "gene_id"),
                     filter=list(tx_chrom=chrom))
  txs <- subsetByOverlaps(txs, ranges=GRanges(chrom, IRanges(start, end),
                                              strand=strand))
  txs$gene_id <- vapply(txs$gene_id, function(.ele) .ele[1], character(1L))
  exons <- exonsBy(txdb, by='tx')
  cds <- cdsBy(txdb, by='tx')
  utr5 <- fiveUTRsByTranscript(txdb)
  utr3 <- threeUTRsByTranscript(txdb)
  exons <- subsetByOverlaps(exons, ranges=txs)
  cds <- subsetByOverlaps(cds, ranges=txs)
  utr5 <- subsetByOverlaps(utr5, ranges = txs)
  utr3 <- subsetByOverlaps(utr3, ranges = txs)
  gene_model <- lapply(names(exons), function(tx_id){
    tx <- GRanges()
    if(tx_id %in% names(utr5)){
      ele <- utr5[[tx_id]]
      ele$feature <- 'utr5'
      tx <- c(tx, ele)
    }
    if(tx_id %in% names(cds)){
      ele <- cds[[tx_id]]
      ele$feature <- 'CDS'
      tx <- c(tx, ele)
    }
    if(tx_id %in% names(utr3)){
      ele <- utr3[[tx_id]]
      ele$feature <- 'utr3'
      tx <- c(tx, ele)
    }
    if(length(tx)==0){
      tx <- exons[[tx_id]]
      tx$feature <- 'ncRNA'
    }
    tx
  })
  
  
  gene_model_ul <- unlist(GRangesList(gene_model))
  gene_model_ul$tx_id <- rep(names(exons), lengths(gene_model))
  idx <- match(gene_model_ul$tx_id, txs$tx_id)
  gene_model_ul$transcript <- txs$tx_name[idx]
  gene_model_ul$gene <- txs$gene_id[idx]
  gene_model_ul$symbol <- gene_model_ul$transcript
  gene_model_ul$exon_id <- NULL
  gene_model_ul$exon_name <- NULL
  gene_model_ul$cds_id <- NULL
  gene_model_ul$cds_name <- NULL
  return(gene_model_ul)
}