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
#' @importFrom Gviz GeneRegionTrack
#' @import GenomicRanges
#' @import GenomicFeatures
#' @importFrom GenomeInfoDb seqnames
#' @importFrom BiocGenerics strand
#' @importFrom AnnotationDbi mapIds columns
#' @return An object of \code{\link{track}}
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
    if(missing(txdb)||missing(orgDb))
        stop("txdb and orgDb are required!")
    stopifnot(is(txdb, "TxDb"))
    stopifnot(is(orgDb, "OrgDb"))
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
    start <- gr@ranges@start
    end <- start + gr@ranges@width - 1
    if(strand=="*"){
        genes <- transcripts(txdb, columns="exon_id", 
                       filter=list(tx_chrom=chrom))
    }else{
        genes <- transcripts(txdb, columns="exon_id", 
                       filter=list(tx_chrom=chrom, tx_strand=strand))
    }
    ignore.strand <- strand=="*"
    ol <- findOverlaps(gr, genes, ignore.strand=ignore.strand)
    if(length(ol)>0){
        exon <- genes[subjectHits(ol)]
        r <- range(exon)
        ##subset txdb
        if(is.null(txdump)) txdump <- as.list(txdb)
        txdump$chrominfo <- txdump$chrominfo[txdump$chrominfo$chrom==chrom, , drop=FALSE]
        txdump$transcripts <- txdump$transcripts[txdump$transcripts$tx_chrom==chrom &
                                                 txdump$transcripts$tx_start<end(r) &
                                                 txdump$transcripts$tx_end>start(r), , drop=FALSE]
        txdump$genes <- txdump$genes[txdump$genes$tx_id %in% txdump$transcripts$tx_id, , drop=FALSE]
        txdump$splicings <- txdump$splicings[txdump$splicings$tx_id %in% txdump$transcripts$tx_id, , drop=FALSE]
        txdb <- do.call(makeTxDb, txdump)
        exons <- GeneRegionTrack(txdb, chromosome=chrom,
                                 start=start(r), end=end(r), strand=strand)##time comsuming
        exons <- exons@range
        if(all(exons$symbol==exons$transcript) && "SYMBOL" %in% columns(orgDb)){
          suppressMessages(symbol <- tryCatch(mapIds(x=orgDb, keys=exons$gene, 
                                                     column="SYMBOL", keytype="ENTREZID",
                                                     multiVals="first"), 
                                              error=function(e) NA))
            if(!is.na(symbol[1]) && length(symbol)==length(exons)){
                exons$symbol <- symbol
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
