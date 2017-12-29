#' plot gene track
#' @description different from transcript track, gene track is a container.
#' It can contain multiple genes in same chromosome. The input should be
#' a GRanges with strand information. The gene will not show difference between
#' UTR and CDS region. It will be showed as connected exons with TSS arrow.
#' Gene track will be plotted in multiple layers if needed. The algorithm is
#' like this: 1) if there is no overlaps of different genes, one layer.
#' 2) if there are overlaps, first arrage all the un-overlaped genes. for the
#' overlaped genes, try to locate as much as possible gene in upper layer.
#' then alocate the rest. If two layers can not handle, use three or more.
#' Gene track may be a connection track for RNAseq/ChIPseq, etc (low) with lollipops.
#' So, the extension line may be applied to Gene track.
#' @param track a gene track
#' @param gr the range to be plotted. GRanges object 
#' @param txdb if track is not supplied, txdb will be used to extract the genes
#' @param org org database
#' @examples 
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' library(org.Hs.eg.db)
#' gr <- GRanges("chr11", IRanges(122800000, 123000000))
#' plotGeneTrack(gr=gr, txdb=TxDb.Hsapiens.UCSC.hg19.knownGene, org=org.Hs.eg.db)

# plotGeneTrack <- function(track, gr, txdb, org){
#   if(missing(track)){
#     if(missing(txdb) || missing(gr)){
#       stop("txdb and gr is required if missing track")
#     }
#     stopifnot(is(txdb, "TxDb"))
#     stopifnot(is(gr, "GRanges"))
#     stopifnot(is(org, "OrgDb"))
#     ## extract genes
#     exons <- exons(txdb, columns=c("exon_id", "gene_id"))
#     exons <- subsetByOverlaps(exons, gr, ignore.strand=FALSE)
#     exons$gene_id <- sapply(exons$gene_id, `[`, 1)
#     suppressMessages(
#       exons$symbol <- mapIds(org, exons$gene_id, "SYMBOL", 
#                            "ENTREZID", multiVals = "first"))
#     exons <- split(exons, exons$symbol)
#     exons <- lapply(exons, reduce)
#   }else{
#     stopifnot(is(track, "track"))
#     if(track@type!="gene"){
#       stop("track must be a gene track.")
#     }
#     if(missing(gr)){
#       gr <- granges(track@dat)[1]
#     }
#     exons <- () ## change track into a GRangesList named by genenames
#   }
#   ## determin how many tracks should be plotted.
#   genes <- unlist(GRangesList(lapply(exons, range)))
#   ol <- findOverlaps(genes, drop.self=TRUE, drop.redundant=TRUE)
#   if(length(ol)>0){
#     
#   }
# }