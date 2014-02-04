geneModelFromTxdb <- function(txdb, chrom, start, end, strand=c("*", "+", "-"), txdump=NULL){
    if(missing(txdb)||missing(chrom)||missing(start)||missing(end))
        stop("All inputs are required!")
    if(class(txdb)!="TranscriptDb")
        stop("txdb must be an object of TranscriptDb")
    strand <- match.arg(strand)
#     activeSeq <- isActiveSeq(txdb)
#     namesActiveSeq <- names(activeSeq)
#     activeSeq <- namesActiveSeq==chrom
#     names(activeSeq) <- namesActiveSeq
#     isActiveSeq(txdb) <- activeSeq   
#     
#     genes <- transcripts(txdb, columns=c("exon_id", "tx_name"))
#     exons <- exons(txdb, columns=c("exon_id", "tx_name"), 
#                    vals=list(exon_id=unlist(genes$exon_id)))
#     gr <- GRanges(chrom, IRanges(start, end), strand=strand)
#     ignore.strand <- strand=="*"
#     ol <- findOverlaps(gr, exons, ignore.strand=ignore.strand)
#     if(length(ol)>0){
#         exon <- exons[subjectHits(ol)]
#     }
#     ##get transcript UTRs
#     ### Note that, for each transcript, exons that have a 5' UTR are all the exons
#     ### before the first exon with a CDS, including the first exon with a CDS (even
#     ### though this one might actually have a 0-width 5' UTR but we will take care
#     ### of this later
#     cds <- cds(txdb, columns=c("exon_id", "tx_name"))
#     r <- range(exon)
#     cds <- cds[start(cds)<end(r) & end(cds)>start(r)]
    
    if(strand=="*"){
        genes <- transcripts(txdb, columns="exon_id", 
                       vals=list(tx_chrom=chrom))
    }else{
        genes <- transcripts(txdb, columns="exon_id", 
                       vals=list(tx_chrom=chrom, tx_strand=strand))
    }
    
    gr <- GRanges(chrom, IRanges(start, end), strand=strand)
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
        txdb <- do.call(makeTranscriptDb, txdump)
        exons <- GeneRegionTrack(txdb, chromosome=chrom,
                                 start=start(r), end=end(r), strand=strand)##time comsuming
        exons <- exons@range
        trs <- split(exons, as.character(exons$transcript))
        idx <- sapply(trs, function(.ele){
            .r <- range(.ele)
            if(end(.r)>=start && start(.r)<=end) return(TRUE)
            return(FALSE)
        })
        trs <- trs[idx]
        trs <- lapply(trs, function(.ele){
            new("track", dat=.ele, type="gene", 
                style=new("trackStyle", color="lightblue"))
        })
        return(trs)
    }else{
        stop("No transcripts in the given range.")
    }
}