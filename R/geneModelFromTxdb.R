geneModelFromTxdb <- function(txdb, orgDb, gr, 
                              chrom, start, end, 
                              strand=c("*", "+", "-"), 
                              txdump=NULL){
    if(missing(txdb)||missing(orgDb))
        stop("txdb and orgDb are required!")
    stopifnot(class(txdb)=="TxDb")
    stopifnot(class(orgDb)=="OrgDb")
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
    stopifnot(class(gr)=="GRanges")
    chrom <- as.character(seqnames(gr))
    strand <- as.character(strand(gr))
    start <- gr@ranges@start
    end <- start + gr@ranges@width - 1
    if(strand=="*"){
        genes <- transcripts(txdb, columns="exon_id", 
                       vals=list(tx_chrom=chrom))
    }else{
        genes <- transcripts(txdb, columns="exon_id", 
                       vals=list(tx_chrom=chrom, tx_strand=strand))
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
            symbol <- tryCatch(mapIds(x=orgDb, keys=exons$gene, 
                                      column="SYMBOL", keytype="ENTREZID",
                                      multiVals="first"), error=function(e) NA)
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
            new("track", dat=.ele, type="gene", 
                name=as.character(.ele$symbol)[1],
                style=new("trackStyle", color="lightblue"))
        })
        return(trs)
    }else{
        stop("No transcripts in the given range.")
    }
}
