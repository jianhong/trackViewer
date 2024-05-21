library(FLAMINGOrLite)
res = flamingo_main(hic_data='4DNFI1UEG1HD.hic',
                    file_format='hic',
                    domain_res=1e6,
                    frag_res=5e3,
                    chr_name='chr21',
                    normalization='KR',
                    nThread=20)
library(GenomicRanges)
gr <- with(res, GRanges(seqnames = chr, ranges=IRanges(start+1, end),
                        x=x, y=y, z=z))
range <- GRanges('chr21:31375000-32985000')
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
p <- subsetByOverlaps(gr, range)
feature.gr <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
feature.gr <- subsetByOverlaps(feature.gr, range)
symbols <- mget(feature.gr$gene_id, org.Hs.egSYMBOL, ifnotfound=NA)
feature.gr$label[lengths(symbols)==1] <- unlist(symbols[lengths(symbols)==1])
feature.gr$col <- sample(1:7, length(feature.gr), replace=TRUE)
feature.gr$type <- "gene"
saveRDS(feature.gr, '4DNFI1UEG1HD.feature.gr.rds')
view3dStructure(p, k=3, feature.gr=feature.gr, length.arrow=unit(0.000006, 'native'))