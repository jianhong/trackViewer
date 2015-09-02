importData <- 
    function(files, 
             format=NA,
             ranges=GRanges()){
    stopifnot(class(ranges)=="GRanges")
    ranges <- reduce(ranges)
    formats <- c("BAM", "BED", "bedGraph", "WIG", "BigWig")
    if(is.na(format[1])){
        format <- sapply(files, file_ext) ## need tools package
    }
    format <- rep(format, length(files))[1:length(files)]
    format <- formats[match(tolower(format), tolower(formats))]
    if(any(is.na(format)))
        stop("Can not determine the file format.",
             "Format must be BAM, BED, bedGraph, WIG, or BigWig")
    covBAM <- NULL
    covBW <- NULL
    covBED <- NULL
    covBedGraph <- NULL
    covWIG <- NULL
    bams <- files[format=="BAM"]
    BigWigs <- files[format=="BigWig"]
    BEDs <- files[format=="BED"]
    bedGraphs <- files[format=="bedGraph"]
    WIGs <- files[format=="WIG"]
    if(length(bams)>0){
        ## import bams, need Rsamtools, need GenomicAlignments
        param <- ScanBamParam(which=ranges)
        bam <- lapply(bams, readGAlignments, param=param)
        names(bam) <- bams
        covBAM <- sapply(bam, coverage)
    }
    if(length(BigWigs)>0){
        ## import by rtracklayer
        ### if format == "BigWig"
        covBW <- sapply(BigWigs, import, 
                       format="BigWig", 
                       which=ranges, 
                       as="RleList")
    }
    if(length(BEDs)>0){
        covBED <- sapply(BEDs, importScore, 
                         format="BED",
                         ranges=ranges,
                         ignore.strand=TRUE)
        covBED <- sapply(covBED, function(.ele){
            coverage(.ele@dat, weight=score(.ele@dat))
        })
    }
    if(length(bedGraphs)>0){
        covBedGraph <- sapply(bedGraphs, importScore, 
                              format="bedGraph",
                              ranges=ranges,
                              ignore.strand=TRUE)
        covBedGraph <- sapply(covBedGraph, function(.ele){
            coverage(.ele@dat, weight=score(.ele@dat))
        })
    }
    if(length(WIGs)>0){
        ranges <- reduce(ranges, ignore.strand=TRUE)
        covWIG <- sapply(WIGs, importScore, 
                              format="WIG",
                              ranges=ranges,
                              ignore.strand=TRUE)
        covWIG <- sapply(covWIG, function(.ele){
            .tracks <- lapply(seq_len(length(ranges)), function(id){
                .gr <- ranges[id]
                .out <- parseWIG(trackScore=.ele,
                                 chrom=as.character(seqnames(.gr)),
                                 from=start(.gr),
                                 to=end(.gr))
                .out@dat
            })
            .tracks <- unlist(GRangesList(.tracks))
            coverage(.tracks, weight=score(.tracks))
        })
    }
    covs <- c(covBAM, covBW, covBED, covBedGraph, covWIG)
    covs[files]
}