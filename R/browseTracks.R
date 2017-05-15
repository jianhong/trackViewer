#' interactive trackViewer
#'
#' @description plot tracks via d3.js
#'
#' @import htmlwidgets
#' 
browseTracks <- function(trackList, 
                         gr=GRanges(),
                         ignore.strand=TRUE,
                         width=NULL, height=NULL, 
                         ...){
    if(missing(trackList)){
        stop("trackList is required.")
    }
    if(class(trackList)!="trackList" && 
       !((is.list(trackList) && all(sapply(trackList, class)=="track")))){
        stop("trackList must be an object of \"trackList\"
             (See ?trackList) or a list of track")
    }
    if(any(is.na(names(trackList)))){
        stop("trackList must have names")
    }
    if(any(duplicated(names(trackList)))){
        stop("trackList must have unqiue names")
    }
    if(missing(gr)){
        stop("gr is required.")
    }else{
        if(!is(gr, "GRanges")){
            stop("gr must be an object of GRanges.")
        }
        if(length(gr)!=1){
            stop("the length of gr must be 1.")
        }
    }
    chromosome <- as.character(GenomicRanges::seqnames(gr))[1]
    start <- GenomicRanges::start(gr)[1]
    end <- GenomicRanges::end(gr)[1]
    strand <- as.character(GenomicRanges::strand(gr))[1]
    if(ignore.strand) strand <- "*"
    if(end < start) stop("end should be greater than start.")
    
    trackList <- trackViewer:::filterTracks(trackList, chromosome, 
                              start, end, strand)
    
    viewerStyle=trackViewerStyle()
    ## trackList split into two parts: gene and data
    irl <- IRangesList(ranges(gr)[1])
    names(irl) <- chromosome
    object2list <- function(obj){
        if(!isS4(obj)) return(obj)
        s <- slotNames(obj)
        x <- lapply(s, function(.ele){
            .ele <- slot(obj, .ele)
            if(isS4(.ele)) object2list(.ele) else .ele
        })
        names(x) <- s
        x
    }
    ZERO <- numeric(0)
    getData <- function(.dat){
        if(length(.dat)>0){
            .dat <- .dat[c(1, seq_along(.dat))]
            start(.dat[1]) <- end(.dat[1]) <- end+1
            .dat[1]$score <- 0
            strand(.dat) <- "*"
            co <- coverage(.dat, weight=.dat$score)
            as.numeric(co[irl][[1]])
        }else{
            ZERO
        }
    }
    col2Hex <- function(x){
        rgb <- col2rgb(x,)
        rgb(rgb[1, ], rgb[2, ], rgb[3, ], maxColorValue = 255)
    }
    trackdat <- function(.ele){
        dat <- .ele$dat
        dat2 <- .ele$dat2
        if(.ele$type=="data"){
            dat <- getData(.ele$dat)
            dat2 <- getData(.ele$dat2)
            ylim <- range(c(dat, -1*dat2))
        }else{
            dat <- as.list(as.data.frame(.ele$dat))
            dat$seqnames <- chromosome
            dat$strand <- as.character(dat$strand)[1]
            dat2 <- ZERO
            ylim <- c(0, 1)
        }
        style <- object2list(.ele$style)
        style$color <- col2Hex(style$color)
        list(dat=dat,
             dat2=dat2,
             style=style,
             ylim=ylim)
    }
    x <- list(
        tracklist = lapply(trackList, trackdat),
        type = lapply(trackList, function(.ele) .ele$type),
        name = names(trackList),
        chromosome = chromosome,
        start = start,
        end = end,
        strand =strand,
        height = lapply(trackList, function(.ele) .ele$style@height)
    )
    
    htmlwidgets::createWidget(
        name = 'browseTracks',
        x = x,
        width = width,
        height = height,
        package = getPackageName()
    )
}
