viewTracks <- function(trackList, chromosome, start, end, strand, gr=GRanges(),
                       ignore.strand=TRUE,
                       viewerStyle=trackViewerStyle(), autoOptimizeStyle=FALSE,
                       newpage=TRUE, operator=NULL){
    if(!is.null(operator)){
        if(!operator %in% c("+", "-", "*", "/", "^", "%%")){
            stop('operator must be one of "+", "-", "*", "/", "^", "%%"')
        }
    }
    if(missing(trackList)){
        stop("trackList is required.")
    }
    if(class(trackList)!="trackList" && 
           !((is.list(trackList) && all(sapply(trackList, class)=="track")))){
        stop("trackList must be an object of \"trackList\"
             (See ?trackList) or a list of track")
    }
    filterTracksFlag <- TRUE
    if(missing(gr)){
        if(missing(chromosome) || missing(start) || missing(end))
            stop("Please input the coordinate.")
        if(missing(strand) || !strand %in% c("+", "-"))
            strand <- "*"
        wavyLine <- FALSE
    }else{
        if(!is(gr, "GRanges")){
            stop("gr must be an object of GRanges.")
        }
        if(length(gr)==0){
            stop("the length of gr must greater than 0.")
        }
        if(length(gr$filterTracks)>0){
            filterTracksFlag <- as.logical(gr$filterTracks[1])
        }
        grs <- range(gr)
        chromosome <- as.character(GenomicRanges::seqnames(grs))[1]
        start <- GenomicRanges::start(grs)[1]
        end <- GenomicRanges::end(grs)[1]
        strand <- as.character(GenomicRanges::strand(grs))[1]
        if(ignore.strand) strand <- "*"
        if(length(gr)>1){ ## more than one region
            if(length(unique(as.character(seqnames(gr))))>1){
                stop("gr has multiple seqnames.")
            }
            if(length(gr$percentage)==0){
                gr$percentage <- 1/length(gr)
            }else{
                gr$percentage <- gr$percentage/sum(gr$percentage, na.rm = TRUE)
                gr$percentage[is.na(gr$percentage)] <- 0
            }
            gr$wavyLine <- TRUE
            gr$wavyLine[length(gr)] <- FALSE
            gr$filterTracks <- FALSE
            if(newpage) grid.newpage()
            
            if(autoOptimizeStyle){
                opt <- optimizeStyle(trackList, viewerStyle)
                trackList <- opt$tracks
                viewerStyle <- opt$style
            }
            trackList <- filterTracks(trackList, chromosome, start, end, strand)
            for(i in seq_along(gr)){
                current.viewerStyle <- viewerStyle
                current.trackList <- trackList
                hgap <- convertWidth(unit(0.25, "lines"), 
                                     unitTo = "npc", 
                                     valueOnly = TRUE)
                if(i==1){
                    setTrackViewerStyleParam(current.viewerStyle, 
                                             "margin", 
                                             c(current.viewerStyle@margin[1:3], 
                                               hgap))
                    for(j in seq_along(current.trackList)){
                        tLjStyle <- current.trackList[[j]]@style
                        if(tLjStyle@yaxis@draw && !tLjStyle@yaxis@main){
                            setTrackYaxisParam(current.trackList[[j]],
                                               "draw",
                                               FALSE)
                        }
                        if(!tLjStyle@ylabpos %in% 
                           c("upstream", "downstream") &&
                           grepl("right", tLjStyle@ylabpos)){
                            names(current.trackList)[[j]] <- ""
                        }
                    }
                }else{
                    if(i==length(gr)){
                        setTrackViewerStyleParam(
                            current.viewerStyle, 
                            "margin", 
                            c(current.viewerStyle@margin[1],
                              hgap,
                              current.viewerStyle@margin[3:4])
                        )
                        for(j in seq_along(current.trackList)){
                            tLjStyle <- current.trackList[[j]]@style
                            if(tLjStyle@yaxis@draw && tLjStyle@yaxis@main){
                                setTrackYaxisParam(current.trackList[[j]],
                                                   "draw",
                                                   FALSE)
                            }
                            if(!tLjStyle@ylabpos %in% 
                               c("upstream", "downstream") &&
                               grepl("left", tLjStyle@ylabpos)){
                                names(current.trackList)[[j]] <- ""
                            }
                        }
                    }else{
                        setTrackViewerStyleParam(
                            current.viewerStyle, 
                            "margin", 
                            c(current.viewerStyle@margin[1],
                              hgap,
                              current.viewerStyle@margin[3],
                              hgap)
                        )
                        for(j in seq_along(current.trackList)){
                            setTrackYaxisParam(current.trackList[[j]],
                                               "draw",
                                               FALSE)
                            if(!current.trackList[[j]]@style@ylabpos %in% 
                               c("upstream", "downstream")){
                                names(current.trackList)[[j]] <- ""
                            }
                        }
                    }
                }
                pushViewport(vp=viewport(
                    x = unit(ifelse(i==1, 
                                    0, 
                                    sum(gr$percentage[seq_len(i-1)]))+
                                 gr$percentage[i]/2-hgap/2, 
                             "npc"),
                    width = unit(gr$percentage[i]-hgap, "npc"),
                    name = paste0("vp.track.v.", i)
                ))
                vp <- viewTracks(trackList = current.trackList,
                                 gr = gr[i],
                                 ignore.strand = ignore.strand,
                                 viewerStyle = current.viewerStyle,
                                 autoOptimizeStyle = FALSE,
                                 newpage = FALSE, operator = operator)
                vp$name <- paste0("panel.", i)
                pushViewport(vp)
                upViewport()
                upViewport()
            }
            return(invisible(current.vpTree()))
        }
        wavyLine <- ifelse(length(gr$wavyLine)>0, as.logical(gr$wavyLine[1]),
                           FALSE)
    }
    if(ignore.strand) strand <- "*"
    if(end < start) stop("end should be greater than start.")
    
    if(class(viewerStyle)!="trackViewerStyle"){
        stop("viewerStyle must be an object of 'trackViewerStyle'.")
    }
    if(!is.logical(newpage)) stop("newpage should be a logical vector.")
    
    if(filterTracksFlag) {
        trackList <- filterTracks(trackList, chromosome, start, end, strand)
    }
    
    if(!is.null(operator)){
        ##change dat as operator(dat, dat2)
        ##dat2 no change
        trackList <- lapply(trackList, function(.ele){
            if(.ele@type=="data" && 
                   length(.ele@dat2)>0 &&
                   length(.ele@dat)>0){
                .ele@dat <- GRoperator(.ele@dat, 
                                       .ele@dat2, 
                                       col="score", 
                                       operator=operator)
                if(operator!="+") .ele@dat2 <- GRanges()
            }
            .ele
        })
    }
    
    if(newpage) grid.newpage()
    
    if(autoOptimizeStyle){
        opt <- optimizeStyle(trackList, viewerStyle)
        trackList <- opt$tracks
        viewerStyle <- opt$style
    }
    
    margin <- viewerStyle@margin
    
    ##xscale, yscale
    xscale <- c(start, end)
    
    yscales <- getYlim(trackList, operator)
    yHeights <- getYheight(trackList)
    
    pushViewport(viewport(x=margin[2], y=margin[1], width=1-margin[2]-margin[4],
                          height=1-margin[1]-margin[3], just=c(0,0)))
    if(viewerStyle@xaxis){##draw x axis in bottom
        drawXaxis(xscale, viewerStyle)
    }
    popViewport()
    
    pushViewport(viewport(x=0, y=margin[1], width=1,
                          height=1-margin[1]-margin[3], just=c(0,0)))
    total <- length(trackList)
#    if(interactive()) pb <- txtProgressBar(min=0, max=total, style=3)
    ht <- 0
    xy <- vector(mode="list", length=total)
    yHeightBottom <- yHeightTop <- rep(0.01, total)
    for(i in 1:total){
        xy[[i]] <- plotTrack(names(trackList)[i], trackList[[i]], 
                             viewerStyle, ht,
                             yscales[[i]], yHeights[i], xscale,
                             chromosome, strand, operator, wavyLine)
        ht <- ht + yHeights[i]
        if(length(trackList[[i]]@style@marginBottom)>0){
            yHeightBottom[i] <- trackList[[i]]@style@marginBottom
        }
        if(length(trackList[[i]]@style@marginTop)>0){
            yHeightTop[i] <- trackList[[i]]@style@marginTop
        }
#        if(interactive()) setTxtProgressBar(pb, i)
    }
#    if(interactive()) close(pb)
    popViewport()
    if(viewerStyle@flip) xscale <- rev(xscale)
    options(LastTrackViewer=list(vp=viewport(x=margin[2], y=margin[1], 
                                             height=1 - margin[1]- margin[3], 
                                             width=1 -margin[2] - margin[4],
                                             just=c(0,0), xscale=xscale, yscale=c(0,1)),
                                 xy=xy,
                                 yHeights=yHeights,
                                 yscales=yscales,
                                 xscale=xscale,
                                 yHeightBottom=yHeightBottom,
                                 yHeightTop=yHeightTop,
                                 viewerStyle=viewerStyle))
    return(invisible(viewport(x=margin[2], y=margin[1], 
                              height=1 - margin[1]- margin[3], 
                              width=1 -margin[2] - margin[4],
                              just=c(0,0), xscale=xscale, yscale=c(0,1))))
}