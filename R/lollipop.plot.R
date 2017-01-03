lolliplot <- function(SNP.gr, features=NULL, ranges=NULL,
                      type=c("circle", "pie", "pin", 
                             "pie.stack"),
                      newpage=TRUE, ylab=TRUE, yaxis=TRUE,
                      xaxis=TRUE, legend=NULL, cex=1, 
                      dashline.col="gray80", 
                      jitter=c("node", "label"), ...){
    stopifnot(inherits(SNP.gr, c("GRanges", "GRangesList", "list")))
    stopifnot(inherits(features, c("GRanges", "GRangesList", "list")))
    jitter <- match.arg(jitter)
    if(type!="circle"&&jitter=="label"){
      jitter <- "node"
      warning("if jitter set to label, type must be cirle.")
      message("jitter is set to node.")
    }
    SNP.gr.name <- deparse(substitute(SNP.gr))
    if(class(SNP.gr)=="GRanges"){
        SNP.gr <- GRangesList(SNP.gr)
        names(SNP.gr) <- SNP.gr.name
    }
    len <- length(SNP.gr)
    for(i in 1:len){
        stopifnot(class(SNP.gr[[i]])=="GRanges")
    }
    
    TYPES <- c("circle", "pie", "pin", "pie.stack")
    if(any(!type %in% TYPES)){
        stop("Error in match argument: ",
             paste0("'type' should be one of '",  
                    paste(TYPES, collapse="', '"), "'."))
    }
    types <- rep(type, length=len)[1:len]
    rm(type)
    ############### handle legend ####################
    ## set the legend as a list, 
    ## if all the legend for different tracks is same
    ## set draw legend for last track later
    if(length(legend)>0){
        if(!is.list(legend)){
            tmp <- legend
            legend <- vector(mode = "list", length = len)
            legend[[len]] <- tmp
            rm(tmp)
        }else{
            if(length(legend)==1){
                tmp <- legend[[1]]
                legend <- vector(mode = "list", length = len)
                legend[[len]] <- tmp
                rm(tmp)
            }else{
                if("labels" %in% names(legend)){
                    tmp <- legend
                    legend <- vector(mode = "list", length = len)
                    legend[[len]] <- tmp
                    rm(tmp)
                }else{
                    if(length(legend)<len){
                        length(legend) <- len
                    }
                }
            }
        }
    }
    
    features.name <- deparse(substitute(features))
    
    ################ handle ranges #####################
    ## if !missing(ranges) set ranges as feature ranges
    if(length(ranges)>0){
        stopifnot(class(ranges)=="GRanges")
        ranges <- rep(ranges, length(SNP.gr))[1:length(SNP.gr)]
        stopifnot(length(ranges)==length(SNP.gr))
    }else{
        if(class(features)=="GRanges"){
            ranges <- range(features)[rep(1, len)]
        }else{
            if(length(features)!=len){
                stop("if both SNP.gr and features is GRangesList,",
                     " the lengthes of them should be identical.")
            }
            ranges <- unlist(GRangesList(lapply(features, range)))
        }
    }
    if(class(ranges)=="GRanges"){
        ##cut all SNP.gr by the range
        for(i in len){
            range <- ranges[i]
            stopifnot(all(width(SNP.gr[[i]])==1))
            ol <- findOverlaps(SNP.gr[[i]], range)
            SNP.gr[[i]] <- SNP.gr[[i]][queryHits(ol)]
        }
    }
    
    ################## plot ############################
    ## total height == 1
    height <- 1/len
    height0 <- 0
    if(newpage) grid.newpage()
    for(i in 1:len){
        type <- match.arg(types[i], TYPES)
        if(type=="pin"){ ## read the pin shape file
            pinpath <- system.file("extdata", "map-pin-red.xml", package="trackViewer")
            pin <- readPicture(pinpath)
        }else{
            pin <- NULL
        }
        ## Here we don't know the real height of each tracks
        vp <- viewport(x=.5, y=height0 + height*0.5, width=1, height=height)
        pushViewport(vp)
        LINEW <- as.numeric(convertX(unit(1, "line"), "npc"))
        LINEH <- as.numeric(convertY(unit(1, "line"), "npc"))
        ## GAP the gaps between any elements
        GAP <- .2 * LINEH
        ratio.yx <- 1/as.numeric(convertX(unit(1, "snpc"), "npc"))
        
        ## prepare the feature
        if(inherits(features, c("GRangesList", "list"))){
            feature <- features[[i]]
            stopifnot(class(feature)=="GRanges")
        }else{
            feature <- features
        }
        ## convert height to npc number
        feature$height <- convertHeight2NPCnum(feature$height)
        ## multiple transcripts in one gene could be separated by featureLayerID
        if(length(feature$featureLayerID)!=length(feature)){
            feature$featureLayerID <- rep("1", length(feature))
        }
        feature <- feature[end(feature)>=start(ranges[i]) & 
                               start(feature)<=end(ranges[i])]
        feature$featureLayerID <- as.character(feature$featureLayerID)
        start(feature)[start(feature)<start(ranges[i])] <- start(ranges[i])
        end(feature)[end(feature)>end(ranges[i])] <- end(ranges[i])
        feature.splited <- split(feature, feature$featureLayerID)
        
        ## bottomblank, the transcripts legend height
        bottomblank <- 4
        if(length(names(feature))>0){ ## features legend
            feature.s <- feature[!duplicated(names(feature))]
            ncol <- getColNum(names(feature.s))
            bottomblank <- max(ceiling(length(names(feature.s)) / ncol), 4)
            pushViewport(viewport(x=.5, y=bottomblank*LINEH/2, 
                                  width=1,
                                  height=bottomblank*LINEH,
                                  xscale=c(start(ranges[i]), end(ranges[i]))))
            color <- if(length(unlist(feature.s$color))==length(feature.s)) 
                unlist(feature.s$color) else "black"
            fill <- if(length(unlist(feature.s$fill))==length(feature.s)) 
                unlist(feature.s$fill) else "black"
            pch <- if(length(unlist(feature.s$pch))==length(feature.s)) 
                unlist(feature.s$pch) else 22
            grid.legend(label=names(feature.s), ncol=ncol,
                        byrow=TRUE, vgap=unit(.2, "lines"),
                        pch=pch,
                        gp=gpar(col=color, fill=fill))
            popViewport()
        }else{
          if(length(xaxis)>1 || as.logical(xaxis[1])){
            bottomblank <- 2
          }else{
            bottomblank <- 0
          }
        }
        
        SNPs <- SNP.gr[[i]]
        strand(SNPs) <- "*"
        SNPs <- sort(SNPs)
        ## get the max score and scoreType
        scoreMax0 <- scoreMax <- 
            if(length(SNPs$score)>0) ceiling(max(c(SNPs$score, 1), na.rm=TRUE)) else 1
        if(type=="pie.stack") scoreMax <- length(unique(SNPs$stack.factor))
        if(!type %in% c("pie", "pie.stack")){
            if(length(yaxis)>1 && is.numeric(yaxis)){
                if(length(names(yaxis))!=length(yaxis)){
                    names(yaxis) <- yaxis
                }
                scoreMax0 <- max(yaxis, scoreMax0)
            }
            if(scoreMax>10) {
                SNPs$score <- 10*SNPs$score/scoreMax
                scoreMax <- 10*scoreMax0/scoreMax
            }else{
                scoreMax <- scoreMax0
            }
            scoreType <- 
                if(length(SNPs$score)>0) all(floor(SNPs$score)==SNPs$score) else FALSE
        }else{
            scoreType <- FALSE
        }
        
        ## if the type is caterpillar, there are lollipop in both sides
        ## plot the bottom lollipops first. And push a new viewport
        
        IsCaterpillar <- length(SNPs$SNPsideID) > 0
        if(IsCaterpillar){
            if(any(is.na(SNPs$SNPsideID)) || 
               !all(SNPs$SNPsideID %in% c('top', 'bottom'))){
                warning("Not all SNPsideID is top or bottom")
                IsCaterpillar <- FALSE
            }
        }
        
        if(IsCaterpillar){
            SNPs.top <- SNPs[SNPs$SNPsideID=='top']
            SNPs.bottom <- SNPs[SNPs$SNPsideID=='bottom']
        }else{
            SNPs.top <- SNPs
            SNPs.bottom <- GRanges()
        }
        if(length(SNPs.bottom)<1) IsCaterpillar <- FALSE
        ## viewport of plot region
        if(!IsCaterpillar){
            bottomblank <- bottomblank + 2 ## the height of xaxis
        }
        pushViewport(viewport(x=LINEW + .5, y=bottomblank*LINEH/2 + .5, 
                              width= 1 - 7*LINEW,
                              height= 1 - bottomblank*LINEH,
                              xscale=c(start(ranges[i]), end(ranges[i])),
                              clip="off"))
        
        plot.grid.xaxis <- function(col="black"){
            ## axis, should be in the bottom of transcripts
            if(length(xaxis)==1 && as.logical(xaxis)) {
                grid.xaxis(gp=gpar(col=col))
            }
            if(length(xaxis)>1 && is.numeric(xaxis)){
                xaxisLabel <- names(xaxis)
                if(length(xaxisLabel)!=length(xaxis)) xaxisLabel <- TRUE
                grid.xaxis(at=xaxis, label=xaxisLabel, gp=gpar(col=col))
            }
        }
        ## plot xaxis
        bottomHeight <- 0
        if(IsCaterpillar){
            ## total height == maxscore + extension + gap + labels
            bottomHeight <- getHeight(SNPs=SNPs.bottom, 
                                      ratio.yx=ratio.yx, 
                                      LINEW=LINEW, 
                                      GAP=GAP, 
                                      cex=cex, 
                                      type=type,
                                      scoreMax=scoreMax,
                                      level="data&labels")
            vp <- viewport(y=bottomHeight, just="bottom",
                           xscale=c(start(ranges[i]), end(ranges[i])))
            pushViewport(vp)
            plot.grid.xaxis("gray")
            popViewport()
        }else{
            plot.grid.xaxis()
        }
        
        ## the baseline, the center of the first transcript
        baseline <- 
            max(c(feature.splited[[1]]$height/2, 
                  .0001)) + 0.2 * LINEH
        baselineN <- 
            max(c(feature.splited[[length(feature.splited)]]$height/2, 
                .0001)) + 0.2 * LINEH
        
        ##plot features
        feature.height <- plotFeatures(feature.splited, LINEH, bottomHeight)
        
        if(length(SNPs.bottom)>0){
            plotLollipops(SNPs.bottom, feature.height, bottomHeight, baselineN, 
                          type, ranges[i], yaxis, scoreMax, scoreMax0, scoreType, 
                          LINEW, cex, ratio.yx, GAP, pin, dashline.col,
                          side="bottom", jitter=jitter)
        }
        feature.height <- feature.height + 2*GAP
        if(length(SNPs.top)>0){
            plotLollipops(SNPs.top, feature.height, bottomHeight, baseline, 
                          type, ranges[i], yaxis, scoreMax, scoreMax0, scoreType, 
                          LINEW, cex, ratio.yx, GAP, pin, dashline.col,
                          side="top", jitter=jitter)
        }
        
        ## legend
        this.height <- getHeight(SNPs.top, 
                                 ratio.yx, LINEW, GAP, cex, type,
                                 scoreMax=scoreMax,
                                 level="data&labels")
        this.height <- this.height + bottomHeight + feature.height
        this.height <- plotLegend(legend[[i]], this.height, LINEH)
        
        popViewport()
        
        this.height <- bottomblank*LINEH + 
            this.height * (1 - bottomblank*LINEH)
        
        ## ylab
        vp <- viewport(x=.5, y=this.height*0.5, 
                       width=1, height=this.height)
        pushViewport(vp)
        if(is.logical(ylab)){
            if(ylab && length(names(SNP.gr))>0){
                grid.text(names(SNP.gr)[i], x = LINEW, 
                          y = .5, rot = 90)
            }
        }
        if(is.character(ylab)){
            if(length(ylab)==1) ylab <- rep(ylab, len)
            grid.text(ylab[i], x = LINEW,
                      y = .5, rot = 90)
        }
        popViewport()
        
        popViewport()
        height0 <-  height0 + this.height*height
    }
}

