#' dandelion.plots
#' @description Plot variants and somatic mutations
#' @param SNP.gr A object of \link[GenomicRanges:GRanges-class]{GRanges} or 
#' \link[GenomicRanges:GRangesList-class]{GRangesList}. 
#' All the width of GRanges must be 1.
#' @param features A object of \link[GenomicRanges:GRanges-class]{GRanges} or
#' \link[GenomicRanges:GRangesList-class]{GRangesList}.
#' @param ranges A object of \link[GenomicRanges:GRanges-class]{GRanges} or 
#' \link[GenomicRanges:GRangesList-class]{GRangesList}.
#' @param type Character. Could be fan, circle, pie or pin.
#' @param newpage plot in the new page or not.
#' @param ylab plot ylab or not. If it is a character vector, 
#' the vector will be used as ylab.
#' @param xaxis,yaxis plot xaxis/yaxis or not. If it is a numeric vector with length 
#' greater than 1, the vector will be used as the 
#' points at which tick-marks are to be drawn. And the names of the vector will be
#' used to as labels to be placed at the tick points if it has names. 
#' @param ylab.gp,xaxis.gp,yaxis.gp An object of class gpar for ylab, xaxis or yaxis.
#' @param legend If it is a list with named color vectors, a legend will be added.
#' @param cex cex will control the size of circle.
#' @param maxgaps maxgaps between the stem of dandelions. 
#' It is calculated by the width of plot region divided by maxgaps. 
#' If a GRanges object is set, the dandelions stem will be clustered
#' in each genomic range.
#' @param heightMethod A function used to determine the height of stem of 
#' dandelion. eg. Mean. Default is length.
#' @param label_on_feature Labels of the feature directly on them. 
#' Default FALSE.
#' @param ... not used.
#' @details In SNP.gr and features, metadata of the GRanges object will be used to 
#' control the color, fill, border, height, data source of pie if the type is pie.
#' @return NULL
#' @import GenomicRanges
#' @import grid
#' @importClassesFrom grImport Picture
#' @importFrom grImport readPicture grid.picture
#' @export
#' @examples
#' SNP <- c(10, 100, 105, 108, 400, 410, 420, 600, 700, 805, 840, 1400, 1402)
#' SNP.gr <- GRanges("chr1", IRanges(SNP, width=1, names=paste0("snp", SNP)), 
#'                   score=sample.int(100, length(SNP))/100)
#' features <- GRanges("chr1", IRanges(c(1, 501, 1001), 
#'                                     width=c(120, 500, 405),
#'                                     names=paste0("block", 1:3)),
#'                     color="black",
#'                     fill=c("#FF8833", "#51C6E6", "#DFA32D"),
#'                     height=c(0.1, 0.05, 0.08))
#' dandelion.plot(SNP.gr, features, type="fan")

dandelion.plot <- function(SNP.gr, features=NULL, ranges=NULL,
                      type=c("fan", "circle", "pie", "pin"),
                      newpage=TRUE, ylab=TRUE, ylab.gp=gpar(col="black"),
                      xaxis=TRUE, xaxis.gp=gpar(col="black"), 
                      yaxis=FALSE, yaxis.gp=gpar(col="black"), 
                      legend=NULL, cex=1, maxgaps=1/50, heightMethod=NULL, 
                      label_on_feature=FALSE, ...){
    stopifnot(inherits(SNP.gr, c("GRanges", "GRangesList")))
    stopifnot(inherits(features, c("GRanges", "GRangesList")))
    type <- match.arg(type)
    if(length(heightMethod)>0){
      stopifnot(is.function(heightMethod))
    }else{
      heightMethod <- length
    }
    if(is(maxgaps, "GRanges")){
      ol <- findOverlaps(maxgaps, drop.self=TRUE, drop.redundant=TRUE)
      if(length(ol)>0){
        stop("If maxgaps is an object of GRanges, maxgaps could not have overlaps.")
      }
    }
    if(type=="pin"){
        pinpath <- system.file("extdata", "map-pin-red.xml", package="trackViewer")
        pin <- readPicture(pinpath)
    }else{
        pin <- NULL
    }
    SNP.gr.name <- deparse(substitute(SNP.gr))
    if(is(SNP.gr, "GRanges")){
        SNP.gr <- list(SNP.gr)
        if(length(SNP.gr.name)==length(SNP.gr)) {
          names(SNP.gr) <- SNP.gr.name
        }
    }
    len <- length(SNP.gr)
    
    ## prepare the feature
    if(inherits(features, c("GRangesList", "list"))){
      for(i in seq_along(features)){
        stopifnot("features must be a GRanges or GRangesList object"=
                    is(features[[i]], "GRanges"))
      }
      features <- features[seq.int(len)]
    }else{
      stopifnot("features must be a GRanges or GRangesList object"=
                  is(features, "GRanges"))
      features <- list(features)[seq.int(len)]
    }
    
    ############### handle legend ####################
    ## set the legend as a list, 
    ## if all the legend for different tracks is same
    ## set draw legend for last track later
    legend <- handleLegend(legend, len)
    
    features.name <- deparse(substitute(features))
    
    ################ handle ranges #####################
    ## if !missing(ranges) set ranges as feature ranges
    ranges <- handleRanges(ranges, SNP.gr, features, len)
    
    ##cut all SNP.gr by the range
    SNP.gr <- cutSNP(SNP.gr, ranges, len)
    
    ################## plot ############################
    ## total height == 1
    height <- 1/sum(lengths(ranges))
    if(newpage) grid.newpage()
    for(i in seq.int(len)){
      if(length(ranges[[i]])>1){## split into multilayers
        args <- as.list(match.call())
        args$newpage <- FALSE
        vp <- viewport(x=.5, y=height*(i-0.5), width=1, height=height)
        pushViewport(vp)
        for(j in rev(seq_along(ranges[[i]]))){
          args$ranges <- ranges[[i]][j]
          args$SNP.gr <- SNP.gr[[i]]
          args$features <- features[[i]]
          args$legend <- legend[[i]]
          do.call(what = dandelion.plot, args = args)
        }
        popViewport()
      }else{## is GRanges with length==1
        vp <- viewport(x=.5, y=height*(i-0.5), width=1, height=height)
        pushViewport(vp)
        LINEW <- as.numeric(convertX(unit(1, "line"), "npc"))
        LINEH <- as.numeric(convertY(unit(1, "line"), "npc")) # LINEH is used for gap
        totalH <- as.numeric(unit(1, "npc"))
        if(LINEH > totalH/20){
          LINEH <- totalH/20
        }
        
        ## ylab
        if(length(yaxis)>1 && is.numeric(yaxis)){
          x <- LINEW
        }else{
          x <- unit(3, "lines")
          if(yaxis){
            x <- LINEW
          }
        }
        
        if(is.logical(ylab)){
          if(ylab && length(names(SNP.gr))>0){
            grid.text(names(SNP.gr)[i], x = x, 
                      y = .5, rot = 90, gp=ylab.gp)
          }
        }
        if(is.character(ylab)){
          if(length(ylab)==1) ylab <- rep(ylab, len)
          grid.text(ylab[i], x = x,
                    y = .5, rot = 90, gp=ylab.gp)
        }
        
        
        feature <- features[[i]]

        ## convert height to npc number
        feature$height <- convertHeight2NPCnum(feature$height)
        ## multiple transcripts in one gene could be separated by featureLayerID
        feature <- setFeatureLayerID(feature, ranges[[i]])
        feature.splited <- split(feature, feature$featureLayerID)
        
        baseline <- max(c(unlist(feature$height)/2, .0001)) + 0.2 * LINEH
        gap <- .2 * LINEH
        
        ## bottomblank, the transcripts legend height
        bottomblank <- plotFeatureLegend(feature, as.numeric(convertY(unit(1, "line"), "npc")), 
                                         ranges[[i]], xaxis, xaxis.gp, label_on_feature)
        
        pushViewport(viewport(x=LINEW + .5, y= bottomblank/2 + .5, 
                              width= 1 - 7*LINEW,
                              height= 1 - bottomblank,
                              xscale=c(start(ranges[[i]]), end(ranges[[i]])),
                              clip = "off"))
        ## axis
        plot_grid_xaxis(xaxis, gp=xaxis.gp)
        
        ## the baseline, the center of the first transcript
        baseline <- 
          max(c(feature.splited[[1]]$height/2, 
                .0001)) + 0.2 * LINEH
        baselineN <- 
          max(c(feature.splited[[length(feature.splited)]]$height/2, 
                .0001)) + 0.2 * LINEH
        
        ##plot features
        bottomHeight <- 0
        feature.height <- plotFeatures(feature.splited, LINEH, bottomHeight,
                                       label_on_feature)
        
        SNPs <- SNP.gr[[i]]
        if(length(SNPs)>0){
          strand(SNPs) <- "*"
          SNPs <- sort(SNPs)
          width <- feature.height + 2*gap
          SNPs.groups <- SNPs
          mcols(SNPs.groups) <- NULL
          SNPs.groups$w <- 0
          SNPs.groups$idx <- seq_along(SNPs)
          if(length(SNPs)==length(SNPs$score)){
            SNPs.groups$score <- SNPs$score
          }else{
            SNPs.groups$score <- 0
          }
          if(is(maxgaps, "GRanges")){
            SNPs.ol <- countOverlaps(maxgaps, SNPs)
            maxgaps <- maxgaps[SNPs.ol>0]
            if(length(maxgaps)<1){
              stop("Can not cluster the SNPs by given maxgaps")
            }
            SNPs.ol <- findOverlaps(maxgaps, SNPs)
            SNPs.groups$gps <- length(maxgaps) + seq_along(SNPs)
            SNPs.groups$gps[subjectHits(SNPs.ol)] <- queryHits(SNPs.ol)
          }else{
            if(!is.numeric(maxgaps)){
              stop("maxgaps must be a number or GRanges.")
            }
            SNPs.gap <- gaps(SNPs)
            SNPs.gap <- SNPs.gap[as.character(seqnames(SNPs.gap)) %in% as.character(seqnames(ranges[[i]])) & as.character(strand(SNPs.gap))=="*" & start(SNPs.gap)>=start(ranges[[i]]) & end(SNPs.gap)<=end(ranges[[i]])]
            SNPs.gap$w <- width(SNPs.gap)
            range.width <- floor(width(ranges[[i]])*maxgaps)
            SNPs.gap$idx <- rep(0, length(SNPs.gap))
            SNPs.gap$score <- 0
            SNPs.groups <- sort(c(SNPs.gap, SNPs.groups))
            SNPs.groups$gps <- cumsum(SNPs.groups$w >=range.width)
            SNPs.groups <- SNPs.groups[SNPs.groups$idx>0]
            SNPs.groups <- SNPs.groups[order(SNPs.groups$idx)]
          }
          if(length(names(SNPs))>0){
            maxStrHeight <- 
              max(as.numeric(
                convertX(stringWidth("W"), "npc")
              )*nchar(names(SNPs)))+LINEW/2
          }else{
            maxStrHeight <- 0
          }
          ratio.yx <- 1/as.numeric(convertX(unit(1, "snpc"), "npc"))
          ypos <- LINEW*max(ratio.yx, 1.2) + maxStrHeight*cex + LINEH*cex 
          if(length(legend[[i]])>0){
            if(!is.null(legend[[i]]$gp$cex)){
              ypos <- ypos + 1.5*LINEH*legend[[i]]$gp$cex
            }else{
              ypos <- ypos + 1.5*LINEH
            }
          }
          SNPs.groups <- Y1pos(SNPs.groups, c(start(ranges[[i]]), end(ranges[[i]])), LINEW, width, cex, 
                               ypos, 
                               length(yaxis) > 1 || (length(yaxis)==1 && as.logical(yaxis)), heightMethod)
          
          # yaxis
          yyscale <- c(0, SNPs.groups$yyscaleMax[1])
          if(length(yaxis)==1 && as.logical(yaxis)) {
            pushViewport(viewport(y= width + (1-width-ypos-cex*LINEW*ratio.yx)/2, 
                                  height= 1 - width-ypos-cex*LINEW*ratio.yx,
                                  yscale = yyscale))
            grid.yaxis(gp=yaxis.gp)
            popViewport()
          }
          if(length(yaxis)>1 && is.numeric(yaxis)){
            yaxisLabel <- names(yaxis)
            if(length(yaxisLabel)!=length(yaxis)) yaxisLabel <- TRUE
            pushViewport(viewport(y= width + (1-width-ypos-cex*LINEW*ratio.yx)/2, 
                                  height= 1 - width-ypos-cex*LINEW*ratio.yx,
                                  yscale = yyscale))
            grid.yaxis(at=yaxis, label=yaxisLabel, gp=yaxis.gp)
            popViewport()
          }
          ## legend
          scoreMax <- max(SNPs.groups$Y2, na.rm = TRUE)
          if(length(legend[[i]])>0){
            ypos <- LINEW*max(ratio.yx, 1.2) + 
              scoreMax + maxStrHeight*cex
            if(is.list(legend[[i]])){
              thisLabels <- legend[[i]][["labels"]]
              gp <- legend[[i]][names(legend[[i]])!="labels"]
              if(is.null(gp$cex)) gp$cex <- 1
              class(gp) <- "gpar"
            }else{
              thisLabels <- names(legend[[i]])
              gp <- gpar(fill=legend[[i]], cex=1) 
            }
            if(length(thisLabels)>0){
              ncol <- getColNum(thisLabels, cex=gp$cex)
              topblank <- ceiling(length(thisLabels) / ncol) * gp$cex
              pushViewport(viewport(x=.5, y=ypos+topblank*LINEH/2, 
                                    width=1,
                                    height=topblank*LINEH,
                                    just="bottom"))
              grid.legend(label=thisLabels, ncol=ncol,
                          byrow=TRUE, vgap=unit(.1*gp$cex, "lines"),
                          hgap=unit(.5*gp$cex, "lines"),
                          pch=21,
                          gp=gp)
              popViewport()
            }
          }
          for(m in 1:length(SNPs)){
            this.dat <- SNPs[m]
            this.dat.grp <- SNPs.groups[m]
            color <- if(is.list(this.dat$color)) this.dat$color[[1]] else this.dat$color
            border <- if(is.list(this.dat$border)) this.dat$border[[1]] else this.dat$border
            fill <- if(is.list(this.dat$fill)) this.dat$fill[[1]] else this.dat$fill
            lwd <- if(is.list(this.dat$lwd)) this.dat$lwd[[1]] else this.dat$lwd
            id <- if(is.character(this.dat$label)) this.dat$label else NA
            id.col <- if(length(this.dat$label.col)>0) this.dat$label.col else "black"
            this.dat.mcols <- mcols(this.dat)
            this.dat.mcols <- cleanDataMcols(this.dat.mcols, type)
            grid.dandelion(x0=(start(this.dat)-start(ranges[[i]]))/width(ranges[[i]]), 
                           y0=baseline,
                           x1=(start(this.dat)-start(ranges[[i]]))/width(ranges[[i]]),
                           y1=this.dat.grp$Y1,
                           x2=(this.dat.grp$X2-start(ranges[[i]]))/width(ranges[[i]]), 
                           y2=this.dat.grp$Y2, 
                           radius=cex*LINEW/2,
                           col=color,
                           border=border,
                           percent=this.dat.mcols,
                           edges=100,
                           alpha=this.dat.grp$alpha,
                           type=type,
                           ratio.yx=ratio.yx,
                           pin=pin,
                           scoreMax=scoreMax,
                           id=id, id.col=id.col, 
                           name=names(this.dat), 
                           cex=cex, lwd=lwd)
          }
        }
        
        popViewport()
        popViewport()
      }
    }
}