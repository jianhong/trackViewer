lolliplot <- function(SNP.gr, features=NULL, ranges=NULL,
                      type=c("circle", "pie", "pin", "pie.stack"),
                      newpage=TRUE, ylab=TRUE, yaxis=TRUE,
                      xaxis=TRUE, legend=NULL, cex=1, dashline.col="gray80", ...){
    stopifnot(inherits(SNP.gr, c("GRanges", "GRangesList")))
    stopifnot(inherits(features, c("GRanges", "GRangesList")))
    type <- match.arg(type)
    if(type=="pin"){
        pinpath <- system.file("extdata", "map-pin-red.xml", package="trackViewer")
        pin <- readPicture(pinpath)
    }else{
        pin <- NULL
    }
    SNP.gr.name <- deparse(substitute(SNP.gr))
    if(class(SNP.gr)=="GRanges"){
        SNP.gr <- GRangesList(SNP.gr)
        names(SNP.gr) <- SNP.gr.name
    }
    len <- length(SNP.gr)
    if(length(legend)>0){
        if(!is.list(legend)){
            tmp <- legend
            legend <- list()
            length(legend) <- len
            legend[[len]] <- tmp
            rm(tmp)
        }else{
            if(length(legend)==1){
                tmp <- legend[[1]]
                legend <- list()
                length(legend) <- len
                legend[[len]] <- tmp
                rm(tmp)
            }else{
                if("labels" %in% names(legend)){
                    tmp <- legend
                    legend <- list()
                    length(legend) <- len
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
    if(length(ranges)>0){
        stopifnot(class(ranges)=="GRanges"&length(ranges)==length(SNP.gr))
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
    height <- 1/len
    if(newpage) grid.newpage()
    for(i in 1:len){
        vp <- viewport(x=.5, y=height*(i-0.5), width=1, height=height)
        pushViewport(vp)
        lineW <- as.numeric(convertX(unit(1, "line"), "npc"))
        lineH <- as.numeric(convertY(unit(1, "line"), "npc"))
        ## ylab
        if(is.logical(ylab)){
            if(ylab && length(names(SNP.gr))>0){
                grid.text(names(SNP.gr)[i], x = lineW, 
                          y = .5, rot = 90)
            }
        }
        if(is.character(ylab)){
            if(length(ylab)==1) ylab <- rep(ylab, len)
            grid.text(ylab[i], x = lineW,
                      y = .5, rot = 90)
        }
        
        if(class(features)=="GRangesList"){
            feature <- features[[i]]
        }else{
            feature <- features
        }
        start(feature)[start(feature)<start(ranges[i])] <- start(ranges[i])
        end(feature)[end(feature)>end(ranges[i])] <- end(ranges[i])
        baseline <- max(c(unlist(feature$height)/2, .0001)) + 0.2 * lineH
        gap <- .2 * lineH
        bottomblank <- 4
        if(length(names(feature))>0){
            feature.s <- feature[!duplicated(names(feature))]
            ncol <- getColNum(names(feature.s))
            bottomblank <- max(ceiling(length(names(feature.s)) / ncol), 4)
            pushViewport(viewport(x=.5, y=bottomblank*lineH/2, 
                                  width=1,
                                  height=bottomblank*lineH,
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
        }
        pushViewport(viewport(x=lineW + .5, y= (bottomblank+2)*lineH/2 + .5, 
                              width= 1 - 7*lineW,
                              height= 1 - (bottomblank+2)*lineH,
                              xscale=c(start(ranges[i]), end(ranges[i])),
                              clip="off"))
        ## axis
        if(length(xaxis)==1 && as.logical(xaxis)) {
            grid.xaxis()
        }
        if(length(xaxis)>1 && is.numeric(xaxis)){
            xaxisLabel <- names(xaxis)
            if(length(xaxisLabel)!=length(xaxis)) xaxisLabel <- TRUE
            grid.xaxis(at=xaxis, label=xaxisLabel)
        }
        
        ##baseline
        grid.lines(x=c(0, 1), y=c(baseline, baseline)) #baseline
        
        for(m in 1:length(feature)){
            this.dat <- feature[m]
            color <- if(is.list(this.dat$color)) this.dat$color[[1]] else this.dat$color
            fill <- if(is.list(this.dat$fill)) this.dat$fill[[1]] else this.dat$fill
            this.cex <- if(length(this.dat$cex)>0) this.dat$cex[[1]][1] else 1
            width <- if(length(this.dat$height)>0) this.dat$height[[1]][1] else 2*baseline
            rot <- if(length(this.dat$rot)>0) this.dat$rot[[1]][1] else 45
            grid.rect(x=start(this.dat), y=baseline, width=width(this.dat), height=width,
                      just="left", gp=gpar(col=color, fill=fill), default.units = "native")
        }
        SNPs <- SNP.gr[[i]]
        strand(SNPs) <- "*"
        SNPs <- sort(SNPs)
        width <- 2 * baseline + 2*gap
        if(type=="pie.stack" && length(SNPs$stack.factor)>0){
            stopifnot(is.vector(SNPs$stack.factor, mode="character"))
            if(length(SNPs$stack.factor.order)>0 || 
               length(SNPs$stack.factor.first)>0){
                warning("stack.factor.order and stack.factor.first are used by this function!",
                        "The values in these column will be removed.")
            }
            ## condense the SNPs
            stack.factors <- unique(as.character(SNPs$stack.factor))
            stack.factors <- sort(stack.factors)
            stack.factors.order <- 1:length(stack.factors)
            names(stack.factors.order) <- stack.factors
            SNPs <- SNPs[order(as.character(seqnames(SNPs)), start(SNPs), 
                               as.character(SNPs$stack.factor))]
            SNPs$stack.factor.order <- stack.factors.order[SNPs$stack.factor]
            SNPs$stack.factor.first <- !duplicated(SNPs)
            SNPs.condense <- SNPs
            SNPs.condense$oid <- 1:length(SNPs)
            SNPs.condense$factor <- paste(as.character(seqnames(SNPs)), start(SNPs), end(SNPs))
            SNPs.condense <- split(SNPs.condense, SNPs.condense$factor)
            SNPs.condense <- lapply(SNPs.condense, function(.ele){
                .oid <- .ele$oid
                .gr <- .ele[1]
                mcols(.gr) <- NULL
                .gr$oid <- NumericList(.oid)
                .gr
            })
            SNPs.condense <- unlist(GRangesList(SNPs.condense), use.names = FALSE)
            SNPs.condense <- sort(SNPs.condense)
            lab.pos.condense <- jitterLables(start(SNPs.condense), 
                                    xscale=c(start(ranges[i]), end(ranges[i])), 
                                    lineW=lineW*cex)
            condense.ids <- SNPs.condense$oid
            lab.pos <- rep(lab.pos.condense, elementNROWS(condense.ids))
            lab.pos <- lab.pos[order(unlist(condense.ids))]
        }else{
            lab.pos <- jitterLables(start(SNPs), 
                                    xscale=c(start(ranges[i]), end(ranges[i])), 
                                    lineW=lineW*cex)
        }
        
        if(length(SNPs)>0){
            scoreMax0 <- scoreMax <- if(length(SNPs$score)>0) ceiling(max(c(SNPs$score, 1), na.rm=TRUE)) else 1
            if(type=="pie.stack") scoreMax <- length(stack.factors)
            if(!type %in% c("pie", "pie.stack")){
                if(scoreMax>10) {
                    SNPs$score <- 10*SNPs$score/scoreMax 
                    scoreMax <- ceiling(max(c(SNPs$score, 1), na.rm=TRUE))
                }
                scoreType <- if(length(SNPs$score)>0) all(floor(SNPs$score)==SNPs$score) else FALSE
            }else{
                scoreType <- FALSE
            }
            
            ratio.yx <- 1/as.numeric(convertX(unit(1, "snpc"), "npc"))
            if(yaxis && scoreMax>1 && !type %in% c("pie", "pie.stack")){
                grid.yaxis(vp=viewport(x=.5-lineW,
                                       y=width+5.25*gap*cex+scoreMax*lineW*ratio.yx/2*cex,
                                       width=1,
                                       height=scoreMax*lineW*ratio.yx*cex,
                                       yscale=c(0, scoreMax0+.5)))
            }
            for(m in 1:length(SNPs)){
                this.dat <- SNPs[m]
                color <- if(is.list(this.dat$color)) this.dat$color[[1]] else this.dat$color
                border <- if(is.list(this.dat$border)) this.dat$border[[1]] else this.dat$border
                fill <- if(is.list(this.dat$fill)) this.dat$fill[[1]] else this.dat$fill
                lwd <- if(is.list(this.dat$lwd)) this.dat$lwd[[1]] else this.dat$lwd
                id <- if(is.character(this.dat$label)) this.dat$label else NA
                id.col <- if(length(this.dat$label.col)>0) this.dat$label.col else "black"
                this.dat.mcols <- mcols(this.dat)
                this.dat.mcols <- this.dat.mcols[, !colnames(this.dat.mcols) %in% c("color", "fill", "lwd", "id", "id.col", "stack.factor"), drop=FALSE]
                if(type!="pie.stack"){
                    this.dat.mcols <- this.dat.mcols[, !colnames(this.dat.mcols) %in% c("stack.factor.order", "stack.factor.first"), drop=FALSE]
                }
                this.dat.mcols <- this.dat.mcols[, !grepl("^label.parameter", colnames(this.dat.mcols)), drop=FALSE]
                grid.lollipop(x1=(start(this.dat)-start(ranges[i]))/width(ranges[i]), 
                              y1=baseline,
                              x2=(lab.pos[m]-start(ranges[i]))/width(ranges[i]), y2=width,
                              y3=4*gap*cex, y4=2.5*gap*cex, 
                              radius=cex*lineW/2,
                              col=color,
                              border=border,
                              percent=this.dat.mcols,
                              edges=100,
                              type=type,
                              ratio.yx=ratio.yx,
                              pin=pin,
                              scoreMax=(scoreMax-0.5) * lineW * cex,
                              scoreType=scoreType,
                              id=id, id.col=id.col,
                              cex=cex, lwd=lwd, dashline.col=dashline.col)
            }
            labels.rot <- 90
            if(length(names(SNPs))>0){
                if(type=="pie.stack"){
                    ## unique lab.pos and SNPs
                    idx <- !duplicated(names(SNPs))
                    lab.pos <- lab.pos[idx]
                    SNPs <- SNPs[idx]
                }
                labels.x <- lab.pos
                labels.text <- names(SNPs)
                labels.just <- "left"
                labels.hjust <- NULL
                labels.vjust <- NULL
                labels.check.overlap <- FALSE
                labels.default.units <- "native"
                labels.gp <- gpar(cex=cex)
                switch(type,
                       circle={
                           labels.y <- width + lineW*max(ratio.yx, 1.2) + 
                               6.5*gap*cex + 
                               (scoreMax-0.5) * lineW * ratio.yx*cex
                       },
                       pin={
                           this.scores <- if(length(SNPs$score)>0) ceiling(SNPs$score) else .5
                           this.scores[is.na(this.scores)] <- .5
                           labels.y <- width + lineW*max(ratio.yx, 1.2) + 
                               6.5*gap*cex + 
                               (this.scores-0.5) * lineW * ratio.yx*cex
                       },
                       pie={
                           labels.y <- width + lineW*max(ratio.yx, 1.2) + 
                               6.5*gap*cex + 0.5 * lineW * ratio.yx * cex
                       },
                       pie.stack={
                           labels.y <- width + lineW*max(ratio.yx, 1.2) + 
                               6.5*gap*cex + 
                               (scoreMax-0.5) * lineW * ratio.yx*cex
                       })
                ## change the parameter by use definations.
                for(label.parameter in c("x", "y", "just", "hjust", "vjust",
                                         "rot", "check.overlap", "default.units",
                                         "gp")){
                    label.para <- paste0("label.parameter.", label.parameter)
                    if(label.para %in% colnames(mcols(SNPs))){
                        assign(paste0("labels.", label.parameter), 
                               mcols(SNPs)[, label.para])
                    }
                }
                labels.gp <- c(labels.gp, cex=cex)
                labels.gp[duplicated(names(labels.gp))] <- NULL
                labels.gp <- do.call(gpar, labels.gp)
                
                grid.text(x=labels.x, y=labels.y, 
                          label = labels.text,  
                          just=labels.just, 
                          hjust = labels.hjust,
                          vjust = labels.vjust,
                          rot=labels.rot,
                          check.overlap = labels.check.overlap,
                          default.units = labels.default.units,
                          gp=labels.gp)
            }
            ## legend
            labels.length.rate <- max(cospi((labels.rot-90)/180), 0)
            if(length(legend[[i]])>0){
                switch(type,
                       circle={
                           if(length(names(SNPs))>0){
                               maxStrHeight <- 
                                   max(as.numeric(
                                       convertY(stringWidth(names(SNPs)), "npc")
                                       ))+lineW/2
                           }else{
                               maxStrHeight <- 0
                           }
                           maxStrHeight <- maxStrHeight * labels.length.rate
                           ypos <- width + lineW*max(ratio.yx, 1.2) + 6.5*gap*cex + 
                               (scoreMax-0.5) * lineW * ratio.yx*cex + maxStrHeight*cex
                       },
                       pin={
                           if(length(names(SNPs))>0){
                               thisStrHeight <- as.numeric(
                                   convertY(stringWidth(names(SNPs)), "npc"))+lineW/2
                           }else{
                               thisStrHeight <- 0
                           }
                           thisStrHeight <- thisStrHeight * labels.length.rate
                           if(length(SNPs$score)>0){
                               ypos <- 
                                   max(width + lineW*max(ratio.yx, 1.2) + 
                                           6.5*gap*cex + 
                                           (SNPs$score-0.5) * lineW * ratio.yx*cex + 
                                           thisStrHeight*cex)
                           }else{
                               ypos <- max(width + lineW*max(ratio.yx, 1.2) + 
                                               6.5*gap*cex + thisStrHeight*cex)
                           }
                       },
                       pie={
                           if(length(names(SNPs))>0){
                               maxStrHeight <- 
                                   max(as.numeric(
                                       convertY(stringWidth(names(SNPs)), "npc")
                                   ))+lineW/2
                           }else{
                               maxStrHeight <- 0
                           }
                           maxStrHeight <- maxStrHeight * labels.length.rate
                           ypos <- width + lineW*max(ratio.yx, 1.2) + 
                               6.5*gap*cex + maxStrHeight*cex
                       },
                       pie.stack={
                           if(length(names(SNPs))>0){
                               maxStrHeight <- 
                                   max(as.numeric(
                                       convertY(stringWidth(names(SNPs)), "npc")
                                   ))+lineW/2
                           }else{
                               maxStrHeight <- 0
                           }
                           maxStrHeight <- maxStrHeight * labels.length.rate
                           ypos <- width + lineW*max(ratio.yx, 1.2) + 
                               6.5*gap*cex + maxStrHeight*cex +
                               (scoreMax-0.5) * lineW * ratio.yx*cex
                       }
                       )
                if(is.list(legend[[i]])){
                    thisLabels <- legend[[i]][["labels"]]
                    gp <- legend[[i]][names(legend[[i]])!="labels"]
                    class(gp) <- "gpar"
                }else{
                    thisLabels <- names(legend[[i]])
                    gp <- gpar(fill=legend[[i]]) 
                }
                if(length(thisLabels)>0){
                    ncol <- getColNum(thisLabels)
                    topblank <- ceiling(length(thisLabels) / ncol)
                    pushViewport(viewport(x=.5, y=ypos+topblank*lineH/2, 
                                          width=1,
                                          height=topblank*lineH,
                                          just="bottom"))
                    grid.legend(label=thisLabels, ncol=ncol,
                                byrow=TRUE, vgap=unit(.2, "lines"),
                                pch=21,
                                gp=gp)
                    popViewport()
                }
            }
        }
        
        popViewport()
        popViewport()
    }
}