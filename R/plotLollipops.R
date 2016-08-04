convertHeight2NPCnum <- function(.ele){
    switch(class(.ele),
           "unit"=convertHeight(.ele, unitTo="npc", valueOnly=TRUE),
           "list"={
               .ele <- sapply(.ele, function(.e){
                   if(class(.e)=="unit"){
                       .e <- convertHeight(.e, unitTo="npc", valueOnly=TRUE)
                   }
                   .e[1]
               })
               unlist(.ele)
           },
           "numeric"=.ele,
           "integer"=.ele,
           .ele)
}
plotFeatures <- function(feature.splited, LINEH, bottomHeight){
    feature.height <- 0
    for(n in 1:length(feature.splited)){
        this.feature.height <- 
            max(c(feature.splited[[n]]$height/2, 
                  .0001)) + 0.2 * LINEH
        feature.height <- feature.height + this.feature.height
        ##baseline
        grid.lines(x=c(0, 1), y=c(bottomHeight+feature.height, 
                                  bottomHeight+feature.height))
        for(m in 1:length(feature.splited[[n]])){
            this.dat <- feature.splited[[n]][m]
            color <- if(is.list(this.dat$color)) this.dat$color[[1]] else 
                this.dat$color
            fill <- if(is.list(this.dat$fill)) this.dat$fill[[1]] else 
                this.dat$fill
            this.cex <- if(length(this.dat$cex)>0) this.dat$cex[[1]][1] else 1
            this.feature.height.m <- 
                if(length(this.dat$height)>0) 
                    this.dat$height[[1]][1] else 
                        2*this.feature.height
            grid.rect(x=start(this.dat), y=bottomHeight+feature.height, 
                      width=width(this.dat), 
                      height=this.feature.height.m,
                      just="left", gp=gpar(col=color, fill=fill), 
                      default.units = "native")
        }
        feature.height <- feature.height + this.feature.height
    }
    feature.height
}

plotLollipops <- function(SNPs, feature.height, bottomHeight, baseline, 
                          type, ranges, yaxis, scoreMax, scoreMax0, scoreType,
                          LINEW, cex, ratio.yx, GAP, pin, dashline.col,
                          side=c("top", "bottom")){
    side <- match.arg(side)
    if(side=="top"){
        pushViewport(viewport(y=bottomHeight,
                              height=1,
                              just="bottom",
                              xscale=c(start(ranges), 
                                       end(ranges)),
                              clip="off"))
    }else{
        pushViewport(viewport(y=bottomHeight+feature.height,
                              height=1,
                              just="top",
                              xscale=c(start(ranges), 
                                       end(ranges)),
                              yscale=c(1, 0),
                              clip="off"))
    }
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
                                         xscale=c(start(ranges), end(ranges)), 
                                         lineW=LINEW*cex)
        condense.ids <- SNPs.condense$oid
        lab.pos <- rep(lab.pos.condense, elementNROWS(condense.ids))
        lab.pos <- lab.pos[order(unlist(condense.ids))]
    }else{
        lab.pos <- jitterLables(start(SNPs), 
                                xscale=c(start(ranges), end(ranges)), 
                                lineW=LINEW*cex)
    }
    
    if(length(SNPs)>0){
        if(yaxis && scoreMax>1 && !type %in% c("pie", "pie.stack")){
            if(side=="top"){
                grid.yaxis(vp=viewport(x=.5-LINEW,
                                       y=feature.height+5.25*GAP*cex+
                                           scoreMax*LINEW*ratio.yx/2*cex,
                                       width=1,
                                       height=scoreMax*LINEW*ratio.yx*cex,
                                       yscale=c(0, scoreMax0+.5)))
            }else{
                grid.yaxis(vp=viewport(x=.5-LINEW,
                                       y=1-(feature.height+5.25*GAP*cex+
                                           scoreMax*LINEW*ratio.yx/2*cex),
                                       width=1,
                                       height=scoreMax*LINEW*ratio.yx*cex,
                                       yscale=c(scoreMax0+.5, 0)))
            }
        }
        for(m in 1:length(SNPs)){
            this.dat <- SNPs[m]
            color <- if(is.list(this.dat$color)) this.dat$color[[1]] else this.dat$color
            border <- 
                if(is.list(this.dat$border)) this.dat$border[[1]] else this.dat$border
            fill <- if(is.list(this.dat$fill)) this.dat$fill[[1]] else this.dat$fill
            lwd <- if(is.list(this.dat$lwd)) this.dat$lwd[[1]] else this.dat$lwd
            id <- if(is.character(this.dat$label)) this.dat$label else NA
            id.col <- if(length(this.dat$label.col)>0) this.dat$label.col else "black"
            this.dat.mcols <- mcols(this.dat)
            this.dat.mcols <- 
                this.dat.mcols[, 
                               !colnames(this.dat.mcols) %in% 
                                   c("color", "fill", "lwd", "id", 
                                     "id.col", "stack.factor", "SNPsideID"), 
                               drop=FALSE]
            if(type!="pie.stack"){
                this.dat.mcols <- 
                    this.dat.mcols[, !colnames(this.dat.mcols) %in% 
                                       c("stack.factor.order", 
                                         "stack.factor.first"), 
                                   drop=FALSE]
            }
            this.dat.mcols <- 
                this.dat.mcols[, !grepl("^label.parameter",
                                        colnames(this.dat.mcols)), 
                               drop=FALSE]

            grid.lollipop(x1=(start(this.dat)-start(ranges))/width(ranges), 
                          y1=baseline,
                          x2=(lab.pos[m]-start(ranges))/width(ranges), 
                          y2=feature.height,
                          y3=4*GAP*cex, y4=2.5*GAP*cex, 
                          radius=cex*LINEW/2,
                          col=color,
                          border=border,
                          percent=this.dat.mcols,
                          edges=100,
                          type=type,
                          ratio.yx=ratio.yx,
                          pin=pin,
                          scoreMax=(scoreMax-0.5) * LINEW * cex,
                          scoreType=scoreType,
                          id=id, id.col=id.col,
                          cex=cex, lwd=lwd, dashline.col=dashline.col,
                          side=side)

        }
        this.height <- getHeight(SNPs, 
                                 ratio.yx, LINEW, GAP, cex, type,
                                 level="data")
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
            labels.just <- ifelse(side=="top", "left", "right")
            labels.hjust <- NULL
            labels.vjust <- NULL
            labels.check.overlap <- FALSE
            labels.default.units <- "native"
            labels.gp <- gpar(cex=cex)
            
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
            
            grid.text(x=labels.x, y=this.height + feature.height, 
                      label = labels.text,  
                      just = labels.just, 
                      hjust = labels.hjust,
                      vjust = labels.vjust,
                      rot=labels.rot,
                      check.overlap = labels.check.overlap,
                      default.units = labels.default.units,
                      gp=labels.gp)
        }
    }
    popViewport()
}

plotLegend <- function(legend, this.height, LINEH){
    ypos <- this.height
    if(length(legend)>0){
        if(is.list(legend)){
            thisLabels <- legend[["labels"]]
            gp <- legend[names(legend)!="labels"]
            class(gp) <- "gpar"
        }else{
            thisLabels <- names(legend)
            gp <- gpar(fill=legend) 
        }
        if(length(thisLabels)>0){
            ncol <- getColNum(thisLabels)
            topblank <- ceiling(length(thisLabels) / ncol)
            pushViewport(viewport(x=.5, 
                                  y=ypos+(topblank+.5)*LINEH/2, 
                                  width=1,
                                  height=topblank*LINEH,
                                  just="bottom"))
            this.height <- ypos + (topblank+1)*LINEH 
            grid.legend(label=thisLabels, ncol=ncol,
                        byrow=TRUE, vgap=unit(.2, "lines"),
                        pch=21,
                        gp=gp)
            popViewport()
        }
    }
    this.height
}