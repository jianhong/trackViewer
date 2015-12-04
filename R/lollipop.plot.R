lolliplot <- function(SNP.gr, features=NULL, ranges=NULL,
                          type=c("circle", "pie", "pin"),
                          newpage=TRUE, ylab=TRUE, yaxis=TRUE){
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
            sW <- max(as.numeric(convertX(stringWidth(names(feature)), "line"))) + 1
            ncol <- floor(as.numeric(convertX(unit(1, "npc"), "line")) / 
                              sW / as.numeric(convertX(stringWidth("W"), "line")))
            bottomblank <- max(ceiling(length(names(feature)) / ncol), 4)
            pushViewport(viewport(x=.5, y=bottomblank*lineH/2, 
                                  width=1,
                                  height=bottomblank*lineH,
                                  xscale=c(start(ranges[i]), end(ranges[i]))))
            color <- if(length(unlist(feature$color))==length(feature)) 
                unlist(feature$color) else "black"
            fill <- if(length(unlist(feature$fill))==length(feature)) 
                unlist(feature$fill) else "black"
            pch <- if(length(unlist(feature$pch))==length(feature)) 
                unlist(feature$pch) else 22
            grid.legend(label=names(feature), ncol=ncol,
                      byrow=TRUE, vgap=unit(.2, "lines"),
                      pch=pch,
                      gp=gpar(col=color, fill=fill))
            popViewport()
        }
        pushViewport(viewport(x=lineW + .5, y= (bottomblank+2)*lineH/2 + .5, 
                              width= 1 - 7*lineW,
                              height= 1 - (bottomblank+2)*lineH,
                              xscale=c(start(ranges[i]), end(ranges[i]))))
        ## axis
        grid.xaxis()
        
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
        lab.pos <- jitterLables(start(SNPs), 
                                xscale=c(start(ranges[i]), end(ranges[i])), 
                                lineW=lineW)
        if(length(SNPs)>0){
            scoreMax0 <- scoreMax <- if(length(SNPs$score)>0) ceiling(max(c(SNPs$score, 1), na.rm=TRUE)) else 1
            if(scoreMax>10) {
                SNPs$score <- 10*SNPs$score/scoreMax 
                scoreMax <- ceiling(max(c(SNPs$score, 1), na.rm=TRUE))
            }
            scoreType <- if(length(SNPs$score)>0) all(floor(SNPs$score)==SNPs$score) else FALSE
            ratio.yx <- 1/as.numeric(convertX(unit(1, "snpc"), "npc"))
            if(yaxis && 
               scoreMax>1 && ((!scoreType) || type=="pin") && type!="pie"){
                grid.yaxis(vp=viewport(x=.5-lineW,
                                       y=width+5.25*gap+scoreMax*lineW*ratio.yx/2,
                                       width=1,
                                       height=scoreMax*lineW*ratio.yx,
                                       yscale=c(0, scoreMax0+.5)))
            }
            for(m in 1:length(SNPs)){
                this.dat <- SNPs[m]
                color <- if(is.list(this.dat$color)) this.dat$color[[1]] else this.dat$color
                border <- if(is.list(this.dat$border)) this.dat$border[[1]] else this.dat$border
                fill <- if(is.list(this.dat$fill)) this.dat$fill[[1]] else this.dat$fill
                id <- if(is.character(this.dat$label)) this.dat$label else NA
                id.col <- if(length(this.dat$label.col)>0) this.dat$label.col else "black"
                grid.lollipop(x1=(start(this.dat)-start(ranges[i]))/width(ranges[i]), 
                              y1=baseline,
                              x2=(lab.pos[m]-start(ranges[i]))/width(ranges[i]), y2=width,
                              y3=4*gap, y4=2.5*gap, 
                              radius=lineW/2,
                              col=color,
                              border=border,
                              percent=mcols(this.dat),
                              edges=100,
                              type=type,
                              ratio.yx=ratio.yx,
                              pin=pin,
                              scoreMax=(scoreMax-0.5) * lineW,
                              scoreType=scoreType,
                              id=id, id.col=id.col)
            }
            if(length(names(SNPs))>0){
                switch(type,
                       circle={
                           grid.text(x=lab.pos, 
                                     y=width + lineW*max(ratio.yx, 1.2) + 6.5*gap + (scoreMax-0.5) * lineW * ratio.yx, 
                                     label = names(SNPs), rot=90, just="left", 
                                     default.units = "native")
                       },
                       pin={
                           this.scores <- if(length(SNPs$score)>0) ceiling(SNPs$score) else .5
                           this.scores[is.na(this.scores)] <- .5
                           grid.text(x=lab.pos, 
                                     y=width + lineW*max(ratio.yx, 1.2) + 6.5*gap + (this.scores-0.5) * lineW * ratio.yx, 
                                     label = names(SNPs), rot=90, just="left", 
                                     default.units = "native")
                       },
                       pie={
                           grid.text(x=lab.pos, y=width + lineW*max(ratio.yx, 1.2) + 6.5*gap, 
                                     label = names(SNPs), rot=90, just="left", 
                                     default.units = "native")
                       })
            }
        }
        
        popViewport()
        popViewport()
    }
}