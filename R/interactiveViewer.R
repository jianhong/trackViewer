
interactiveViewer <- function(trackList, chromosome, start, end, strand, 
                              viewerStyle=trackViewerStyle(), autoOptimizeStyle=FALSE){
    savepar <- par(ask=FALSE)
    on.exit(par(savepar))
    txdb <- NULL
    txdump <- NULL
    xscaleold <- NULL
    if(missing(trackList)){
        #stop("trackList is required.")
        importData <- function(){
            checkformat <- function(filename){
                format <- c(BAM="bam", BED="bed", BigWig="(bw)|(bigWig)", WIG="wig", bedGraph="bedGraph")
                for(i in 1:length(format)){
                    if(grepl(paste(format[i], '$', sep=""), filename, ignore.case=TRUE))
                        return(names(format)[i])
                }
            }
            window <- gwindow("import data", width=300, height=300)
            group <- ggroup(container=window, horizontal=FALSE)
            info.group <- ggroup(horizontal=FALSE, container=group)
            tmp <- ggroup(horizontal=TRUE, container=info.group)
            flab <- glabel("file path: ", container =tmp)
            filename <- NULL
            gfb<- gfilebrowse("Select a file", 
                                    type="open", 
                                    filter = list("Track files" = list(patterns = 
                                                                           c("*.bam","*.bed", "*.wig", "*.bw", "*.bedGraph")),
                                                  "All files" = list(patterns = c("*"))),
                                    multi=FALSE, 
                                    container=tmp,
                              handler = function(h, ...){
                                  if(!is.na(svalue(h$obj))){
                                      if(length(filename)==0){
                                          filename <<- svalue(gfb)
                                          svalue(selectedfilesFormat) <- checkformat(svalue(gfb))
                                          selectedfiles[1] <- filename
                                      }else{
                                          if(!svalue(gfb) %in% filename){
                                              filename <<- c(filename, svalue(gfb))
                                              tmpselgp <- ggroup(horizontal=TRUE, container=tmpsel)
                                              format <- checkformat(svalue(gfb))
                                              selectedfilesFormat <<- c(selectedfilesFormat,
                                                                        gedit(format, container=tmpselgp, width=6))
                                              selectedfiles <<- c(selectedfiles, 
                                                                  gcheckbox(svalue(gfb), 
                                                                            checked=TRUE,
                                                                            container=tmpselgp))
                                          }
                                      }
                                  }
                              })
            tmpsel <- gframe("import files", container=info.group, expand=TRUE, horizontal=FALSE)
            tmpselgp <- ggroup(horizontal=TRUE, container=tmpsel)
            selectedfilesFormat <- gedit("format", container=tmpselgp, width=6)
            selectedfiles <- gcheckbox("NULL", checked=TRUE, container=tmpselgp)
            
            tmp <- gframe("genomic coordinates", container=info.group, expand=TRUE, horizontal=FALSE)
            availGenome <- c("TxDb.Athaliana.BioMart.plantsmart10", "TxDb.Athaliana.BioMart.plantsmart12", 
                             "TxDb.Athaliana.BioMart.plantsmart14", "TxDb.Athaliana.BioMart.plantsmart16", 
                             "TxDb.Athaliana.BioMart.plantsmart19", "TxDb.Celegans.UCSC.ce6.ensGene", 
                             "TxDb.Dmelanogaster.UCSC.dm3.ensGene", "TxDb.Hsapiens.UCSC.hg18.knownGene", 
                             "TxDb.Hsapiens.UCSC.hg19.knownGene", "TxDb.Hsapiens.UCSC.hg19.lincRNAsTranscripts", 
                             "TxDb.Mmusculus.UCSC.mm10.ensGene", "TxDb.Mmusculus.UCSC.mm10.knownGene", 
                             "TxDb.Mmusculus.UCSC.mm9.knownGene", "TxDb.Rnorvegicus.UCSC.rn4.ensGene", 
                             "TxDb.Rnorvegicus.UCSC.rn5.refGene", "TxDb.Scerevisiae.UCSC.sacCer2.sgdGene", 
                             "TxDb.Scerevisiae.UCSC.sacCer3.sgdGene")
            genome <- gradio(availGenome, selected=9, container=tmp, horizontal=FALSE)
            tmp <- gframe("genome", container=info.group, expand=TRUE, horizontal=FALSE) 
            ltmp <- ggroup(horizontal=TRUE, container=tmp)
            ltmpl <- ggroup(horizontal=FALSE, container=ltmp)
            label <- glabel("chromosome", container =ltmpl)
            label <- glabel("start: ", container =ltmpl)
            label <- glabel("end: ", container =ltmpl)
            label <- glabel("strand: ", container =ltmpl)
            ltmpr <- ggroup(horizontal=FALSE, container=ltmp)
            chr <- gedit(chromosome, container=ltmpr, width=8)
            xaxislim1 <- gedit(start, container=ltmpr, width=16)
            xaxislim2 <- gedit(end, container=ltmpr, width=16)
            availStrand <- c("+", "-", "*")
            strand <- gradio(availStrand, selected=which(availStrand==strand),
                             container=ltmpr, width=8, horizontal=TRUE)
            
            
            button.group <- ggroup(horizontal=TRUE, container=group)
            addSpring(button.group)
            gbutton("ok", handler=function(h, ...){
                chromosome <<- svalue(chr)
                start <<- as.numeric(svalue(xaxislim1))
                end <<- as.numeric(svalue(xaxislim2))
                xscaleold <<- c(start, end)
                strand <<- svalue(strand)
                if(length(filename)==1) {
                    filename <- filename[svalue(selectedfiles)]
                    selectedfilesFormat <- svalue(selectedfilesFormat)[svalue(selectedfiles)]
                }else {
                    filename <- filename[sapply(selectedfiles, svalue)]
                    selectedfilesFormat <- sapply(selectedfilesFormat, svalue)
                    selectedfilesFormat <- selectedfilesFormat[sapply(selectedfiles, svalue)]
                }
                gr <- GRanges(chromosome, IRanges(1, 249250621), strand)
 #               pb <- tkProgressBar("import progress bar", "imported", 
 #                                   max=length(filename), width=300)
                
                data <- mapply(function(.ele, format, i){
                    if(format=="BAM"){
                        importBam(.ele, ranges=gr)
                    }else{
                        .dat <- importScore(.ele, format=format)
                        if(format!="WIG"){
                            .dat@dat <- coverageGR(.dat@dat)
                        }
                        .dat
                    }
#                    setTkProgressBar(pb, i)
                }, filename, selectedfilesFormat, 1:length(filename))
   #             close(pb)
                require(svalue(genome), character.only=TRUE)
                txdb <<- get(svalue(genome))
                txdump <<- as.list(txdb)
                filename <- gsub("\\..*?$", "", filename)
                while(any(grepl("\\/|\\\\", filename))){
                    filename <- gsub("^.*(\\/|\\\\)", "", filename)
                }
                names(data) <- make.names(filename)
                trackList <<- trackList(data)
                autoOptimizeStyle <<- TRUE
                dispose(window)
            }, container=button.group)
            gbutton("cancel", handler = function(h,...) {
                dispose(window)
                stopped <<- TRUE
                },
                    container=button.group)
            return()
        }
        stopped <- FALSE
        trackList <- NULL
        if(missing(chromosome)) chromosome <- ""
        if(missing(start)) start <- 122929275
        if(missing(end)) end <- 122930122 ##249250621
        if(missing(strand)) strand <- "*"
        importData()
        while(is.null(trackList)){
            if(stopped) return()
            Sys.sleep(1)
        }
    }
    
    if(.Platform$OS.type=="unix"){
        if(!capabilities("X11")) stop("X11 is not available.")
        else{
            if(names(dev.cur())!="X11"){
                x11(type="Xlib")
            }
        }
    }
    
    if(class(trackList)!="trackList" || 
           !((is.list(trackList) && all(sapply(trackList, class)=="track")))){
        stop("trackList must be an object of \"trackList\"
             (See ?trackList) or a list of track")
    }
    if(missing(chromosome) || missing(start) || missing(end))
        stop("Please input the coordinate.")
    if(missing(strand) || !strand %in% c("+", "-"))
        strand <- "*"
    if(class(viewerStyle)!="trackViewerStyle"){
        stop("viewerStyle must be an object of 'trackViewerStyle'.")
    }
    if(autoOptimizeStyle){
        opt <- optimizeStyle(trackList, viewerStyle)
        trackList <- opt$tracks
        viewerStyle <- opt$style
    }
    yHeights <- NULL
    yscales <- NULL
    trackBrowser <- function(trackList, chromosome, start, end, strand, 
                             viewerStyle){
        grid.newpage()
        xscale <- c(start, end)
        margin <- viewerStyle@margin
        if(!is.null(txdb)){
            trs <- trackList[sapply(trackList, function(.ele) .ele@type=="gene")]
            if(length(trs)>0){
                trsRange <- lapply(trs, function(.ele) range(.ele@dat))
                trsRange <- range(unlist(GRangesList(trsRange)))
            }else{
                trsRange <- GRanges("unknown", IRanges(1,2))
            }
            if(as.character(seqnames(trsRange))[1]!=chromosome 
               || start < min(c(start(trsRange), xscaleold)) 
               || end > max(c(end(trsRange), xscaleold))){
                trackList <- trackList[sapply(trackList, function(.ele) .ele@type=="data")]
                if(is.null(txdump)) txdump <<- as.list(txdb)
                trs <- tryCatch(geneModelFromTxdb(txdb, chromosome, start, end, "*", txdump),
                                error=function(e) return(NULL))
                if(length(trs)>0){
                    trackList <<- 
                        trackList <- trackList(trackList, trs, 
                                               heightDist=c(length(trackList)/(length(trackList)+1), 
                                                            1/(length(trackList)+1)))
                }
            }
        }
        xscaleold <- xscale
        trackList <- filterTracks(trackList, chromosome, start, end, strand)
        yscales <<- getYlim(trackList)
        yHeights <<- getYheight(trackList)
        pushViewport(viewport(x=margin[2], y=margin[1], width=1-margin[2]-margin[4],
                              height=1-margin[1]-margin[3], just=c(0,0)))
        if(viewerStyle@xaxis){##draw x axis in bottom
            drawXaxis(xscale, viewerStyle)
        }
        popViewport()
        pushViewport(viewport(x=0, y=margin[1], width=1,
                              height=1-margin[1]-margin[3], just=c(0,0)))
        total <- length(trackList)
        ht <- 0
        for(i in 1:total){
            trackList[[i]]@style@height <<- yHeights[i]
            plotTrack(names(trackList)[i], trackList[[i]], 
                      viewerStyle, ht,
                      yscales[[i]], yHeights[i], xscale,
                      chromosome, strand)
            ht <- ht + yHeights[i]
        }
        popViewport()
    }
#    args <- as.list(match.call())[-1]
#    vp <- do.call(viewTracks, args)
#    args$start <- NULL
#    args$end <- NULL
    trackBrowser(trackList, chromosome, start, end, strand, viewerStyle)
    
    checkCursorY <- function(y){
        y <- grconvertX(y, "ndc", "npc")
        if(y < viewerStyle@margin[1]){
            return(0)
        }else{
            if(y > 1 - viewerStyle@margin[3]){
                return(length(yHeights)+1)
            }else{
                h <- cumsum(yHeights)
                h <- rescale(h, to=c(viewerStyle@margin[1], 
                                     1 - viewerStyle@margin[3]),
                             from=c(0, 1))
                return(which(y<=h)[1])
            }
        }                
    }
    checkCursorX <- function(x, y, n){
        h <- cumsum(yHeights)
        h <- rescale(h, to=c(viewerStyle@margin[1], 
                             1 - viewerStyle@margin[3]),
                     from=c(0, 1))
        h0 <- ifelse(n>1, h[n-1], viewerStyle@margin[1])
        top <- h[n] - yHeights[n]*trackList[[n]]@style@marginTop
        bottom <- h0 + yHeights[n]*trackList[[n]]@style@marginBottom
        if(y>top || y<bottom){
            return("panel")
        }
        x <- grconvertX(x, "ndc", "npc")
        left <- viewerStyle@margin[2]
        right <- 1- viewerStyle@margin[4]
        if(x>left && x<right){
            return("track")
        }
        return("panel")
    }
    saveFile <- function(h, ...){
        window <- gwindow("save figure", width=300, height=300)
        group <- ggroup(container=window, horizontal=FALSE)
        info.group <- ggroup(horizontal=FALSE, container=group)
        tmp <- ggroup(horizontal=TRUE, container=info.group)
        flab <- glabel("filename: ", container =tmp)
        filename <- gfilebrowse("trackViewer.pdf", type="save", 
                                initialfilename = "trackViewerOutput.pdf",
                                container=tmp, width=25)
        tmp <- ggroup(horizontal=TRUE, container=info.group)
        wlab <- glabel("width: ", container =tmp)
        width <- gedit(dev.size()[1], container=tmp, width=25)
        tmp <- ggroup(horizontal=TRUE, container=info.group)
        hlab <- glabel("height: ", container =tmp)
        height <- gedit(dev.size()[2], container=tmp, width=25)
        button.group <- ggroup(horizontal=TRUE, container=group)
        addSpring(button.group)
        gbutton("ok", handler=function(h, ...){
            dispose(window)
            pdf(svalue(filename), 
                width=as.numeric(svalue(width)), 
                height=as.numeric(svalue(height)))
            trackBrowser(trackList, chromosome, start, end, strand, viewerStyle)
            dev.off()
        }, container=button.group)
        gbutton("cancel", handler = function(h,...) dispose(window),
                container=button.group)
        return()
    }
    xaxisDialog <- function(){
        updateData <- function(h, ...){
            viewerStyle@xaxis <<- svalue(xaxisdraw)=="Yes"
            viewerStyle@xlas <<- ifelse(svalue(xaxislas)=="horizontal", 0, 2)
            viewerStyle@xgp$cex <<- as.numeric(svalue(xaxiscex))
            viewerStyle@xgp$col <<- svalue(xaxiscolor)
            chromosome <<- svalue(chr)
            start <<- as.numeric(svalue(xaxislim1))
            end <<- as.numeric(svalue(xaxislim2))
            strand <<- svalue(strand)
            viewerStyle@margin[2] <<- as.numeric(svalue(mleft))
            viewerStyle@margin[4] <<- as.numeric(svalue(mright))
            viewerStyle@margin[1] <<- as.numeric(svalue(mbottom))
            viewerStyle@margin[3] <<- as.numeric(svalue(mtop))
        }
        window <- gwindow("x-axis dialog", width=300, height=300)
        group <- ggroup(horizontal=FALSE, container=window)
        ##xaxis font size
        ##xaxis lab direction
        ##xaxis height
        ##xaxis font color
        tmp <- gframe("x axis", container=group, expand=TRUE, horizontal=FALSE) 
        ltmp <- ggroup(horizontal=TRUE, container=tmp)
        yaxis <- glabel("draw? ", container =ltmp)
        xaxisdraw <- gradio(c("Yes", "No"), horizontal = TRUE, 
                            selected=ifelse(viewerStyle@xaxis, 1, 2),
                            container=ltmp)
        ltmp <- ggroup(horizontal=TRUE, container=tmp)
        xaxis <- glabel("label direction: ", container =ltmp)
        xaxislas <- gradio(c("horizontal", "vertical"), horizontal = TRUE,
                            selected=ifelse(viewerStyle@xlas %in% c(0, 1), 1, 2),
                            container=ltmp)
        ltmp <- ggroup(horizontal=TRUE, container=tmp)
        xaxis <- glabel("cex: ", container =ltmp)
        xaxiscex <- gedit(ifelse(is.null(viewerStyle@xgp$cex), 
                                 optFontSize("x", viewerStyle), 
                                 viewerStyle@xgp$cex), 
                          container=ltmp, width=25)
        ltmp <- ggroup(horizontal=TRUE, container=tmp)
        xaxis <- glabel("color: ", container =ltmp)
        xaxiscolor <- gedit(ifelse(is.null(viewerStyle@xgp$col), "black", 
                                   viewerStyle@xgp$col), 
                            container=ltmp, width=8)
        tmp <- gframe("genomic coordinates", container=group, expand=TRUE, horizontal=FALSE) 
        ltmp <- ggroup(horizontal=TRUE, container=tmp)
        ltmpl <- ggroup(horizontal=FALSE, container=ltmp)
        label <- glabel("chromosome", container =ltmpl)
        label <- glabel("right: ", container =ltmpl)
        label <- glabel("top: ", container =ltmpl)
        label <- glabel("bottom: ", container =ltmpl)
        ltmpr <- ggroup(horizontal=FALSE, container=ltmp)
        chr <- gedit(chromosome, container=ltmpr, width=8)
        xaxislim1 <- gedit(start, container=ltmpr, width=16)
        xaxislim2 <- gedit(end, container=ltmpr, width=16)
        availStrand <- c("+", "-", "*")
        strand <- gradio(availStrand, selected=which(availStrand==strand),
                         container=ltmpr, width=8, horizontal=TRUE)
        
        tmp <- gframe("viewer margin", container=group, expand=TRUE, horizontal=FALSE) 
        ltmp <- ggroup(horizontal=TRUE, container=tmp)
        ltmpl <- ggroup(horizontal=FALSE, container=ltmp)
        margin <- glabel("left: ", container =ltmpl)
        margin <- glabel("right: ", container =ltmpl)
        margin <- glabel("top: ", container =ltmpl)
        margin <- glabel("bottom: ", container =ltmpl)
        ltmpr <- ggroup(horizontal=FALSE, container=ltmp)
        mleft <- gedit(viewerStyle@margin[2], container=ltmpr, width=8)
        mright <- gedit(viewerStyle@margin[4], container=ltmpr, width=8)
        mtop <- gedit(viewerStyle@margin[3], container=ltmpr, width=8)
        mbottom <- gedit(viewerStyle@margin[1], container=ltmpr, width=8)
        
        button.group <- ggroup(horizontal=TRUE, container=group)
        addSpring(button.group)
        gbutton("update", handler=function(h, ...){
            updateData()
            trackBrowser(trackList, chromosome, start, end, strand, viewerStyle)
        }, container=button.group)
        gbutton("return", handler=function(h, ...){
            dispose(window)
            eval.parent(getGraphicsEvent(), n=2)
        }, container=button.group)
        gbutton("save", handler=saveFile, container=button.group)
        
        return()
    }
    margintopDialog <- function(){
        updateData <- function(h, ...){
            viewerStyle@margin[2] <<- as.numeric(svalue(mleft))
            viewerStyle@margin[4] <<- as.numeric(svalue(mright))
            viewerStyle@margin[1] <<- as.numeric(svalue(mbottom))
            viewerStyle@margin[3] <<- as.numeric(svalue(mtop))
        }
        
        window <- gwindow("margin dialog", width=300, height=300)
        group <- ggroup(horizontal=FALSE, container=window)
        ##margin
        tmp <- gframe("viewer margin", container=group, expand=TRUE, horizontal=FALSE) 
        ltmp <- ggroup(horizontal=TRUE, container=tmp)
        ltmpl <- ggroup(horizontal=FALSE, container=ltmp)
        margin <- glabel("left: ", container =ltmpl)
        margin <- glabel("right: ", container =ltmpl)
        margin <- glabel("top: ", container =ltmpl)
        margin <- glabel("bottom: ", container =ltmpl)
        ltmpr <- ggroup(horizontal=FALSE, container=ltmp)
        mleft <- gedit(viewerStyle@margin[2], container=ltmpr, width=8)
        mright <- gedit(viewerStyle@margin[4], container=ltmpr, width=8)
        mtop <- gedit(viewerStyle@margin[3], container=ltmpr, width=8)
        mbottom <- gedit(viewerStyle@margin[1], container=ltmpr, width=8)
        
        button.group <- ggroup(horizontal=TRUE, container=group)
        addSpring(button.group)
        gbutton("update", handler=function(h, ...){
            updateData()
            trackBrowser(trackList, chromosome, start, end, strand, viewerStyle)
        }, container=button.group)
        gbutton("return", handler=function(h, ...){
            dispose(window)
            eval.parent(getGraphicsEvent(), n=2)
        }, container=button.group)
        gbutton("save", handler=saveFile, container=button.group)
        
        return()
    }
    trackLabelDialog <- function(n){
        updateData <- function(h, ...){
            trackList[[n]]@style@yaxis@draw <<- svalue(yaxisdraw)=="Yes"
            trackList[[n]]@style@yaxis@main <<- svalue(yaxismain)=="left"
            trackList[[n]]@style@yaxis@label <<- svalue(yaxislabel)=="Yes"
            trackList[[n]]@style@yaxis@gp$cex <<- as.numeric(svalue(yaxiscex))
            trackList[[n]]@style@yaxis@gp$col <<- svalue(yaxiscolor)
            trackList[[n]]@style@ylim <<- as.numeric(c(svalue(yaxislim1), svalue(yaxislim2)))
            if(as.numeric(svalue(height))>0 && 
                   as.numeric(svalue(height))<trackList[[n]]@style@height){
                trackList[[n]]@style@height <<- as.numeric(svalue(height))
            }
            viewerStyle@margin[2] <<- as.numeric(svalue(mleft))
            viewerStyle@margin[4] <<- as.numeric(svalue(mright))
            trackList[[n]]@style@marginTop <<- as.numeric(svalue(mtop))
            trackList[[n]]@style@marginBottom <<- as.numeric(svalue(mbottom))
            trackList[[n]]@style@ylablas <<- ifelse(svalue(ylabellas)=="vertical", 0, 1)
            trackList[[n]]@style@ylabgp$cex <<- as.numeric(svalue(ylabelcex))
            trackList[[n]]@style@ylabpos <<- svalue(ylabelpos)
            trackList[[n]]@style@ylabgp$col <<- svalue(ylabelcolor)
            names(trackList)[n] <<- svalue(ylabelname)
        }
        window <- gwindow("track label dialog", width=300, height=300)
        group <- ggroup(horizontal=FALSE, container=window)
        info.group <- ggroup(hroizontal=FALSE, container=group)
        ##yaxis font size
        ##yaxis position
        ##ylim
        tmp <- gframe("y axis", container=group, expand=TRUE, horizontal=FALSE) 
        ltmp <- ggroup(horizontal=TRUE, container=tmp)
        yaxis <- glabel("draw? ", container =ltmp)
        yaxisdraw <- gradio(c("Yes", "No"), horizontal = TRUE, 
                         selected=ifelse(trackList[[n]]@style@yaxis@draw, 1, 2),
                         container=ltmp)
        ltmp <- ggroup(horizontal=TRUE, container=tmp)
        yaxis <- glabel("position: ", container =ltmp)
        yaxismain <- gradio(c("left", "right"), horizontal = TRUE,
                         selected=ifelse(trackList[[n]]@style@yaxis@main, 1, 2),
                         container=ltmp)
        ltmp <- ggroup(horizontal=TRUE, container=tmp)
        yaxis <- glabel("labels on the tick marks: ", container =ltmp)
        yaxislabel <- gradio(c("Yes", "No"), horizontal = TRUE,
                            selected=ifelse(trackList[[n]]@style@yaxis@label, 1, 2),
                            container=ltmp)
        ltmp <- ggroup(horizontal=TRUE, container=tmp)
        yaxis <- glabel("cex: ", container =ltmp)
        yaxiscex <- gedit(ifelse(is.null(trackList[[n]]@style@yaxis@gp$cex), 
                                 optFontSize("y", viewerStyle, trackList[[n]]@style@height), 
                             trackList[[n]]@style@yaxis@gp$cex), 
                          container=ltmp, width=25)
        ltmp <- ggroup(horizontal=TRUE, container=tmp)
        yaxis <- glabel("color: ", container =ltmp)
        yaxiscolor <- gedit(ifelse(is.null(trackList[[n]]@style@yaxis@gp$col), "black", 
                               trackList[[n]]@style@yaxis@gp$col), 
                            container=ltmp, width=8)
        ltmp <- ggroup(horizontal=TRUE, container=tmp)
        yaxis <- glabel("ylim: from ", container =ltmp)
        yaxislim1 <- gedit(yscales[[n]][1], 
                            container=ltmp, width=6)
        yaxis <- glabel("to ", container =ltmp)
        yaxislim2 <- gedit(yscales[[n]][2], 
                           container=ltmp, width=6)
        
        ##track height
        tmp <- gframe("track", container=group, expand=TRUE, horizontal=FALSE) 
        ltmp <- ggroup(horizontal=TRUE, container=tmp)
        yaxis <- glabel("track height: ", container =ltmp)
        height <- gedit(trackList[[n]]@style@height, 
                           container=ltmp, width=16)
        ##margin
        tmp <- gframe("track margin", container=group, expand=TRUE, horizontal=FALSE) 
        ltmp <- ggroup(horizontal=TRUE, container=tmp)
        ltmpl <- ggroup(horizontal=FALSE, container=ltmp)
        margin <- glabel("left: ", container =ltmpl)
        margin <- glabel("right: ", container =ltmpl)
        margin <- glabel("top: ", container =ltmpl)
        margin <- glabel("bottom: ", container =ltmpl)
        ltmpr <- ggroup(horizontal=FALSE, container=ltmp)
        mleft <- gedit(viewerStyle@margin[2], container=ltmpr, width=8)
        mright <- gedit(viewerStyle@margin[4], container=ltmpr, width=8)
        mtop <- gedit(trackList[[n]]@style@marginTop, container=ltmpr, width=8)
        mbottom <- gedit(trackList[[n]]@style@marginBottom, container=ltmpr, width=8)
        
        ##ylabel
        ##ylabel direction
        ##ylabel color
        ##ylabel font size
        ##ylabel position
        tmp <- gframe("y label", container=group, expand=TRUE, horizontal=FALSE) 
        availpos<- c("left", "right", "topleft", 
                         "bottomleft", "topright", "bottomright")
        ltmp <- ggroup(horizontal=TRUE, container=tmp)
        ylabel <- glabel("label direction: ", container =ltmp)
        ylabellas <- gradio(c("horizontal", "vertical"), horizontal = TRUE,
                            selected=ifelse(trackList[[n]]@style@ylablas %in% c(0, 3), 2, 1),
                            container=ltmp)
        ltmp <- ggroup(horizontal=TRUE, container=tmp)
        ylabel <- glabel("position: ", container =ltmp)
        ylabelpos <- gradio(availpos, horizontal = FALSE,
                            selected=which(availpos==trackList[[n]]@style@ylabpos),
                            container=ltmp)
        ltmp <- ggroup(horizontal=TRUE, container=tmp)
        ylabel <- glabel("cex: ", container =ltmp)
        ylabelcex <- gedit(ifelse(is.null(trackList[[n]]@style@ylabgp$cex), 
                                  optFontSize("z", viewerStyle, trackList[[n]]@style@height), 
                                 trackList[[n]]@style@ylabgp$cex), 
                          container=ltmp, width=25)
        ltmp <- ggroup(horizontal=TRUE, container=tmp)
        ylabel <- glabel("color: ", container =ltmp)
        ylabelcolor <- gedit(ifelse(is.null(trackList[[n]]@style@ylabgp$col), "black", 
                                   trackList[[n]]@style@ylabgp$col), 
                            container=ltmp, width=8)
        ltmp <- ggroup(horizontal=TRUE, container=tmp)
        ylabel <- glabel("label: ", container =ltmp)
        ylabelname <- gedit(ifelse(is.null(names(trackList[n])), "NA", 
                                    names(trackList)[n]), 
                             container=ltmp, width=25)
        
        
        button.group <- ggroup(horizontal=TRUE, container=group)
        addSpring(button.group)
        gbutton("update", handler=function(h, ...){
            updateData()
            trackBrowser(trackList, chromosome, start, end, strand, viewerStyle)
        }, container=button.group)
        gbutton("return", handler=function(h, ...){
            dispose(window)
            eval.parent(getGraphicsEvent(), n=2)
        }, container=button.group)
        gbutton("save", handler=saveFile, container=button.group)
        return()
    }
    trackDialog<- function(n){
        updateData <- function(h, ...){
            if(svalue(zoomOp)=="zoom out"){
                interval <- round((svalue(zoomX) - 1) * interval / 2)
                start <<- start - interval
                end <<- end + interval
                if(start < 0){
                    end <<- end - start
                    start <<- 0
                }
                interval <<- abs(end - start)
            }else{
                interval <- round((svalue(zoomX) - 1)/svalue(zoomX) * interval / 2)
                start <<- start + interval
                end <<- end - interval
                interval <<- abs(end - start)
            }
            trackList[[n]]@style@color[1] <<- svalue(color1)
            trackList[[n]]@style@color[2] <<- svalue(color2)
        }
        window <- gwindow("track dialog", width=300, height=300)
        group <- ggroup(horizontal=FALSE, container=window)
        zoom <- c("zoom in", "zoom out")
        tmp <- gframe("Zoom", container=group, expand=TRUE, horizontal=FALSE)
        zoomOp <- gradio(zoom, horizontal=TRUE,
                         container=tmp)
        tmp <- gframe("Zoom times", container=tmp, expand=TRUE)
        zoomX <- gradio(c(1, 1.5, 3, 10), container=tmp)
        
        tmp <- gframe("color", container=group, expand=TRUE, horizontal=FALSE)
        ltmp <- ggroup(hroizontal=TRUE, container=tmp)
        clabel1 <- glabel("channel 1: ", container =ltmp)
        color1 <- gedit(trackList[[n]]@style@color[1], container=ltmp, width=8)
        ltmp <- ggroup(hroizontal=TRUE, container=tmp)
        clabel1 <- glabel("channel 2: ", container =ltmp)
        color2 <- gedit(trackList[[n]]@style@color[2], container=ltmp, width=8)
        
        button.group <- ggroup(horizontal=TRUE, container = group)
        addSpring(button.group)
        gbutton("update", handler=function(h, ...){
            updateData()
            trackBrowser(trackList, chromosome, start, end, strand, viewerStyle)
        }, container=button.group)
        gbutton("return", handler=function(h, ...){
            dispose(window)
            eval.parent(getGraphicsEvent(), n=2)
        }, container=button.group)
        gbutton("save", handler=saveFile, container=button.group)
        
        return()
    }
    
    interval <- abs(end - start)
    startx <- NULL
    starty <- NULL
    usr <- NULL
    
    devset <- function()
        if (dev.cur() != eventEnv$which) dev.set(eventEnv$which)
    
    dragmousedown <- function(buttons, x, y) {
        startx <<- x
        starty <<- y
        devset()
        usr <<- par("usr")
#        eventEnv$onMouseMove <- dragmousemove
        eventEnv$onMouseUp <- mouseup
        NULL
    }
    
    dragmousemove <- function(buttons, x, y) {
        NULL
    }
    
    mouseup <- function(buttons, x, y) {
        eventEnv$onMouseUp <- NULL
        devset()
        deltax <- diff(grconvertX(c(startx, x), "ndc", "npc"))
#        deltay <- diff(grconvertY(c(starty, y), "ndc", "npc"))
        if(abs(deltax)>0.005){
            delta <- round(interval*deltax)
            start <- start - delta
            end <- end - delta
            if(start<0){
                end <- end - start
                start <- 0
            }
            #do.call(viewTracks, c(args, start, end))
            start <<- start
            end <<- end
        }else{
            ##check x, y position
            ##check y to determine the track number or margin
            n <- checkCursorY(y)
            if(n==0){
                xaxisDialog()
            }else{
                if(n>length(trackList)){
                    margintopDialog()
                }else{
                    ##check x to determine the target region
                    if(checkCursorX(x, y, n)=="panel")
                        trackLabelDialog(n)
                    else trackDialog(n)
                }
            }
            return(invisible(1))
        }
        trackBrowser(trackList, chromosome, start, end, strand, viewerStyle)
#        eventEnv$onMouseMove <- NULL
        NULL
    }	
    
    keydown <- function(key) {
        if (key == "q") {
            ##save
            dev.off()
            return(invisible(1))
        }
 #       eventEnv$onMouseMove <- NULL
        NULL
    }
    
    setGraphicsEventHandlers(prompt = "Click and drag, hit q to quit",
                             onMouseDown = dragmousedown,
                             onKeybd = keydown)
    eventEnv <- getGraphicsEventEnv()
    
    devset()
    
    eval.parent(getGraphicsEvent(), n=2)
}