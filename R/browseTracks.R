#' browse tracks
#'
#' @description browse tracks by a web browser.
#' @param trackList an object of \code{\link{trackList}}
#' @param gr an object of \code{\link[GenomicRanges:GRanges-class]{GRanges}}
#' @param ignore.strand ignore the strand or not when do filter. default TRUE
#' @param width width of the figure
#' @param height height of the figure
#' @param ... parameters not used
#' @import htmlwidgets
#' @import GenomicRanges
#' @importFrom utils installed.packages
#' @export
#' @return An object of class htmlwidget that will intelligently print itself 
#' into HTML in a variety of contexts including the R console, 
#' within R Markdown documents, and within Shiny output bindings.
#' @examples 
#' extdata <- system.file("extdata", package="trackViewer", mustWork=TRUE)
#' files <- dir(extdata, "-.wig")
#' tracks <- lapply(paste(extdata, files, sep="/"), 
#'                  importScore, format="WIG")
#' tracks <- lapply(tracks, function(.ele) {strand(.ele@dat) <- "-"; .ele})
#' names(tracks) <- c("trackA", "trackB")
#' fox2 <- importScore(paste(extdata, "fox2.bed", sep="/"), format="BED")
#' dat <- coverageGR(fox2@dat)
#' fox2@dat <- dat[strand(dat)=="+"]
#' fox2@dat2 <- dat[strand(dat)=="-"]
#' gr <- GRanges("chr11", IRanges(122929275, 122930122))
#' browseTracks(trackList(tracks, fox2), gr=gr)

browseTracks <- function(trackList, 
                         gr=GRanges(),
                         ignore.strand=TRUE,
                         width=NULL, height=NULL, 
                         ...){
    if(missing(gr)){
      if(interactive()){
        ## ask to input the region
        chr <- readline(prompt = "Enter the chromosome, for example chr1: ")
        start <- as.integer(readline(prompt = "Enter the start coordinate: "))
        end <- as.integer(readline(prompt = "Enter the end coordinate: "))
        gr <- GRanges(chr, IRanges(start, end))
      }else{
        stop("gr is required.")
      }
    }else{
      if(!is(gr, "GRanges")){
        stop("gr must be an object of GRanges.")
      }
      if(length(gr)!=1){
        stop("the length of gr must be 1.")
      }
    }
    if(missing(trackList)){
      if(interactive()){
        ## ask to import files
        trackList <- list()
        n <- 1
        while(n>0){
          message("Please select \n\t0\t stop import.\n\t1\t get gene track. \n\t2\t import BED/bedGraph/WIG/BigWig file.")
          filename <- readline(prompt = "Enter the number: ")
          switch(filename, 
                 "0"={
                   n <- 0
                 },
                 "1"={
                   ip <- as.data.frame(installed.packages()[, c(1, 3:4)])
                   ip <- unique(as.character(ip$Package))
                   txdbs <- ip[grepl("^TxDb", ip)]
                   orgs <- ip[grepl("^org", ip)]
                   txdbFlag <- TRUE
                   message("available txdb\n", paste0(paste0(seq_along(txdbs), ": ", txdbs, "\n")))
                   while(txdbFlag){
                     txdb <- readline(prompt = "select txdb: ")
                     txdb <- as.numeric(txdb)
                     if(!is.na(txdb)){
                       txdb <- txdbs[txdb]
                       txdbFlag <- FALSE
                     }
                   }
                   orgFlag <- TRUE
                   message("available org\n", paste0(paste0(seq_along(orgs), ": ", orgs, "\n")))
                   while(orgFlag){
                     org <- readline(prompt = "select org: ")
                     org <- as.numeric(org)
                     if(!is.na(org)){
                       org <- orgs[org]
                       orgFlag <- FALSE
                     }
                   }
                   require(txdb, character.only = TRUE)
                   require(org, character.only = TRUE)
                   trs <- geneModelFromTxdb(get(txdb), get(org), gr=gr)
                   trackList <- c(trackList, trs)
                 },
                 "2"={
                   filename <- file.choose()
                   if(file.exists(filename)){
                     format <- sub("^.*\\.(.*)$", "\\1", basename(filename))
                     if(tolower(format) %in% tolower(c("BED", "bedGraph", "WIG", "BigWig"))){
                       format <- switch(tolower(format), 
                                        bed="BED",
                                        bedgraph="bedGraph",
                                        wig="WIG", 
                                        bigwig="BigWig")
                       sampleName <- make.names(sub(paste0(".", format), "", basename(filename)))
                       if(sampleName %in% names(trackList)){
                         sampleName <- paste0(sampleName, "_copy")
                       }
                       trackList[[sampleName]] <- importScore(filename, format = format, ranges = gr)
                     }else{
                       message("Can not detect the format of the file. Please rename the file with end of .BED, .bedGraph, .WIG or .BigWig")
                     }
                   }else{
                     message("Can not open the file. Please double check the file path.")
                   }
                 })
        }
        trackList <- optimizeStyle(trackList(trackList), theme="col")$tracks
      }else{
        stop("trackList is required.")
      }
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
    chromosome <- as.character(GenomicRanges::seqnames(gr))[1]
    start <- GenomicRanges::start(gr)[1]
    end <- GenomicRanges::end(gr)[1]
    if(ignore.strand) GenomicRanges::strand(gr) <- "*"
    strand <- as.character(GenomicRanges::strand(gr))[1]
    if(end < start) stop("end should be greater than start.")
    
    trackList <- filterTracks(trackList, chromosome, 
                              start, end, strand)
    trackList <- rev(trackList)
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
        if(is(x, "AsIs")){
          return(lapply(x, col2Hex))
        }
        rgb <- col2rgb(x)
        rgb(rgb[1, ], rgb[2, ], rgb[3, ], maxColorValue = 255)
    }
    LINEW <- as.numeric(convertX(unit(1, "line"), "npc"))
    getLolliplotData <- function(.dat){
      if(length(.dat)==0) return(ZERO);
      .dat <- sort(.dat)
      .dat$textlabel <- names(.dat)
      names(.dat) <- NULL
      if(length(.dat$score)!=length(.dat)){
        .dat$score <- 1
      }
      
      cex <- .5
      if(length(.dat$cex)>0) cex <- cex*sapply(cex, `[`, 1)
      pushViewport(viewport(xscale=c(start(gr), end(gr))))
      
      if(length(.dat$stack.factor)>0){
        stopifnot(is.vector(.dat$stack.factor, mode="character"))
        if(length(.dat$stack.factor.order)>0 || 
           length(.dat$stack.factor.first)>0){
          warning("stack.factor.order and stack.factor.first are used by this function!",
                  "The values in these column will be removed.")
        }
        ## condense the SNPs
        stack.factors <- unique(as.character(.dat$stack.factor))
        stack.factors <- sort(stack.factors)
        stack.factors.order <- 1:length(stack.factors)
        names(stack.factors.order) <- stack.factors
        .dat <- .dat[order(as.character(seqnames(.dat)), start(.dat), 
                           as.character(.dat$stack.factor))]
        .dat$stack.factor.order <- stack.factors.order[.dat$stack.factor]
        .dat$stack.factor.first <- !duplicated(.dat)
        .dat.condense <- .dat
        .dat.condense$oid <- 1:length(.dat)
        .dat.condense$factor <- paste(as.character(seqnames(.dat)), start(.dat), end(.dat))
        .dat.condense <- split(.dat.condense, .dat.condense$factor)
        .dat.condense <- lapply(.dat.condense, function(.ele){
          .oid <- .ele$oid
          .gr <- .ele[1]
          mcols(.gr) <- NULL
          .gr$oid <- NumericList(.oid)
          .gr
        })
        .dat.condense <- unlist(GRangesList(.dat.condense), use.names = FALSE)
        .dat.condense <- sort(.dat.condense)
        lab.pos.condense <- jitterLables(start(.dat.condense), 
                                         xscale=c(start(gr), end(gr)), 
                                         lineW=LINEW*cex)
        lab.pos.condense <- reAdjustLabels(lab.pos.condense, 
                                           lineW=LINEW*cex)
        condense.ids <- .dat.condense$oid
        lab.pos <- rep(lab.pos.condense, elementNROWS(condense.ids))
        lab.pos <- lab.pos[order(unlist(condense.ids))]
      }else{
        lab.pos <- jitterLables(start(.dat), 
                                xscale=c(start(gr), end(gr)), 
                                lineW=LINEW*cex)
        lab.pos <- reAdjustLabels(lab.pos, 
                                  lineW=LINEW*cex)
      }
      popViewport()
      
      dat2 <- as.list(as.data.frame(.dat))
      dat2$seqnames <- chromosome
      dat2$border <- col2Hex(dat2$border)
      dat2$color <- col2Hex(dat2$color)
      dat2$labpos <- lab.pos
      dat2
    }
    trackdat <- function(.ele){
        dat <- .ele$dat
        dat2 <- .ele$dat2
        style <- object2list(.ele$style)
        style$color <- col2Hex(style$color)
        if(.ele$type=="data"){
            dat <- getData(.ele$dat)
            dat2 <- getData(.ele$dat2)
            if(length(.ele$style@ylim)==2) {
                ylim <- .ele$style@ylim
              }else{
                ylim <- range(c(dat, -1*dat2))
              }
            if(length(style$color)==1){
                style$color <- rep(style$color, 2)
            }
        }else{
          if(.ele$type=="lollipopData"){
            dat <- getLolliplotData(.ele$dat)
            if(length(.ele$dat2)==0){
              dat2 <- ZERO
            }else{
              dat2 <- getLolliplotData(.ele$dat2)
            }
            ylim <- c(0, 1)
            if(length(style$color)==1){
              style$color <- rep(style$color, 2)
            }
          }else{
            if(.ele$type=="scSeq"){
              stop("Not supported.")
            }else{
              ##gene
              dat$textlabel <- names(.ele$dat)
              names(dat) <- NULL
              if(length(dat$featureLayerID)==0){
                dat$featureLayerID <- 1
              }
              if(length(dat$fill)==0){
                dat$fill <- style$color[1]
              }else{
                dat$fill <- col2Hex(dat$fill)
              }
              dat <- as.list(as.data.frame(dat))
              dat$seqnames <- chromosome
              if(length(.ele$dat2)==0){
                dat2 <- ZERO
              }else{
                dat2 <- getLolliplotData(.ele$dat2)
              }
              ylim <- c(0, 1)
            }
          }
        }
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


#' Shiny bindings for browseTracks
#'
#' Output and render functions for using browseTracks within Shiny
#' applications and interactive Rmd documents.
#'
#' @param outputId output variable to read from
#' @param width,height Must be a valid CSS unit (like \code{'100\%'},
#'   \code{'400px'}, \code{'auto'}) or a number, which will be coerced to a
#'   string and have \code{'px'} appended.
#' @param expr An expression that generates a browseTracks
#' @param env The environment in which to evaluate \code{expr}.
#' @param quoted Is \code{expr} a quoted expression (with \code{quote()})? This
#'   is useful if you want to save an expression in a variable.
#'
#' @name browseTracks-shiny
#'
#' @export
browseTracksOutput <- function(outputId, width = '100%', height = '600px'){
  htmlwidgets::shinyWidgetOutput(outputId, 'browseTracks', width, height, 
                                 package = 'trackViewer')
}

#' @rdname browseTracks-shiny
#' @export
renderbrowseTracks <- function(expr, env = parent.frame(), quoted = FALSE) {
  if (!quoted) { expr <- substitute(expr) } # force quoted
  htmlwidgets::shinyRenderWidget(expr, browseTracksOutput, env, quoted = TRUE)
}

