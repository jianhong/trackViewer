#' Reading data from a BED or WIG file
#' @description Read a \code{\link{track}} object from a BED, bedGraph, WIG or BigWig file
#' @param file The path to the file to read.
#' @param file2 The path to the second file to read.
#' @param format The format of import file. Could be BED, bedGraph, WIG or BigWig
#' @param ranges An object of \code{\link[GenomicRanges:GRanges-class]{GRanges}} to indicate
#' the range to be imported
#' @param ignore.strand ignore the strand or not when do filter. default TRUE
#' @return a \code{\link{track}} object 
#' @import GenomicRanges
#' @import IRanges
#' @export
#' @seealso See Also as \code{\link{importBam}}, \code{\link{track}}, 
#' \code{\link{viewTracks}}
#' @examples 
#' #import a BED file
#' bedfile <- system.file("tests", "test.bed", package="rtracklayer",
#'                        mustWork=TRUE)
#' dat <- importScore(file=bedfile, format="BED",
#'                    ranges=GRanges("chr7", IRanges(127471197, 127474697)))
#' 
#' ##import a WIG file
#' wigfile <- system.file("tests", "step.wig", package = "rtracklayer",
#'                        mustWork=TRUE)
#' dat <- importScore(file=wigfile, format="WIG")
#' 
#' ##import a BigWig file
#' if(.Platform$OS.type!="windows"){##this is because we are using rtracklayer::import
#'   bwfile <- system.file("tests", "test.bw", package = "rtracklayer",
#'                         mustWork=TRUE)
#'   dat <- importScore(file=bwfile, format="BigWig")
#' }
#' 
#' ##import 2 file
#' wigfile1 <- system.file("extdata", "cpsf160.repA_+.wig", package="trackViewer",
#'                         mustWork=TRUE)
#' wigfile2 <- system.file("extdata", "cpsf160.repA_-.wig", package="trackViewer",
#'                         mustWork=TRUE)
#' dat <- importScore(wigfile1, wigfile2, format="WIG", 
#'                    ranges=GRanges("chr11", IRanges(122817703, 122889073)))

importScore <- function(file, file2, 
                        format=c("BED", "bedGraph", "WIG", "BigWig"),
                        ranges=GRanges(), ignore.strand=TRUE){
    if(missing(file))
        stop("file is required.")
    ##    on.exit(closeAllConnections())
    format <- match.arg(format)
    #    res <- GRanges(score=numeric(0))
    if(!is(ranges, "GRanges")) stop("ranges must be an object of GRanges.")
    gr <- orderedGR(ranges)
    seqn <- unique(as.character(seqnames(gr)))
    filterByRange <- function(r){
        if(length(gr)>0){
            ## find in ranges
            r <- r[r[, 1] %in% seqn, , drop=FALSE]
            nr <- nrow(r)
            if(nr>0){
                ## split r and findOverlaps
                idx <- rep(FALSE, nr)
                l <- floor(nr/1000)
                for(i in 0:l){
                    f <- min(i*1000+1, nr)
                    t <- min((i+1)*1000, nr)
                    x <- r[f:t, , drop=FALSE]
                    xgr <- GRanges(x[,1], IRanges(start=as.numeric(x[,2]),
                                                  end=as.numeric(x[,3])))
                    suppressWarnings(ol <- findOverlaps(xgr, gr, ignore.strand=ignore.strand))
                    if(length(ol)>0) idx[queryHits(ol)+i*1000] <- TRUE
                }
                r <- r[idx, , drop=FALSE]
            }
        }
        r
    }
    
    getWigInfo <- function(firstline){
        firstline <- unlist(strsplit(firstline, "\\s"))
        firstline <- firstline[firstline!=""]
        structure <- firstline[1]
        firstline <- firstline[-1]
        firstline <- do.call(rbind,strsplit(firstline, "=", fixed=TRUE))
        firstline <- firstline[match(c("chrom","span","start","step"), 
                                     firstline[,1]),]
        info <- c(structure, firstline[,2])
        names(info) <- c("structure", "chrom","span","start","step")
        return(info)
    }
    readWIG <- function(buf, lastWigInfo=NULL){
        ##remove browser lines and track lines
        buf <- gsub("^\\s+", "", buf)
        buf <- gsub("\\s+$", "", buf)
        buf <- buf[grepl("^(variableStep|fixedStep|([0-9]+))", buf)]
        infoLine <- grep("Step", buf)
        if(length(infoLine)>0){
            if(infoLine[1]!=1){
                if(is.null(lastWigInfo[1])){
                    ##1 line is the info line
                    stop("WIG file must contain track definition line, 
                     which should start by variableStep or fixedStep.")
                }else{
                    buf <- c(lastWigInfo, buf)
                    infoLine <- grep("Step", buf)
                }
            }
        }else{
            if(is.null(lastWigInfo[1])){
                stop("WIG file must contain track definition line, 
                     which should start by variableStep or fixedStep.")
            }else{
                buf <- c(lastWigInfo, buf)
                infoLine <- grep("Step", buf)
            }
        }
        
        lastWigInfo <- buf[infoLine[length(infoLine)]]
        
        while(infoLine[length(infoLine)]==length(buf)){
            buf <- buf[-length(buf)]
        }
        block <- c(infoLine, length(buf)+1)
        dif <- diff(block)
        block <- rep(infoLine, dif)
        buf <- split(buf, block)
        r <- lapply(buf, function(.ele) {
            wiginfo <- getWigInfo(.ele[1])
            span <- wiginfo["span"]
            step <- as.numeric(wiginfo["step"])
            if(wiginfo["structure"]=="variableStep"){
                start <- as.numeric(strsplit(.ele[2], "\\s+")[[1]][1])
                if(is.na(span)) span <- 1
                lastrow <- strsplit(.ele[length(.ele)], "\\s+")[[1]]
                end <- as.numeric(lastrow)[1] + as.numeric(span)
            }else{##wiginfo["structure"]=="fixedStep"
                start <- as.numeric(wiginfo["start"])
                if(!is.na(span)){
                    end <- start + (length(.ele) - 1) * step + as.numeric(span) - 1
                }else{
                    end <- start + length(.ele) * step - 1
                }
            }
            c(wiginfo["chrom"], start, end, span, step, wiginfo["structure"])
        })
        r <- do.call(rbind, r)
        wiginfo <- getWigInfo(lastWigInfo)
        if(wiginfo["structure"]=="fixedStep"){
            lastWigInfo <- gsub("start=\\d+(\\s)", 
                                paste("start=",
                                      as.numeric(r[nrow(r), 3])+1,
                                      "\\1", sep=""), lastWigInfo)
        }
        buf <- lapply(buf, "[", -1)
        buf <- CharacterList(buf, compress=TRUE)
        ## filter by range
        if(length(gr)>0){
            r <- cbind(r, rid=1:nrow(r))
            r <- filterByRange(r)
            buf <- buf[as.numeric(r[,"rid"])]
        }
        
        list(gr=GRanges(seqnames=r[,1], 
                ranges=IRanges(start=as.numeric(r[,2]),
                               end=as.numeric(r[,3])),
                score=buf,
                span=as.numeric(r[,4]),
                step=as.numeric(r[,5]),
                structure=r[,6]),
             lastWigInfo=lastWigInfo)
    }
    readBED <- function(buf){
        ##c("chrom", "chromStart", "chromEnd", "name", "score", "strand", ...)
        ##remove annotations
        buf <- strsplit(buf, "\t", fixed=TRUE)##
        len <- sapply(buf, length)
        buf <- buf[len>2]
        len <- len[len>2]
        if(length(buf)<1) return(GRanges(score=numeric(0)))
        maxLen <- max(len)
        if(all(len==maxLen)){
            buf <- do.call(rbind, buf)
        }else{
            NAs <- rep("", maxLen)
            buf <- do.call(rbind, 
                           lapply(buf, function(.ele) c(.ele, NAs)[1:maxLen]))
        }
        ##check strand
        if(ncol(buf)==3) buf <- cbind(buf, ".")
        if(ncol(buf)==4) {
            if(all(grepl("^[\\d\\.]+$", buf[,4])) && 
                   length(unique(nchar(buf[,4])))>1){
                buf <- cbind(buf, buf[, 4])
            }else{
                buf <- cbind(buf, 1)
            }
        }
        if(ncol(buf)==5) buf <- cbind(buf, "*")
        buf[!buf[,6] %in% c("+", "-"), 6] <- "*"
        ##filter by range
        buf <- filterByRange(buf)
        if(nrow(buf)>0){
            GRanges(seqnames=buf[,1], 
                    ranges=IRanges(start=as.numeric(buf[,2])+1,##0 based, half open
                                   end=as.numeric(buf[,3])),
                    strand=buf[,6],
                    score=as.numeric(buf[,5]))
        }else{##create a length=0 GRanges
            GRanges(seqnames=buf[,1], 
                    ranges=IRanges(start=as.numeric(buf[,2]),
                                   end=as.numeric(buf[,3])),
                    strand=buf[,6],
                    score=as.numeric(buf[,5]))
        }
    }
    readFourCols <- function(buf){
        buf <- gsub("^\\s+", "", buf)
        buf <- gsub("\\s+$", "", buf)
        buf <- buf[!grepl("^(browser|track|#)", buf)]
        buf <- strsplit(buf, "\t", fixed=TRUE)
        len <- sapply(buf, length)
        buf <- buf[len==4]
        if(length(buf)<1) return(GRanges(score=numeric(0)))
        buf <- do.call(rbind, buf)
        ##filter by range
        buf <- filterByRange(buf) 
        if(nrow(buf)>0){
            GRanges(seqnames=buf[,1], 
                    ranges=IRanges(start=as.numeric(buf[,2])+1,##0 based, half open
                                   end=as.numeric(buf[,3])),
                    score=as.numeric(buf[,4]))
        }else{##create a length=0 GRanges
            GRanges(seqnames=buf[,1], 
                    ranges=IRanges(start=as.numeric(buf[,2]),##0 based, half open
                                   end=as.numeric(buf[,3])),
                    score=as.numeric(buf[,4]))
        }
    }
    readbedGraph <- function(buf){
        readFourCols(buf)
    }
    readBigWig <- function(file){
        if(length(gr)>0){
            import(con=file, format="BigWig", which=gr)
        }else{
            import(con=file, format="BigWig")
        } 
    }
    readFile <- function(file, format, FUN){
        if(format=="WIG"){
            res <- NULL
            con <- file(file, open="r")
            on.exit(close(con))
            lastWigInfo <- NULL
            while(length(buf <- readLines(con, n=1000000, warn=FALSE))>0){
                buf <- FUN(buf, lastWigInfo)
                lastWigInfo <- buf$lastWigInfo
                if(length(res)<1) {
                    res <- buf$gr
                }else{
                    suppressWarnings(res <- c(res, buf$gr))
                }
            }
        }else{
            s <- file.info(file)$size
            if(s<100000000){
                buf <- readChar(file, s, useBytes=TRUE)
                buf <- strsplit(buf, "\n", fixed=TRUE, useBytes=TRUE)[[1]]
                res <- FUN(buf)
            }else{
                message("file is too huge. Please consider to use bedtools or bedops to subset the data.")
                res <- NULL
                con <- file(file, open="r")
                on.exit(close(con))
                while(length(buf <- readLines(con, n=1000000, warn=FALSE))>0){
                    buf <- FUN(buf)
                    if(length(res)<1) {
                        res <- buf
                    }else{
                        suppressWarnings(res <- c(res, buf))
                    }
                }
            }
        }
        res <- unique(res)
        return(res)
    }
    readFiles <- function(file, format){
        FUN <- get(paste("read", format, sep=""))
        if(format=="BigWig"){
            if(.Platform$OS.type=="windows")
                stop("rtracklayer can't read bigWig files on a Windows computer. Type in  ?`BigWigFile-class` to get the help.")
            res <- unique(FUN(file))
        }else{
            res <- readFile(file, format, FUN)
        }
        return(res)
    }
    res <- readFiles(file, format)
    if(!missing(file2)){
        res2 <- readFiles(file2, format)
        return(new("track", dat=orderedGR(res), dat2=orderedGR(res2),
                   type="data", format=format))
    }else{
        return(new("track", dat=orderedGR(res), type="data", format=format))
    }
    }