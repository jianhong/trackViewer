importScore <- function(file, file2, 
                        format=c("BED", "bedGraph", "WIG", "BigWig")){
    if(missing(file))
        stop("file is required.")
    ##    on.exit(closeAllConnections())
    format <- match.arg(format)
    #    res <- GRanges(score=numeric(0))
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
    readWIG <- function(buf){
        ##remove browser lines and track lines
        buf <- gsub("^\\s+", "", buf)
        buf <- gsub("\\s+$", "", buf)
        buf <- buf[grepl("^(variableStep|fixedStep|([0-9]+))", buf)]
        infoLine <- grep("Step", buf)
        if(length(infoLine)>0){
            if(infoLine[1]!=1){
                ##1 line is the info line
                stop("WIG file must contain track definition line, 
                     which should start by variableStep or fixedStep.")
            }
            }else{
                stop("WIG file must contain track definition line, 
                     which should start by variableStep or fixedStep.")
        }
        while(infoLine[length(infoLine)]==length(buf)){
            buf <- buf[-length(buf)]
        }
        block <- c(infoLine, length(buf)+1)
        dif <- diff(block)
        block <- rep(infoLine, dif)
        buf <- split(buf, block)
        r <- pblapply(buf, function(.ele) {
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
        buf <- lapply(buf, "[", -1)
        buf <- CharacterList(buf, compress=TRUE)
        GRanges(seqnames=r[,1], 
                ranges=IRanges(start=as.numeric(r[,2]),
                               end=as.numeric(r[,3])),
                score=buf,
                span=as.numeric(r[,4]),
                step=as.numeric(r[,5]),
                structure=r[,6])
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
        GRanges(seqnames=buf[,1], 
                ranges=IRanges(start=as.numeric(buf[,2])+1,##0 based, half open
                               end=as.numeric(buf[,3])),
                strand=buf[,6],
                score=as.numeric(buf[,5]))
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
        GRanges(seqnames=buf[,1], 
                ranges=IRanges(start=as.numeric(buf[,2])+1,##0 based, half open
                               end=as.numeric(buf[,3])),
                score=as.numeric(buf[,4]))
    }
    readbedGraph <- function(buf){
        readFourCols(buf)
    }
    readBigWig <- function(file){
        import(con=file, format="BigWig")
    }
    readFile <- function(file){
        s <- file.info(file)$size
        buf <- readChar(file, s, useBytes=TRUE)
        buf <- strsplit(buf, "\n", fixed=TRUE, useBytes=TRUE)[[1]]
        return(buf)
    }
    readFiles <- function(file, format){
        FUN <- get(paste("read", format, sep=""))
        if(format=="BigWig"){
            if(.Platform$OS.type=="windows")
                stop("rtracklayer can't read bigWig files on a Windows computer. Type in  ?`BigWigFile-class` to get the help.")
            res <- unique(FUN(file))
        }else{
            buf <- readFile(file)
            res <- unique(FUN(buf))
        }
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