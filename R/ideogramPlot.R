#' plot ideogram with data
#' @description plot ideogram with data for multiple chromosomes 
#' @param ideo output of \link{loadIdeogram}.
#' @param dataList a \link[GenomicRanges:GRangesList-class]{GRangesList} of data to plot.
#' @param layout The layout of chromosomes. Could be a list with chromosome names
#' as its elements.
#' @param horiz a logical value. If FALSE, the ideograms are drawn vertically 
#' to the left. If TRUE, the ideograms are drawn horizontally at the bottom.
#' @param parameterList a list of parameters for each dataset in the dataList. 
#' The elements of the parameters could be xlabs, ylabs, etc. type could be
#' barplot, line, point, heatmap.
#' @param colorSheme A character vector of giemsa stain colors.
#' @param gp parameters used for \link[grid]{grid.roundrect}.
#' @param ... parameters not used.
#' @import grid
#' @importFrom GenomeInfoDb seqlevels seqlengths
#' @import S4Vectors
#' @export
#' @examples 
#' \dontrun{
#' ideo <- loadIdeogram("hg38")
#' library(rtracklayer)
#' library(grid)
#' dataList <- ideo
#' dataList$score <- as.numeric(dataList$gieStain)
#' dataList <- dataList[dataList$gieStain!="gneg"]
#' dataList <- GRangesList(dataList)
#' grid.newpage()
#' ideogramPlot(ideo, dataList, 
#'              layout=list("chr1", "chr2", c("chr3", "chr22"), 
#'                          c("chr4", "chr21"), c("chr5", "chr20"), 
#'                          c("chr6", "chr19"), c("chr7", "chr18"),
#'                          c("chr8", "chr17"), c("chr9", "chr16"),
#'                          c("chr10", "chr15"), c("chr11", "chr14"),
#'                          c("chr12", "chr13"), c("chrX", "chrY")),
#'              parameterList = list(types="heatmap", colorKeyTitle="sample1"))
#' }
#' 
#' 
ideogramPlot <- function(ideo, dataList, layout=NULL,
                         horiz=TRUE,
                         parameterList=
                           list(vp=plotViewport(margins=c(.1, 4.1, .3, .1)),
                                ideoHeight=unit(1/(1+length(dataList)), "npc"), 
                                vgap=unit(.3, "lines"),
                                ylabs="auto", 
                                ylabsRot=ifelse(horiz, 0, 90), 
                                ylabsPos=unit(2.5, "lines"), 
                                xaxis=FALSE, yaxis=FALSE,
                                xlab="", 
                                types="barplot", heights=NULL, 
                                dataColumn="score", 
                                gps=gpar(col="black", fill="gray")), 
                         colorSheme=gieStain(),
                         gp=gpar(fill=NA, lwd=2), ...){
  ## get the seqlengths
  seql <- seqlengths(ideo)
  seql <- seql[names(seql) %in% unique(seqnames(ideo))]
  stopifnot(is.list(parameterList))
  parameterList.default=list(vp=plotViewport(margins=c(.1, 4.1, .3, .1)),
                             ideoHeight=unit(1/(1+length(dataList)), "npc"), 
                             vgap=unit(.3, "lines"),
                             ylabs="auto", 
                             ylabsRot=ifelse(horiz, 0, 90), 
                             ylabsPos=unit(2.5, "lines"), 
                             xaxis=FALSE, yaxis=FALSE,
                             xlab="", 
                             types="barplot", heights=NULL, 
                             dataColumn="score", 
                             gps=gpar(col="black", fill="gray"))
  for(n in names(parameterList)){
    parameterList.default[[n]] <- parameterList[[n]]
  }
  if(parameterList.default$ylabs[1]!="auto"){
    if(!is.list(parameterList.default$ylabs)){
      stop("parameterList$ylabs must be a list with names of seqnames of ideo")
    }
    if(any(!seqlevels(ideo) %in% names(parameterList.default$ylabs))){
      stop("parameterList$ylabs must be a list with names of seqnames of ideo")
    }
  }
  if(all(parameterList.default$types=="heatmap")){
    gps <- parameterList.default$gps
    if(is(gps, "gpar")){
      cols <- list(gps$col)
      if(length(cols)==0){
        cols <- rep(list(c("green", "black", "red")), length(dataList))
      }
      if(length(cols)!=length(dataList)){
        cols <- rep(cols, length(dataList))
      }
    }else{
      if(length(gps)!=length(dataList)){
        stopifnot(all(sapply(gps, is.list)))
        gps <- rep(gps, length(dataList))[seq_along(dataList)]
      }
      cols <- lapply(gps, function(.ele) .ele$col)
    }
    
    ## get range of data points
    rg <- mapply(function(.ele, .cn) 
        range(as.numeric(as.matrix(mcols(.ele)[, .cn])), na.rm = TRUE),
        dataList, parameterList.default$dataColumn, SIMPLIFY = FALSE)
    ## break value of each color range
    v <- lapply(rg, function(.ele){
        if(diff(.ele)==0){
            c(.ele[1]-.ele[1]/10, .ele[2], .ele[2]+.ele[2]/10)
        }else{
            seq(.ele[1]-diff(.ele)/10, .ele[2]+diff(.ele)/10, length.out = 101)
        }
    })
    
    cols <- mapply(function(col, .ele){
      if(length(col)==0){
        col <- colorRampPalette(c("green", "black", "red"))(100)
      }else{
        if(length(col)==1){
          col <- if(col=="white") 
            colorRampPalette(c("black", col))(100) else
              colorRampPalette(c("white", col))(100)
        }else{
          col <- colorRampPalette(col)(100)
        }
      }
      if(diff(.ele)==0){
          col[c(1, 100)]
      }else{
          col
      }
    }, cols, rg, SIMPLIFY = FALSE)
    
    parameterList.default$gps <- mapply(function(col, breaks) 
      list(col=col, breaks=breaks), cols, v, SIMPLIFY = FALSE)
    labels <- parameterList.default$colorKeyTitle
    if(length(labels)==0){
      labels <- paste0("colorKey", 1:length(dataList))
    }
    labels.width <- convertX(stringWidth(labels), "inches", valueOnly = TRUE)
    colorKey.width <- 1
    oneChar.width <- convertX(unit(1, "lines"), "inches", valueOnly = TRUE)
    pushViewport(viewport(y=unit(1.75, "lines"), height=unit(3.5, "lines"),
                          width=unit(sum(labels.width) + 
                                       length(labels) +
                                       length(labels)*2*oneChar.width, 
                                     "inches")))
    labels.x.width <- colorKey.width + 2*oneChar.width + labels.width
    this.labels.width <- 0
    for(i in 1:length(labels)){
      pushViewport(viewport(x=unit(labels.x.width[i]/2 + this.labels.width, 
                                   'inches'),
                            width=unit(labels.x.width[i], "inches")))
      this.labels.width <- this.labels.width + labels.x.width[i]
      this.xscale <- range(v[[i]])
      vp <- viewport(x=unit(.5, "inches"), width=unit(1, "inches"),
                     y=unit(3, "lines"), height=unit(.5, "lines"),
                     xscale=this.xscale)
      pushViewport(vp)
      grid.raster(as.raster(matrix(cols[[i]], nrow=1)), 
                  width = 1, height = 1, interpolate=FALSE)
      grid.xaxis()
      upViewport()
      grid.text(label=labels[i], 
                vp = viewport(x=unit(colorKey.width+oneChar.width+
                                labels.width[i]/2, "inches"),
                              width=unit(labels.width[i], "inches"),
                              y=unit(3, "lines"), height=unit(1, "lines")))
      upViewport()
    }
    upViewport()
    
    pushViewport(plotViewport(margins = c(4.1, 0, 0, 0)))
    on.exit(upViewport())
  }
  parameterList.cp <- parameterList <- parameterList.default
  if(any(is.na(seql))){
    stop("seqlengths contain NAs")
  }
  if(length(layout)==0){
    if(length(seql)==0){
      return()
    }
    if(length(seql)==1){
      plotOneIdeo(ideo, dataList, parameterList, chrom=names(seql), 
                  colorSheme=colorSheme, gp=gp, ...)
      return()
    }
    ## 1 column
    layout <- as.list(names(seql))
  }
  stopifnot(is.list(layout))
  if(!horiz){
    canvas.w <- convertX(unit(1, "npc"), "inch", valueOnly = TRUE)
    canvas.h <- convertY(unit(1, "npc"), "inch", valueOnly = TRUE)
    tmp <- canvas.h
    canvas.h <- convertY(unit(canvas.w, "inch"), "npc")
    canvas.w <- convertX(unit(tmp, "inch"), "npc")
    pushViewport(viewport(width=canvas.w, 
                          height=canvas.h, 
                          angle=-90))
    on.exit(upViewport())
  }
  layout.df <- data.frame(seq=unlist(layout), 
                          group=rep(1:length(layout), sapply(layout, length)),
                          stringsAsFactors = FALSE)
  layout.df$seqlength <- seql[layout.df$seq]
  xscale <- max(sapply(split(layout.df$seqlength, layout.df$group), sum)) * 1.05
  vpH <- 1/length(layout)
  for(i in sort(unique(layout.df$group))){
    this.layout <- layout.df[layout.df$group==i, , drop=FALSE]
    vp <- viewport(y=1-vpH*(i-.5), height = vpH, name=paste0("layout", i))
    pushViewport(vp)
    vpW <- this.layout$seqlength/xscale
    if(length(vpW)==1){
      vpWgap <- 0
    }else{
      vpWgap <- (1/1.05-sum(vpW))/(length(vpW)-1)
    }
    for(j in seq_along(vpW)){
      vp <- viewport(x=sum(vpW[-(j:length(vpW))])+.5*vpW[j]+vpWgap*(j-1), 
                     width=vpW[j], name=paste0("layout", i, "_", j))
      pushViewport(vp)
      parameterList <- parameterList.cp
      if(parameterList$ylabs[1]=="auto"){
        parameterList$ylabs <- this.layout$seq[j]
      }else{
        parameterList$ylabs <- parameterList$ylabs[[this.layout$seq[j]]]
      }
      if(j==length(vpW) && length(vpW)>1){
        parameterList$vp$x <- 
          unit(convertX(unit(1, "npc") - parameterList$vp$x - parameterList$vp$width,
                   "lines", valueOnly = TRUE), "lines")
        parameterList$ylabsPos <- unit(1, "npc") - parameterList$ylabsPos
      }
      plotOneIdeo(ideo, dataList, parameterList, 
                  chrom=this.layout$seq[j], colorSheme=colorSheme, 
                  gp=gp, ...)
      upViewport()
    }
    upViewport()
  }
}