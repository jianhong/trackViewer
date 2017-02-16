#' plot ideogram with data
#' @description plot ideogram with data for multiple chromosomes 
#' @param ideo output of \link{loadIdeogram}.
#' @param dataList a \link[GenomicRanges]{GRangesList} of data to plot.
#' @param layout The layout of chromosomes. Could be a list with chromosome names
#' as its elements.
#' @param horiz a logical value. If FALSE, the ideograms are drawn vertically 
#' to the left. If TRUE, the ideograms are drawn horizontally at the bottom.
#' @param parameterList a list of parameters for each dataset in the dataList. 
#' The elements of the parameters could be xlabs, ylabs, etc. type could be
#' barplot, line, point, heatmap.
#' @param gp parameters used for \link[grid]{grid.roundrect}.
#' @param ... parameters not used.
#' @import grid
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
#'              parameterList = list(types="heatmap"))
#' }
#' 
#' 
ideogramPlot <- function(ideo, dataList, layout=NULL,
                         horiz=TRUE,
                         parameterList=
                           list(vp=plotViewport(margins=c(.1, 4.1, .3, .1)),
                                ideoHeight=unit(.5, "npc"), 
                                vgap=unit(.3, "lines"),
                                ylabs="auto", 
                                ylabsRot=ifelse(horiz, 0, 90), 
                                ylabsPos=unit(2.5, "lines"), 
                                xaxis=FALSE, yaxis=FALSE,
                                xlab="", 
                                types="barplot", heights=NULL, 
                                dataColumn="score", 
                                gps=gpar(col="black", fill="gray")), 
                         gp=gpar(fill=NA, lwd=2), ...){
  ## get the seqlengths
  seql <- seqlengths(ideo)
  seql <- seql[names(seql) %in% unique(seqnames(ideo))]
  stopifnot(is.list(parameterList))
  parameterList.default=list(vp=plotViewport(margins=c(.1, 4.1, .3, .1)),
                             ideoHeight=unit(.5, "npc"), 
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
  parameterList.cp <- parameterList <- parameterList.default
  if(any(is.na(seql))){
    stop("seqlengths contain NAs")
  }
  if(length(layout)==0){
    if(length(seql)==0){
      return()
    }
    if(length(seql)==1){
      plotOneIdeo(ideo, dataList, parameterList, chrom=names(seql), gp=gp, ...)
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
    for(j in 1:length(vpW)){
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
                  chrom=this.layout$seq[j], gp=gp, ...)
      upViewport()
    }
    upViewport()
  }
}