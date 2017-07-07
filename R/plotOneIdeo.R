#' plot ideogram with data for one chromosome
#' @description plot ideogram with data for one chromosome 
#' @param ideo output of \link{loadIdeogram}.
#' @param dataList a \link[GenomicRanges]{GRangesList} of data to plot.
#' @param parameterList a list of parameters for each dataset in the dataList. 
#' The elements of the parameters could be xlabs, ylabs, etc. type could be
#' barplot, line, point, heatmap.
#' @param chrom A length 1 character vector of chromosome name.
#' @param colorSheme A character vector of giemsa stain colors.
#' @param gp parameters used for \link[grid]{grid.roundrect}.
#' @param ... parameters not used.
#' @import grid
#' @examples 
#' \dontrun{
#' ideo <- loadIdeogram("hg38")
#' library(rtracklayer)
#' library(grid)
#' dataList <- ideo[seqnames(ideo) %in% "chr1"]
#' dataList$score <- as.numeric(dataList$gieStain)
#' dataList <- dataList[dataList$gieStain!="gneg"]
#' dataList <- GRangesList(dataList, dataList)
#' grid.newpage()
#' plotOneIdeo(ideo, dataList, chrom="chr1")
#' }
#' 
#' 
plotOneIdeo <- 
  function(ideo, dataList, 
           parameterList=list(vp=plotViewport(margins=c(.1, 4.1, 1.1, .1)),
                              ideoHeight=unit(1/(1+length(dataList)), "npc"), 
                              vgap=unit(1, "lines"),
                              ylabs=seqlevels(ideo)[1], ylabsRot=90, 
                              ylabsPos=unit(2.5, "lines"),
                              xaxis=FALSE, yaxis=FALSE,
                              xlab="", 
                              types="barplot", heights=NULL, 
                              dataColumn="score", 
                              gps=gpar(col="black", fill="gray")), 
           chrom=seqlevels(ideo)[1], colorSheme=gieStain(), 
           gp=gpar(fill=NA, lwd=2), ...){
  stopifnot(class(dataList)=="GRangesList")
  stopifnot(class(parameterList)=="list")
  dataList <- lapply(dataList, function(.ele){
    .ele[seqnames(.ele) %in% chrom]
  })
  if(class(parameterList$xlab)=="list"){
    do.call(grid.text, args = parameterList$xlab)
  }else{
    if(parameterList$xlab!=""){
      grid.text(label=parameterList$xlab, y=unit(1.5, "lines"))
    }
  }
  ideo <- ideo[seqnames(ideo) %in% chrom]
  if(class(parameterList$vp)=="viewport"){
    if(identical(parameterList$vp$xscale, c(0, 1))){
      parameterList$vp$xscale <- c(1, max(end(ideo)))
    }
    pushViewport(parameterList$vp)
    
    if(length(parameterList$xaxis)){
      if(is.logical(parameterList$xaxis)){
        if(parameterList$xaxis[1]) grid.xaxis()
      }else{
        do.call(grid.xaxis, args = parameterList$xaxis)
      }
    }
  }
  ## split the plot region into ideo layer and data layer
  ideoHeight <- convertY(parameterList$ideoHeight, "npc", valueOnly = TRUE)
  ideoVP <- viewport(y=ideoHeight/2,
                     height=ideoHeight, 
                     name=paste0("ideoLayer", chrom))
  pushViewport(ideoVP)
  plotIdeo(ideo, chrom, colorSheme=colorSheme, gp=gp)
  upViewport()
  dataVP <- viewport(y=.5 + ideoHeight/2, height=1-ideoHeight, 
                     name=paste0("dataLayer", chrom))
  pushViewport(dataVP)
  if(length(parameterList$heights)){
    if(length(parameterList$heights)!=length(dataList)){
      vpH <- rep(parameterList$heights, length(dataList))[1:length(dataList)]
    }else{
      vpH <- parameterList$heights
    }
  }else{
    vpH <- rep(1/length(dataList), length(dataList))
  }
  if(length(parameterList$types)){
    if(length(parameterList$types)!=length(dataList)){
      types <- rep(parameterList$types, length(dataList))[1:length(dataList)]
    }else{
      types <- parameterList$types
    }
  }else{
    types <- rep("barplot", length(dataList))
  }
  if(length(parameterList$dataColumn)){
    if(length(parameterList$dataColumn)!=length(dataList)){ 
      dataColumn <- rep(parameterList$dataColumn, 
                        length(dataList))
      if(!is.list(parameterList$dataColumn)){
        dataColumn <- as.list(dataColumn)
      }
      dataColumn <- dataColumn[1:length(dataList)]
    }else{
      if(!is.list(parameterList$dataColumn)){
        dataColumn <- as.list(parameterList$dataColumn)
      }else{
        dataColumn <- parameterList$dataColumn
      }
    }
  }else{
    dataColumn <- as.list(rep("score", length(dataList)))
  }
  if(length(parameterList$gps)){
    if(length(parameterList$gps)!=length(dataList)){
      gps <- rep(list(parameterList$gps), length(dataList))[1:length(dataList)]
    }else{
      if(class(parameterList$gps)=="gpar"){
        gps <- rep(list(parameterList$gps), length(dataList))
      }else{
        gps <- parameterList$gps
      }
    }
  }else{
    gps <- rep(list(gpar()), length(dataList))
  }
  vgap <- convertY(parameterList$vgap, "npc", valueOnly = TRUE)
  for(i in seq_along(dataList)){
    data.gr <- dataList[[i]]
    mcols(data.gr) <- mcols(data.gr)[, dataColumn[[i]]]
    if(length(mcols(data.gr))>0){
        yrg <- range(as.numeric(as.matrix(mcols(data.gr))))
        if(!any(is.infinite(yrg))){
            this.vgap <- vgap
            if(vpH[i]-vgap<=0) this.vgap <- 0
            dataVPi <- viewport(y=sum(vpH[-(i:length(vpH))])+vpH[i]*.5+this.vgap, 
                                height = vpH[i]-this.vgap,
                                xscale=parameterList$vp$xscale, 
                                yscale=c(min(0, yrg[1]), max(yrg[2], 0)),
                                name=paste0("dataLayer", chrom, "sub", i))
            pushViewport(dataVPi)
            gridPlot(data.gr, gp=gps[[i]], types[i], parameterList$vp$xscale)
            upViewport()
        }
    }
  }
  upViewport()
  if(class(parameterList$vp)=="viewport"){
    upViewport()
  }
  if(length(parameterList$ylabs)){
    ylabs <- parameterList$ylabs
    if(ylabs[1]!=""){
      autofitFontSize <- 0.5 * 
          convertY(unit(vpH[1], "npc"), "inch", valueOnly=TRUE)/
            convertY(unit(1, "lines"), "inch", valueOnly = TRUE)
      grid.text(label=ylabs[1], x=parameterList$ylabsPos, 
                rot = parameterList$ylabsRot, 
                gp = gpar(cex=autofitFontSize),
                vp=viewport(y=ideoHeight/2,
                            height=ideoHeight))
      if(length(ylabs)>1){
        pushViewport(viewport(y=.5 + ideoHeight/2, height=1-ideoHeight))
        for(i in 2:length(ylabs)){
          this.vgap <- vgap
          if(vpH[i-1]-vgap<=0) this.vgap <- 0
          grid.text(label=ylabs[i], x=parameterList$ylabsPos, 
                    y=unit(.5, "npc"), rot = parameterList$ylabsRot, 
                    gp=gpar(cex=autofitFontSize),
                    vp=viewport(y=sum(vpH[-((i-1):length(vpH))])+vpH[i-1]*.5, 
                                height = vpH[i-1]-this.vgap))
        }
        popViewport()
      }
    }
  }
}

#' plot GRanges metadata 
#' @description plot GRanges metadata for different types
#' @param gr an object of \link[GenomicRanges]{GRanges} with metadata. All 
#' metadata must be numeric.
#' @param gp an object of \link[grid]{gpar}
#' @param type type of the figure, could be barplot, line, point and heatmap
#' @param xscale x scale of the viewport
#' @import grid
#' 

gridPlot <- function(gr, gp, type, xscale){
  if(length(gr)==0) return()
  x0 <- start(gr)
  x1 <- end(gr)
  x <- as.numeric(rbind(x0, x1))
  xc <- (x0+x1)/2
  y0 <- mcols(gr)
  if(ncol(y0)<1){
    return()
  }
  y <- as.matrix(y0[rep(1:nrow(y0), each=2),])
  yc <- as.numeric(y0[, 1])
  w <- width(gr)
  switch(type,
         barplot=grid.rect(x=xc, y=yc/2, 
                           width=w, height=yc,
                           gp=gp, default.units = "native"),
         line=grid.polyline(x=x, y=y, gp=gp, 
                            default.units = "native"),
         point=grid.points(x=xc, y=yc, gp=gp, 
                           default.units = "native"),
         heatmap={
           div <- 10000
           if(diff(xscale)<div) div <- diff(xscale)
           ww <- ceiling(diff(xscale)/div)
           img <- matrix(NA, nrow=ncol(y0), ncol=div)
           wt <- tile(GRanges(as.character(seqnames(gr)[1]),
                              IRanges(xscale[1], xscale[2])), 
                      n=div)[[1]]
           ol <- findOverlaps(wt, gr, ignore.strand=TRUE)
           qol <- unique(queryHits(ol))
           if(length(gp$breaks)){
             cols <- gp$col
             v <- gp$breaks
           }else{
             rg <- range(as.numeric(as.matrix(y0)))
             cols <- gp$col
             if(length(cols)==0){
               cols <- colorRampPalette(c("green", "black", "red"))(100)
             }else{
               if(length(cols)==1){
                 cols <- if(cols=="white") 
                   colorRampPalette(c("black", cols))(100) else
                     colorRampPalette(c("white", cols))(100)
               }else{
                 cols <- colorRampPalette(cols)(100)
               }
             }
             v <- seq(rg[1]-1, rg[2]+1, length.out = 101)
           }
           
           group <- merge(as.data.frame(ol), data.frame(queryHits=1:div), 
                          all=TRUE)
           y00 <- as.data.frame(y0)
           y01 <- matrix(nrow=nrow(group), ncol=ncol(y00))
           y01[!is.na(group$subjectHits), ] <- 
             y00[group$subjectHits[!is.na(group$subjectHits)], ]
           y01.rowsum <- rowsum(y01, group = group$queryHits)
           y01.groupcnt <- table(group$queryHits)
           stopifnot(identical(names(y01.groupcnt), rownames(y01.rowsum)))
           y01.rowsum <- y01.rowsum/as.numeric(y01.groupcnt)
           y00 <- apply(y01.rowsum, 2, function(.ele) 
             cols[findInterval(.ele, v)])
           stopifnot(all(dim(y00)==rev(dim(img))))
           img[, qol] <- y00[qol, ]
           grid.raster(img, width=1, height=1)
           })
}

