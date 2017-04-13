addGuideLine <- function(guideLine, col="gray", lty="dashed", lwd=1, vp=NULL){
    if(missing(guideLine) | 
           !(class(guideLine) %in% c("numeric", "integer")) |
           length(guideLine) < 1)
        stop("guideLine is required as a numeric vector of coordinates of genome")
    len <- length(guideLine)
    trimLen <- function(obj, len){
        if(length(obj)<len) 
            obj <- rep(obj, len)[1:len]
        obj
    }
    vpmultiple <- FALSE
    if(length(vp)>0){
        stopifnot(is(vp, "viewport"))
        if(is(vp, "vpTree")){
            ## check the range
            vpmultiple <- TRUE
        }
    }
    selectVP <- function(x, tree){
        xscales <- sapply(tree$children, function(.ele){
            seekViewport(names(.ele$children))
            current.viewport()$xscale
        }, simplify = FALSE)
        xscales <- do.call(cbind, xscales)
        i <- which(x>=xscales[1, ] & x<=xscales[2, ])
        if(length(i)<1) {
            message(x, " out of the range.")
            return(NULL)
        }
        seekViewport(paste0("panel.", i))
        current.viewport()
    }
    col <- trimLen(col, len)
    lty <- trimLen(lty, len)
    lwd <- trimLen(lwd, len)
    for(i in seq_along(guideLine)){
        if(vpmultiple){
            currentVP <- selectVP(guideLine[i], vp)
        }else{
            currentVP <- vp
        }
        grid.lines(x=guideLine[i], y=c(0, 1), 
                   gp=gpar(col=col[i], lty=lty[i], lwd=lwd[i]), 
                   default.units="native", vp=currentVP)
    }
    return(invisible())
}

getCurTrackViewport <- function(curViewerStyle, start, end){
    if(class(curViewerStyle)!="trackViewerStyle")
        stop("curViewerStyle must be an object of trackViewerStyle")
    margin <- curViewerStyle@margin
    xscale <- c(start, end)
    if(curViewerStyle@flip) xscale <- rev(xscale)
    return(viewport(x=margin[2], y=margin[1], 
                    height=1 - margin[1]- margin[3], 
                    width=1 -margin[2] - margin[4],
                    just=c(0,0), 
                    xscale=xscale, 
                    yscale=c(0,1)))
}