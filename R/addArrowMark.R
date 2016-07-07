addArrowMark <- function(pos=grid.locator(), label=NULL, angle=15, 
                         length=unit(.25, "inches"), col="red", cex=1, quadrant=4,
                         type="closed", vp=NULL){
    x <- pos$x
    y <- pos$y
    if(!(class(x) %in% c("numeric", "integer", "unit")) |
           length(x) < 1)
        stop("x is required as a numeric vector of x-coordinates of genome")
    if(!(class(y) %in% c("numeric", "integer", "unit")) |
           length(y) < 1)
        stop("y is required as a numeric vector of y-coordinates")
    if(length(x)!=length(y))
        stop("The length of x and y is not identical.")
    if(!quadrant %in% 1:4)
        stop("quadrant must be a integer in 1, 2, 3, 4")
    len <- length(x)
    lastTrackViewer <- getOption("LastTrackViewer")
    if(all(class(y) %in% c("numeric", "integer"))){
        if(all(floor(y)==y &ceiling(y)==y)){
            xy <- lastTrackViewer$xy
            yHeights <- c(0, lastTrackViewer$yHeights)
            yHeights <- cumsum(yHeights)
            yscales <- lastTrackViewer$yscales
            xscale <- lastTrackViewer$xscale
        }
    }
    
    trimLen <- function(obj, len){
        if(length(obj)<len) 
            obj <- rep(obj, len)[1:len]
        obj
    }
    col <- trimLen(col, len)
    angle <- trimLen(angle, len)
    length <- trimLen(length, len)
    type <- trimLen(type, len)
    cex <- trimLen(cex, len)
    for(i in 1:len){
        if(class(x[i]) %in% c("numeric", "integer")){
            if(findInterval(x[i], c(0, 1))==1){
                xi <- unit(x[i], "npc")
            }else{
                #if(!is.null(vp)) {
                    #r <- sort(vp$xscale)
                    #if(findInterval(x[i], r)==1){
                        xi <- unit(x[i], "native")
                    #}
                #}
            }
            if(class(xi)!="unit") stop("'pos$x' argument must be a unit object.")
        }else{
            xi <- x[i]
        }
        if(class(y[i]) %in% c("numeric", "integer")){
            if(findInterval(y[i], c(0, 1))==1){
                yi <- unit(y[i], "npc")
            }else{
                #if(!is.null(vp)) {
                    #r <- vp$yscale
                    #if(y[i]<=r[2] && y[i]>=r[1]){
                        yi <- unit(y[i], "native")
                    #}
                #}
                 if(length(lastTrackViewer)>0 && attributes(xi)$unit=="native"){
                     curViewStyle <- lastTrackViewer$viewerStyle
                     yvp <- list(vp1=viewport(x=0, y=yHeights[y[i]], 
                                              height=diff(yHeights)[y[i]],
                                              width=1, 
                                              just=c(0,0),
                                              clip="off"),
                                 vp2=viewport(x=0, y=lastTrackViewer$yHeightBottom[y[i]], 
                                              height=1 - lastTrackViewer$yHeightTop[y[i]] - 
                                                  lastTrackViewer$yHeightBottom[y[i]],
                                              width=1, 
                                              just=c(0,0),
                                              clip="off", 
                                              xscale=xscale, 
                                              yscale=yscales[[y[i]]]))
                     yi <- xy[[y[i]]]
                     yi <- yi$y[c(0, yi$x[-length(yi$x)])<=x[i] & c(yi$x[-1], .Machine$integer.max)>=x[i]]
                     if(length(yi)>0){
                         yi <- yi[1]
                     }else{
                         yi <- 0
                     }
                     yi <- unit(yi, "native")
                 }
            }
            if(class(yi)!="unit") stop("'pos$y' argument must be a unit object.")
        }else{
            yi <- y[i]
        }
        ar <- arrow(angle=angle[i], length=length[i], ends="first", type=type[i])
        hjust <- switch(quadrant, 0, 1, 1, 0, 0)
        if(!is.null(vp)){
            pushViewport(vp)
            popv <- FALSE
            if(length(lastTrackViewer)>0 && attributes(xi)$unit=="native" && 
               (class(y[i]) %in% c("numeric", "integer")) & findInterval(y[i], c(0, 1))!=1){
                pushViewport(yvp$vp1)
                pushViewport(yvp$vp2)
                popv <- TRUE
            }
            xi <- as.numeric(convertUnit(xi, "native", axisFrom="x"))
            yi <- as.numeric(convertUnit(yi, "native", axisFrom="y"))
            if(popv){
                xi1 <- switch(quadrant,
                              xi + as.numeric(convertUnit(length[i], "native", axisFrom="x")) - xscale[1],
                              xi - as.numeric(convertUnit(length[i], "native", axisFrom="x")) + xscale[1],
                              xi - as.numeric(convertUnit(length[i], "native", axisFrom="x")) + xscale[1],
                              xi + as.numeric(convertUnit(length[i], "native", axisFrom="x")) - xscale[1],
                              xi + as.numeric(convertUnit(length[i], "native", axisFrom="x")) - xscale[1])
                yi1 <- switch(quadrant,
                              yi + as.numeric(convertUnit(length[i], "native", axisFrom="y")) - yscales[[y[i]]][1],
                              yi + as.numeric(convertUnit(length[i], "native", axisFrom="y")) - yscales[[y[i]]][1],
                              yi - as.numeric(convertUnit(length[i], "native", axisFrom="y")) + yscales[[y[i]]][1],
                              yi - as.numeric(convertUnit(length[i], "native", axisFrom="y")) + yscales[[y[i]]][1],
                              yi - as.numeric(convertUnit(length[i], "native", axisFrom="y")) + yscales[[y[i]]][1])
            }else{
                xi1 <- switch(quadrant,
                              xi + as.numeric(convertUnit(length[i], "native", axisFrom="x")) - vp$xscale[1],
                              xi - as.numeric(convertUnit(length[i], "native", axisFrom="x")) + vp$xscale[1],
                              xi - as.numeric(convertUnit(length[i], "native", axisFrom="x")) + vp$xscale[1],
                              xi + as.numeric(convertUnit(length[i], "native", axisFrom="x")) - vp$xscale[1],
                              xi + as.numeric(convertUnit(length[i], "native", axisFrom="x")) - vp$xscale[1])
                yi1 <- switch(quadrant,
                              yi + as.numeric(convertUnit(length[i], "native", axisFrom="y")) - vp$yscale[1],
                              yi + as.numeric(convertUnit(length[i], "native", axisFrom="y")) - vp$yscale[1],
                              yi - as.numeric(convertUnit(length[i], "native", axisFrom="y")) + vp$yscale[1],
                              yi - as.numeric(convertUnit(length[i], "native", axisFrom="y")) + vp$yscale[1],
                              yi - as.numeric(convertUnit(length[i], "native", axisFrom="y")) + vp$yscale[1])
            }
            
            grid.lines(x=c(xi, xi1), y=c(yi, yi1), 
                       gp=gpar(col=col[i], fill=col[i]), arrow=ar,
                       default.units="native")
            if(!is.null(label[i])) grid.text(label=label[i], x = xi1, y = yi1, hjust = hjust,
                                          default.units = "native", gp = gpar(col=col[i], cex=cex[i]))
            if(popv){
                popViewport()
                popViewport()
            }
            popViewport()
        }else{
            xi <- as.numeric(convertUnit(xi, "native", axisFrom="x"))
            yi <- as.numeric(convertUnit(yi, "native", axisFrom="y"))
            l <- min(as.numeric(convertUnit(length[i], "native", axisFrom="x")),
                     as.numeric(convertUnit(length[i], "native", axisFrom="y")))
            xi1 <- switch(quadrant,
                          xi + l,
                          xi - l,
                          xi - l,
                          xi + l,
                          xi + l)
            yi1 <- switch(quadrant,
                          yi - l,
                          yi - l,
                          yi + l,
                          yi + l,
                          yi - l)
            
            grid.lines(x=c(xi, xi1), y=c(yi, yi1), 
                       gp=gpar(col=col[i], fill=col[i]), arrow=ar,
                       default.units="native")
            if(!is.null(label[i])) grid.text(label=label[i], x = xi1, y = yi1, hjust = hjust, 
                                          default.units = "native", 
                                          gp = gpar(col=col[i], cex=cex[i]))
        }
    }
    return(invisible(pos))
}