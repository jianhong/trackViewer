addArrowMark <- function(pos=grid.locator(), angle=15, 
                         length=unit(.25, "inches"), col="red",
                         type="closed"){
    x <- pos$x
    y <- pos$y
    if(!(class(x) %in% c("numeric", "integer", "unit")) |
           length(x) < 1)
        stop("x is required as a numeric vector of x-coordinates of genome")
    if(!(class(y) %in% c("numeric", "integer", "unit")) |
           length(y) < 1)
        stop("y is required as a numeric vector of x-coordinates of score")
    if(length(x)!=length(y))
        stop("The length of x and y is not identical.")
    len <- length(x)
    trimLen <- function(obj, len){
        if(length(obj)<len) 
            obj <- rep(obj, len)[1:len]
        obj
    }
    col <- trimLen(col, len)
    angle <- trimLen(angle, len)
    length <- trimLen(length, len)
    type <- trimLen(type, len)
    for(i in 1:len){
        ar <- arrow(angle=angle[i], length=length[i], ends="first", type=type[i])
        xi <- as.numeric(convertUnit(x[i], "native", axisFrom="x"))
        yi <- as.numeric(convertUnit(y[i], "native", axisFrom="y"))
        l <- min(as.numeric(convertUnit(length[i], "native", axisFrom="x")), 
                 as.numeric(convertUnit(length[i], "native", axisFrom="y")))
        xi1 <- xi + l
        yi1 <- yi - l
        grid.lines(x=c(xi, xi1), y=c(yi, yi1), 
                   gp=gpar(col=col[i], fill=col[i]), arrow=ar,
                   default.units="native")
    }
    return(invisible(pos))
}