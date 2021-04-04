Y1pos <- function(SNPs.groups, xscale, lineW, base, cex, ypos, plotYaxis, heightMethod){
    ratio.Yy <- 1/diff(xscale)
    sg <- split(SNPs.groups, SNPs.groups$gps)
    sg <- sg[order(sapply(sg, length), decreasing=FALSE)]
    ratio.yx <- 1/as.numeric(convertX(unit(1, "snpc"), "npc"))
    t2xy <- function(t, rX, rY, angle) {
        t2p <- angle * t - (angle-pi)/2
        list(x = -rX * cos(t2p), y = rY * sin(t2p))
    }
    
    radius <- sapply(sg, function(.ele) floor(width(range(.ele))/2) * ratio.Yy *ratio.yx)
    radius <- radius[order(as.numeric(names(radius)))]
    
    stopifnot(length(plotYaxis)==1)
    if(plotYaxis){
      sg.scores <- sapply(sg, function(.ele) heightMethod(.ele$score))
    }else{
      sg.scores <- sapply(sg, function(.ele) heightMethod(.ele$score)) * lineW * ratio.yx
      sg.scores <- sg.scores[order(as.numeric(names(sg.scores)))]
      if(length(sg.scores)>1){
        for(i in 2:length(sg.scores)){
          if(sg.scores[i]>sg.scores[i-1]){
            if(sg.scores[i-1]*1.25+radius[i-1]+lineW>=sg.scores[i]){
              sg.scores[i] <- sg.scores[i-1]*1.25 + radius[i-1]+lineW
            }
          }else{
            if(sg.scores[i]*1.25+radius[i]+lineW>=sg.scores[i-1]){
              sg.scores[i] <- sg.scores[i-1]-radius[i]-sg.scores[i]/4-lineW
              if(sg.scores[i]<lineW) sg.scores[i] <- lineW
            }
          }
        }
      }
    }
    sg.scores <- sg.scores[names(sg)]
    radius <- radius[names(sg)]
    scoreRatio <- scoreMax0 <- max(sg.scores*5/4, na.rm = TRUE)
    yyscaleMax <- 1
    ## reset sg.scores to fit the paper
    if(scoreRatio + ypos + cex*lineW*ratio.yx + base > 1){
      scoreRatio <- (1-ypos-base-cex*lineW*ratio.yx)/scoreRatio
      sg.scores <- sg.scores * scoreRatio
      yyscaleMax <- scoreMax0
    }
    sg <- mapply(function(.ele, .h){
        .ele$Y1 <- .h + base
        if(length(.ele)>1){
            .range <- range(.ele)
            .rX <- floor(width(.range)/2)
            .rY <- .rX * ratio.Yy * ratio.yx
            .center <- start(.range) + .rX
            ## add 1/3 height to radius
            .rY <- .rY + .h / 4
            .rX <- .rX + .h * diff(xscale) / 4
            .percent <- (1:length(.ele)-1)/(length(.ele)-1)
            ## total angle
            perimeter <- 2 * pi * .rY *3/4 # max 3/2 pi
            angle <- pi * 3/2 * min(lineW*length(.ele)*cex/perimeter, 1)
            P <- t2xy(.percent, .rX, .rY, angle)
            .ele$X2 <- .center + P$x
            .ele$Y2 <- .ele$Y1 + P$y
            .ele$alpha <- (.5 - .percent) * angle
        }else{
            .ele$X2 <- start(.ele)
            .ele$Y2 <- .ele$Y1
            .ele$alpha <- 0
        }
        .ele
    }, sg, sg.scores)
    sg <- unlist(GRangesList(unname(sg)))
    sg <- sg[order(sg$idx)]
    sg$yyscaleMax <- yyscaleMax
    sg
}

grid.dandelion <- function(x0, y0, x1, y1, x2, y2, 
                           radius, col=NULL,
                           border=NULL,
                           percent=NULL,
                           edges=100,
                           alpha=0,
                           type=c("fan", "circle", "pie", "pin"),
                           ratio.yx=1,
                           pin=NULL,
                           scoreMax,
                           id=NA, id.col="black",
                           name=NULL, cex=1, lwd=1){
    stopifnot(is.numeric(c(x0, x1, x2, y0, y1, y2, radius, edges)))
    type <- match.arg(type)
    grid.lines(x=c(x0, x1, x2), y=c(y0, y1, y2), 
               gp=gpar(col=border, lwd=lwd))
    if(length(pin)>0){
        if(length(border)>0) pin@paths[[2]]@rgb <- rgb2hex(col2rgb(border[1]))
        if(length(col)>0) pin@paths[[1]]@rgb <- rgb2hex(col2rgb(col[1]))
        if(length(col)>1) pin@paths[[3]]@rgb <- rgb2hex(col2rgb(col[2]))
    }
    switch(type,
           fan={grid.dan(x=x2, y=y2, 
                         radius = radius, 
                         col = col, 
                         border = border, 
                         percent=percent,
                         edges=edges, 
                         alpha=alpha,
                         ratio.yx=ratio.yx, 
                         lwd=lwd)},
           circle={
               if(length(border)==0) border <- "black"
               if(length(col)==0) col <- "white"
               if(length(lwd)==0) lwd <- 1
               grid.circle1(x=x2, y=y2,
                           r=radius*ratio.yx, 
                           gp=gpar(col=border, fill=col, lwd=lwd))
               if(!is.na(id)){
                   grid.text(label=id, x=x2, 
                             y=y2+radius*ratio.yx,
                             just="centre", gp=gpar(col=id.col, cex=.75*cex))
               }
           },
           pie=grid.pie(x=x2, y=y2, 
                        radius = radius, 
                        col = col, 
                        border = border, 
                        percent=percent,
                        edges=edges, 
                        lwd=lwd),
           pin={
               grid.picture(picture=pin, x=x2, 
                            y=y2+1.5*radius*ratio.yx,
                            width=2*radius*ratio.yx,
                            height=3*radius*ratio.yx)
               if(!is.na(id)){
                   grid.text(label=id, x=x2, 
                             y=y2+1.33*radius*ratio.yx,
                             just="centre", gp=gpar(col=id.col, cex=.5*cex))
               }
           },
           grid.pie(x=x2, y=y2, 
                    radius = radius, 
                    col = col, 
                    border = border, 
                    percent=percent,
                    edges=edges,
                    lwd=lwd))
    if(length(name)>0){
        grid.lines(x=c(x2, x2), 
                   y=c(y2+radius*ratio.yx, 
                       y0+scoreMax+radius*ratio.yx), 
                   gp=gpar(col="gray80", lty=3, lwd=lwd))
        grid.text(x=x2, y=y0+scoreMax+radius*ratio.yx + radius*max(ratio.yx, 1.2), 
                  label = name, rot=90, just="left", gp=gpar(cex=cex))
    }
}

grid.dan <- function (x=.5, y=.5, 
                      radius=.8,
                      col=NULL,
                      border=NULL,
                      percent=NULL,
                      edges=100,
                      alpha=0, ratio.yx,
                      lwd=1) {
    if(length(percent$score)>0){
        if(percent$score>1){
            stop("if type=fan, the score must be a value between [0, 1].")
        }
        percent <- percent$score
    }else{
        percent <- 1
    }
    if(length(percent)<1){
        percent <- 1
    }
    if (is.null(col)) 
        col <- "white"
    twopi <- 2 * pi
    ratio.yx <- 1/as.numeric(convertX(unit(1, "snpc"), "npc"))
    t2xy <- function(t) {
        t2p <- twopi * t + pi/2 + alpha
        list(x = radius * cos(t2p), y = radius * sin(t2p) * ratio.yx)
    }
    n <- max(2, edges*percent*100)
    P <- t2xy(seq.int(-percent/2, percent/2, length.out = n))
    grid.polygon(unit(c(P$x, 0)+x,"npc"), unit(c(P$y, 0)+y, "npc"), gp=gpar(col = border[1], fill = col[1], lwd=lwd))
    invisible(NULL)
}