grid.pie <- function (x=.5, y=.5, 
                      radius=.8,
                      col=NULL,
                      border=NULL,
                      percent=NULL,
                      edges=100,
                      lwd=1) {
    if(length(percent)>0) percent <- unlist(percent[, sapply(percent, is.numeric)])
    if(length(percent)<1){
        percent <- 1
    }
    percent <- c(0, cumsum(percent)/sum(percent))
    dx <- diff(percent)
    nx <- length(dx)
    if (is.null(col)) 
        col <- c("white", "lightblue", "mistyrose", "lightcyan", 
                 "lavender", "cornsilk")
    if (!is.null(col)) 
        col <- rep_len(col, nx)
    if (!is.null(border)) 
        border <- rep_len(border, nx)
    twopi <- 2 * pi
    ratio.yx <- 1/as.numeric(convertX(unit(1, "snpc"), "npc"))
    t2xy <- function(t) {
        t2p <- twopi * t + pi/2
        list(x = radius * cos(t2p), y = radius * sin(t2p) * ratio.yx)
    }
    for (i in 1L:nx) {
        n <- max(2, floor(edges * dx[i]))
        P <- t2xy(seq.int(percent[i], percent[i + 1], length.out = n))
        grid.polygon(unit(c(P$x, 0)+x,"npc"), unit(c(P$y, 0)+y, "npc"), gp=gpar(col = border[i], fill = col[i], lwd=lwd))
    }
    invisible(NULL)
}

rgb2hex <- function(x){
    hex <- function(a) format(as.hexmode(a), width=2, upper.case=TRUE)
    paste0("#",hex(x[1]),hex(x[2]),hex(x[3]))
}
grid.lollipop <- function (x1=.5, y1=.5,
                           x2=.5, y2=.75,
                           y3=.04, y4=.02,
                           radius=.8,
                           col=NULL,
                           border=NULL,
                           percent=NULL,
                           edges=100,
                           type=c("circle", "pie", "pin", "pie.stack"),
                           ratio.yx=1,
                           pin=NULL,
                           scoreMax,
                           scoreType,
                           id=NA, id.col="black",
                           cex=1, lwd=1,
                           dashline.col="gray80"){
    stopifnot(is.numeric(c(x1, x2, y1, y2, y3, y4, radius, edges)))
    type <- match.arg(type)
    if(!type %in% c("pie", "pie.stack")){
        this.score <- if(length(percent$score)>0) max(percent$score, 1) else 1
        if(type=="circle"){
            grid.lines(x=c(x1, x1, x2, x2), y=c(y1, y2, y2+y3, y2+y3+y4+(this.score-.5)*2*radius*ratio.yx), 
                       gp=gpar(col=border))
            grid.lines(x=c(x2, x2), 
                       y=c(y2+y3+y4+(this.score-.5)*2*radius*ratio.yx, 
                           y2+y3+y4+scoreMax*ratio.yx), 
                       gp=gpar(col=dashline.col, lty=3))
        }else{
            grid.lines(x=c(x1, x1, x2, x2), y=c(y1, y2, y2+y3, y2+y3+y4+(this.score-.5)*2*radius*ratio.yx), 
                       gp=gpar(col=border))
        }
        
    }else{
        if(type=="pie.stack"){
            if(percent$stack.factor.first){
                grid.lines(x=c(x1, x1, x2, x2), y=c(y1, y2, y2+y3, y2+y3+y4), 
                           gp=gpar(col=border))
                grid.lines(x=c(x2, x2), 
                           y=c(y2+y3+y4, y2+y3+y4+scoreMax*ratio.yx),
                           gp=gpar(col=dashline.col, lty=3))
            }
        }else{
            grid.lines(x=c(x1, x1, x2, x2), y=c(y1, y2, y2+y3, y2+y3+y4), 
                       gp=gpar(col=border))
        }
    }
    
    if(length(pin)>0){
        if(length(border)>0) pin@paths[[2]]@rgb <- rgb2hex(col2rgb(border[1]))
        if(length(col)>0) pin@paths[[1]]@rgb <- rgb2hex(col2rgb(col[1]))
        if(length(col)>1) pin@paths[[3]]@rgb <- rgb2hex(col2rgb(col[2]))
    }
    switch(type,
           circle={
               if(scoreType){
                   for(i in 1:this.score){
                       grid.circle(x=x2, y=y2+y3+y4/2+2*radius*ratio.yx*(i-.5),
                                   r=radius*ratio.yx, 
                                   gp=gpar(col=border, fill=col, lwd=lwd))
                   }
               }else{
                   grid.circle(x=x2, y=y2+y3+(this.score-.5)*2*radius*ratio.yx+y4/2,
                               r=radius*ratio.yx, 
                               gp=gpar(col=border, fill=col, lwd=lwd))
                   if(!is.na(id)){
                       grid.text(label=id, x=x2, 
                                 y=y2+y3+(this.score-.5)*2*radius*ratio.yx+y4/2,
                                 just="centre", gp=gpar(col=id.col, cex=.75*cex))
                   }
               }
               },
           pie=grid.pie(x=x2, y=y2+y3+y4+radius*ratio.yx, 
                        radius = radius, 
                        col = col, 
                        border = border, 
                        percent=percent,
                        edges=edges,
                        lwd=lwd),
           pie.stack=grid.pie(x=x2, 
                              y=y2+y3+y4+(2*percent$stack.factor.order-1)*radius*ratio.yx, 
                              radius = radius, 
                              col = col, 
                              border = border, 
                              percent=percent[, !colnames(percent) %in% 
                                                  c("stack.factor.order", 
                                                    "stack.factor.first")],
                              edges=edges,
                              lwd=lwd),
           pin={
               grid.picture(picture=pin, x=x2, 
                            y=y2+y3+(this.score-.5)*2*radius*ratio.yx+y4/2,
                            width=2*radius*ratio.yx,
                            height=3*radius*ratio.yx+y4)
               if(!is.na(id)){
                   grid.text(label=id, x=x2, 
                             y=y2+y3+(this.score-.25)*2*radius*ratio.yx+2*y4/3,
                             just="centre", gp=gpar(col=id.col, cex=.5*cex))
               }
               },
           grid.pie(x=x2, y=y2+y3+y4+radius*ratio.yx, 
                    radius = radius, 
                    col = col, 
                    border = border, 
                    percent=percent,
                    edges=edges,
                    lwd=lwd))
}

jitterLables <- function(coor, xscale, lineW, weight=1.2){
    if(weight==1.2) stopifnot(order(coor)==1:length(coor))
    if(weight<0.5) return(coor)
    stopifnot(length(xscale)==2)
    pos <- (coor-xscale[1])/diff(xscale)
    pos.diff <- diff(c(0, pos, 1))
    idx <- which(pos.diff < weight*lineW)
    if(length(idx)<1){
        return(coor)
    }
    if(all(idx %in% c(1, length(pos)))){
        return(coor)
    }
    idx.diff <- diff(c(-1, idx))
    idx.grp <- rle(idx.diff)
    idx.grp$values[idx.grp$values==1] <- length(pos) + 1:sum(idx.grp$values==1)
    idx.grp <- inverse.rle(idx.grp)
    idx.grp.w <- which(idx.grp>length(pos))-1
    idx.grp.w <- idx.grp.w[idx.grp.w>0]
    idx.grp[idx.grp.w] <- idx.grp[idx.grp.w+1]
    idx.grp <- split(idx, idx.grp)
    flag <- as.numeric(names(idx.grp))>length(pos)
    idx.grp.mul <- lapply(idx.grp[flag], function(.ele){
        c(.ele[1]-1, .ele)
    })
    idx.grp.sin <- lapply(idx.grp[!flag], function(.ele){
        lapply(as.list(.ele), function(.ele){c(.ele-1, .ele)})
    })
    idx.grp.sin <- unlist(idx.grp.sin, recursive = FALSE)
    idx.grp <- c(idx.grp.mul, idx.grp.sin)
    
    adj.pos <- lapply(idx.grp, function(.ele){
        .ele <- .ele[.ele>0 & .ele<=length(pos)]
        this.pos <- pos[.ele]
        names(this.pos) <- .ele
        if(length(this.pos)%%2==1){
            center <- ceiling(length(this.pos)/2)
        }else{
            center <- length(this.pos)/2 + .5
        }
        if(length(this.pos)>5){ ## too much, how to jitter?
            this.pos <- this.pos + 
                ((1:length(this.pos))-center) * (weight-.1) * 
                lineW/ceiling(log(length(this.pos), 5))
        }else{
            this.pos <- this.pos + 
                ((1:length(this.pos))-center) * (weight-.1) * lineW
        }
        this.pos
    })
    names(adj.pos) <- NULL
    adj.pos <- unlist(adj.pos)
    coor[as.numeric(names(adj.pos))] <- adj.pos*diff(xscale)+xscale[1]
    
    Recall(coor, xscale=xscale, lineW=lineW, weight=weight-0.2)
}
