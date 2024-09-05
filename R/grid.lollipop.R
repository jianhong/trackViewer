grid.pie <- function (x=.5, y=.5, 
                      radius=.8,
                      col=NULL,
                      border=NULL,
                      percent=NULL,
                      edges=100,
                      lwd=1,
                      alpha=1) {
    if(length(percent)>0) percent <- unlist(percent[, sapply(percent, is.numeric)])
    if(length(percent)<1){
        percent <- 1
    }
    percent <- c(0, cumsum(percent)/sum(percent))
    if(any(is.na(percent))){
      warning("There are events with NA number after calculating the percentage.",
              "Please make sure all the events must contain at least one values greater than 0")
      percent[is.na(percent)] <- 0
    }
    if(any(percent<0)){
      stop("There are events with number smaller than 0 which is not expected.",
           "Please double check all the input values.")
    }
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
    ratio.yx <- getYXratio()
    t2xy <- function(t) {
        t2p <- twopi * t + pi/2
        list(x = radius * cos(t2p), y = radius * sin(t2p) * ratio.yx)
    }
    for (i in 1L:nx) {
        n <- max(2, floor(edges * dx[i]))
        P <- t2xy(seq.int(percent[i], percent[i + 1], length.out = n))
        grid.polygon(unit(c(P$x, 0)+x,"npc"), unit(c(P$y, 0)+y, "npc"), gp=gpar(col = border[i], fill = col[i], lwd=lwd, alpha=alpha))
    }
    invisible(NULL)
}

rgb2hex <- function(x){
    hex <- function(a) format(as.hexmode(a), width=2, upper.case=TRUE)
    if(length(x==3))
      paste0("#",hex(x[1]),hex(x[2]),hex(x[3]))
    else
      paste0("#",hex(x[1]),hex(x[2]),hex(x[3]),hex(x[4]))
}
grid.lollipop <- function (x1=.5, y1=.5,
                           x2=.5, y2=.75,
                           y3=.04, y4=.02,
                           radius=.8,
                           col=NULL,
                           border=NULL,
                           percent=NULL,
                           edges=100,
                           type=c("circle", "pie", "pin", "pie.stack", "flag"),
                           ratio.yx=1,
                           pin=NULL,
                           scoreMax,
                           scoreType,
                           id=NA,
                           cex=1, lwd=1,
                           dashline.col="gray80",
                           side=c("top", "bottom"),
                           alpha=NULL,
                           shape=shape){
    side <- match.arg(side)
    stopifnot(is.numeric(c(x1, x2, y1, y2, y3, y4, radius, edges)))
    type <- match.arg(type)
    side <- side!="top"
    if(!type %in% c("pie", "pie.stack")){
        this.score <- if(length(percent$score)>0) percent$score else 1
        if(type=="circle"){
            y0 <- c(y1, y2, y2+y3, y2+y3+y4+(this.score-1)*2*radius*ratio.yx)
            if(scoreType) y0[4] <- y2+y3+y4
            if(side) y0 <- 1 - y0
            grid.lines(x=c(x1, x1, x2, x2), y=y0, 
                       gp=gpar(col=border, lwd=lwd))
            y0 <- c(y2+y3+y4+this.score*2*radius*ratio.yx, 
                    y2+y3+y4+scoreMax*ratio.yx)
            if(scoreType) y0[1] <- y2+y3+y4+this.score*2*radius*ratio.yx
            if(side) y0 <- 1 - y0
            grid.lines(x=c(x2, x2), 
                       y=y0, 
                       gp=gpar(col=dashline.col, lty=3, lwd=lwd))
        }else{
            y0 <- c(y1, y2, y2+y3, y2+y3+y4+(this.score-.5)*2*radius*ratio.yx)
            if(side) y0 <- 1 - y0
            grid.lines(x=c(x1, x1, x2, x2), y=y0, 
                       gp=gpar(col=border, lwd=lwd))
        }
        
    }else{
        if(type=="pie.stack"){
            if(percent$stack.factor.first){
                y0 <- c(y1, y2, y2+y3, y2+y3+y4)
                if(side) y0 <- 1 - y0
                grid.lines(x=c(x1, x1, x2, x2), y=y0, 
                           gp=gpar(col=border, lwd=lwd))
                y0 <- c(y2+y3+y4, y2+y3+y4+scoreMax*ratio.yx)
                if(side) y0 <- 1 - y0
                grid.lines(x=c(x2, x2), 
                           y=y0,
                           gp=gpar(col=dashline.col, lty=3, lwd=lwd))
            }
        }else{
            y0 <- c(y1, y2, y2+y3, y2+y3+y4)
            if(side) y0 <- 1 - y0
            grid.lines(x=c(x1, x1, x2, x2), y=y0, 
                       gp=gpar(col=border, lwd=lwd))
        }
    }
    if(length(pin)>0){
        if(length(border)>0) pin@paths[[2]]@rgb <- rgb2hex(col2rgb(border[1]))
        if(length(col)>0) pin@paths[[1]]@rgb <- rgb2hex(col2rgb(col[1]))
        if(length(col)>1) pin@paths[[3]]@rgb <- rgb2hex(col2rgb(col[2]))
    }
    switch(type,
           circle={
               if(length(border)==0) border <- "black"
               if(length(col)==0) col <- "white"
               if(scoreType){
                   for(i in seq_len(this.score)){
                       y0 <- y2+y3+y4+2*radius*ratio.yx*(i-.5)
                       if(side) y0 <- 1 - y0
                       if(length(shape)==this.score){
                         ## multiple shape for each point
                         this_shape <- shape[i]
                       }else{
                         this_shape <- ifelse(length(shape)>0, shape[1], 'circle')
                       }
                       if(length(border)==this.score){
                         ## multiple border for each point
                         this_border <- border[i]
                       }else{
                         this_border <- border[1]
                       }
                       if(length(col)==this.score){
                         ## multiple color for each point
                         this_col <- col[i]
                       }else{
                         this_col <- col[1]
                       }
                       if(length(lwd)==this.score){
                         ## multiple alpha for each point
                         this_lwd <- lwd[i]
                       }else{
                         this_lwd <- lwd[1]
                       }
                       if(length(alpha)==this.score){
                         ## multiple alpha for each point
                         this_alpha <- alpha[i]
                       }else{
                         this_alpha <- alpha[1]
                       }
                       if(length(this_shape)!=1) this_shape <- 'circle'
                       this_gp <- gpar(col=this_border, fill=this_col, lwd=this_lwd, alpha=this_alpha)
                       switch(this_shape, #"circle", "square", "diamond", "triangle_point_up", "star", or "triangle_point_down"
                              circle=grid.circle1(x=x2, y=y0,
                                                  r=radius*ratio.yx*cex, 
                                                  gp=this_gp),
                              square=grid.square(x=x2, y=y0,
                                                  r=radius*ratio.yx*cex, 
                                                  gp=this_gp),
                              diamond=grid.diamond(x=x2, y=y0,
                                                  r=radius*ratio.yx*cex, 
                                                  gp=this_gp),
                              triangle_point_up=grid.triangle_point_up(x=x2, y=y0,
                                                  r=radius*ratio.yx*cex, 
                                                  gp=this_gp),
                              triangle_point_down=grid.triangle_point_down(x=x2, y=y0,
                                                  r=radius*ratio.yx*cex, 
                                                  gp=this_gp),
                              star=grid.star(x=x2, y=y0,
                                             r=radius*ratio.yx*cex, 
                                             gp=this_gp),
                              grid.circle1(x=x2, y=y0,
                                           r=radius*ratio.yx*cex, 
                                           gp=this_gp))
                       
                   }
               }else{
                   y0 <- y2+y3+y4+(this.score-.5)*2*radius*ratio.yx
                   if(side) y0 <- 1 - y0
                   switch(shape,
                          circle=grid.circle1(x=x2, y=y0,
                                              r=radius*ratio.yx*cex, 
                                              gp=gpar(col=border, fill=col, lwd=lwd, alpha=alpha)),
                          square=grid.square(x=x2, y=y0,
                                             r=radius*ratio.yx*cex, 
                                             gp=gpar(col=border, fill=col, lwd=lwd, alpha=alpha)),
                          diamond=grid.diamond(x=x2, y=y0,
                                               r=radius*ratio.yx*cex, 
                                               gp=gpar(col=border, fill=col, lwd=lwd, alpha=alpha)),
                          triangle_point_up=grid.triangle_point_up(x=x2, y=y0,
                                                                   r=radius*ratio.yx*cex, 
                                                                   gp=gpar(col=border, fill=col, lwd=lwd, alpha=alpha)),
                          triangle_point_down=grid.triangle_point_down(x=x2, y=y0,
                                                                       r=radius*ratio.yx*cex, 
                                                                       gp=gpar(col=border, fill=col, lwd=lwd, alpha=alpha)),
                          star=grid.star(x=x2, y=y0,
                                         r=radius*ratio.yx*cex, 
                                         gp=gpar(col=border, fill=col, lwd=lwd, alpha=alpha)),
                          grid.circle1(x=x2, y=y0,
                                       r=radius*ratio.yx*cex, 
                                       gp=gpar(col=border, fill=col, lwd=lwd, alpha=alpha)))
                   if(!is.na(id$label)){
                       y0 <- y2+y3+(this.score-.5)*2*radius*ratio.yx+y4
                       if(side) y0 <- 1 - y0
                       id$gp$cex <- .75*id$gp$cex
                       id$x <- x2
                       id$y <- y0
                       do.call(grid.text, id)
                   }
               }
               },
           pie={
               y0 <- y2+y3+y4+radius*ratio.yx
               if(side) y0 <- 1 - y0
               grid.pie(x=x2, y=y0, 
                        radius = radius*cex, 
                        col = col, 
                        border = border, 
                        percent=percent,
                        edges=edges,
                        lwd=lwd, alpha=alpha)
               },
           pie.stack={
               y0 <- y2+y3+y4+(2*percent$stack.factor.order-1)*radius*ratio.yx
               if(side) y0 <- 1 - y0
               grid.pie(x=x2, 
                        y=y0, 
                        radius = radius*cex, 
                        col = col, 
                        border = border, 
                        percent=percent[, !colnames(percent) %in% 
                                            c("stack.factor.order", 
                                              "stack.factor.first")],
                        edges=edges,
                        lwd=lwd, alpha=alpha)
               },
           pin={
               y0 <- y2+y3+(this.score-.5)*2*radius*ratio.yx+y4/2
               if(side) y0 <- 1 - y0
               grid.picture(picture=pin, x=x2, 
                            y=y0,
                            width=2*radius*ratio.yx*cex,
                            height=3*radius*ratio.yx*cex+y4)
               if(!is.na(id$label)){
                   y0 <- y2+y3+(this.score-.25)*2*radius*ratio.yx+2*y4/3
                   id$gp$cex <- .5*id$gp$cex
                   id$x <- x2
                   id$y <- y0
                   do.call(grid.text, id)
               }
               },
           flag={
             if(is.na(id$label)){
               id$label <- " "
             }
             this.cex <- id$gp$cex
             LINEH <- as.numeric(convertY(unit(1, "line"), "npc"))*this.cex
             y0 <- y2+y3+(this.score-.5)*2*radius*ratio.yx+y4/2
             if(side) y0 <- 1 - y0
             LINEW <- as.numeric(convertX(stringWidth(paste0("o", id$label, "u")), "npc"))*this.cex
             LINEW <- LINEW * sign(cos(pi*id$rot/180))
             LINEH0 <- LINEW*ratio.yx*tan(pi*id$rot/180)
             grid.polygon(x=c(x2, x2+LINEW, x2+LINEW, x2),
                          y=c(y0, y0+LINEH0, y0+LINEH0+LINEH*1.25, y0+LINEH*1.25),
                          gp=gpar(fill=col, col=border, alpha=alpha))
             id$x <- x2+LINEW*.5
             id$y <- y0 + LINEH*.625+LINEH0*.5
             do.call(grid.text, id)
           },
           grid.pie(x=x2, y=y2+y3+y4+radius*ratio.yx, 
                    radius = radius*cex, 
                    col = col, 
                    border = border, 
                    percent=percent,
                    edges=edges,
                    lwd=lwd, alpha=alpha))
}

jitterLables <- function(coor, xscale, lineW, weight=1.2){
    if(weight==1.2) {
      stopifnot("Please sort your inputs by start position"= 
                  order(coor)==1:length(coor))
    }
    if(weight<0.5) return(coor)
    stopifnot(length(xscale)==2)
    pos <- convertX(unit(coor, "native"), "npc", valueOnly=TRUE)
    pos.diff <- diff(c(0, pos, 1))
    idx <- which(pos.diff < weight*lineW)
    if(length(idx)<1){
        return(coor)
    }
    if(all(idx %in% c(1, length(pos)+1))){
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

reAdjustLabels <- function(coor, lineW){
  # resort
  coor <- sort(coor)
  bins <- ceiling(1/lineW)
  pos <- convertX(unit(coor, "native"), "npc", valueOnly=TRUE)
  pos.bin <- cut(pos, c(-Inf, (0:bins)*lineW, Inf), labels=0:(bins+1), right=FALSE)
  
  ## split the coors by into clusters
  ## give the clusters with more idx more spaces if there are spaces between clusters
  tbl <- table(pos.bin)
  if(all(tbl<2)) return(coor)
  tbl.len <- length(tbl)
  if(tbl.len<3) return(coor)
  loops <- 1000
  loop <- 1
  while(any(tbl==0) && any(tbl>1) && loop < loops){
    tbl.bk <- tbl
    for(i in order(tbl.bk, decreasing=TRUE)){
      if(tbl[i]>1 && tbl.bk[i]==tbl[i]){
        if(i==1){
          if(tbl[2]<tbl[1]){
            half <- sum(tbl[1:2])/2
            tbl[2] <- ceiling(half)
            tbl[1] <- floor(half)
          }
        }else{
          if(i==tbl.len){
            if(tbl[tbl.len]>tbl[tbl.len-1]){
              half <- sum(tbl[(tbl.len-1):tbl.len])/2
              tbl[tbl.len-1] <- ceiling(half)
              tbl[tbl.len] <- floor(half)
            }
          }else{
            if(tbl[i-1]<tbl[i+1]){
              ## i-1 and i should be balanced
              half <- sum(tbl[(i-1):i])/2
              tbl[i-1] <- floor(half)
              tbl[i] <- ceiling(half)
            }else{
              half <- sum(tbl[i:(i+1)])/2
              tbl[i] <- floor(half)
              tbl[i+1] <- ceiling(half)
            }
          }
        }
      }
    }
    loop <- loop + 1
  }
  coef <- unlist(lapply(tbl, function(.ele){
    if(.ele==0) return(0)
    .ele <- seq(from=0, to=1, length.out=.ele+1)
    (.ele[-length(.ele)] + .ele[-1])/2
  }))
  coef <- coef[coef!=0]
  coor <- (rep(as.numeric(names(tbl)), tbl) - 1 + coef) * lineW
  coor <- convertX(unit(coor, "npc"), "native", valueOnly=TRUE)
  coor
}