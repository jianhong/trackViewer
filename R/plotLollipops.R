convertHeight2NPCnum <- function(.ele){
  if(is(.ele, "unit")){
    return(convertHeight(.ele, unitTo="npc", valueOnly=TRUE))
  }else{
    if(is.list(.ele)){
      .ele <- sapply(.ele, function(.e){
        if(is(.e, "unit")){
          .e <- convertHeight(.e, unitTo="npc", valueOnly=TRUE)
        }
        .e[1]
      })
      return(unlist(.ele))
    }else{
      if(is.numeric(.ele)){
        return(.ele)
      }else{
        if(is.integer(.ele)){
          return(.ele)
        }else{
          return(.ele)
        }
      }
    }
  }
}
plotFeatures <- function(feature.splited, LINEH, bottomHeight, 
                         label_on_feature=FALSE){
    feature.height <- 0
    for(n in seq_along(feature.splited)){
        this.feature.height <- 
            max(c(feature.splited[[n]]$height/2, 
                  .0001)) + 0.2 * LINEH
        feature.height <- feature.height + this.feature.height
        ##baseline
        grid.lines(x=c(0, 1), y=c(bottomHeight+feature.height, 
                                  bottomHeight+feature.height))
        for(m in seq_along(feature.splited[[n]])){
            this.dat <- feature.splited[[n]][m]
            color <- if(is.list(this.dat$color)) this.dat$color[[1]] else 
                this.dat$color
            if(length(color)==0) color <- "black"
            fill <- if(is.list(this.dat$fill)) this.dat$fill[[1]] else 
                this.dat$fill
            if(length(fill)==0) fill <- "white"
            this.cex <- if(length(this.dat$cex)>0) this.dat$cex[[1]][1] else 1
            lwd <- if(length(this.dat$lwd)>0) this.dat$lwd[[1]][1] else 1
            this.feature.height.m <- 
                if(length(this.dat$height)>0) 
                    this.dat$height[[1]][1] else 
                        2*this.feature.height
            grid.rect(x=start(this.dat)-.1, y=bottomHeight+feature.height, 
                      width=width(this.dat)-.8, 
                      height=this.feature.height.m,
                      just="left", gp=gpar(col=color, fill=fill, lwd=lwd), 
                      default.units = "native")
            if(!is.null(names(this.dat))){
              if(label_on_feature & !is.na(names(this.dat)[1])){
                grid.text(x=(start(this.dat)+end(this.dat))/2, 
                          y=bottomHeight+feature.height,
                          just = "centre",
                          label = names(this.dat)[1],
                          gp= gpar(list(cex=this.cex * 
                                          this.feature.height.m/
                                          this.feature.height,
                                        color=color)), 
                          default.units = "native")
              }
            }
        }
        feature.height <- feature.height + this.feature.height
    }
    feature.height
}
getDrawLabelParam <- function(SNPs){
  label.parameter.draw <- rep(TRUE, length(SNPs))
  if(length(SNPs$label.parameter.draw)==length(SNPs)){
    label.parameter.draw <- vapply(SNPs$label.parameter.draw, `[`, i=1,
                                   FUN.VALUE = logical(1L))
  }
  label.parameter.draw
}
handleLabelParams <- function(SNPs, prefix="label.parameter.", cex=1,
                              ...){
  labels <- list(...)
  ## change the parameter by use definitions.
  for(label.parameter in names(labels)){
    label.para <- paste0(prefix, label.parameter)
    if(label.para %in% colnames(mcols(SNPs))){
      labels[[label.parameter]] <- mcols(SNPs)[, label.para]
    }
  }
  for(label.parameter in c("col", "fill", "alpha", "lty", "lwd", "lex",
                          "lineend", "linejoin", "linemitre",
                          "fontsize", "cex", "fontfamily", "fontface",
                          "lineheight", "font")){
    label.para <- paste0(prefix, label.parameter)
    if(label.para %in% colnames(mcols(SNPs))){
      labels$gp[[label.parameter]] <- mcols(SNPs)[, label.para]
    }
  }
  if(!"cex" %in% names(labels$gp)){
    labels$gp$cex <- cex
  }
  mergeList <- function(.ele){
    .ele <- do.call(list, .ele)
    .n <- unique(unlist(lapply(.ele, names)))
    .out <- list()
    if(length(.n)>0){
      for(.name in .n){
        .out[[.name]] <- sapply(.ele, function(.e){
          if(.name %in% names(.e)){
            .e[[.name]][1]
          }else{
            NA
          }
        })
      }
    }else{
      .n <- unique(names(.ele))
      for(.name in .n){
        .out[[.name]] <- unlist(.ele[names(.ele) %in% .name])
      }
    }
    .out
  }
  labels$gp <- mergeList(labels$gp)
  labels$gp[duplicated(names(labels$gp))] <- NULL
  labels$gp <- do.call(gpar, labels$gp)
  return(labels)
}
filterThisLabel <- function(this.label){
  na_label <- is.na(this.label$label)
  for(key in c("x", "just", "hjust", "vjust", "rot",
               "check.overlap", "default.units")){
    if(length(this.label[[key]])>1 &&
       length(this.label[[key]]==length(this.label$label))){
      this.label[[key]] <- this.label[[key]][!na_label]
    }
  }
  for(key in names(this.label$gp)){
    if(length(this.label$gp[[key]])>1 &&
       length(this.label$gp[[key]])==length(this.label$label)){
      this.label$gp[[key]] <- this.label$gp[[key]][!na_label]
    }
  }
  this.label$gp <- do.call(gpar, this.label$gp)
  this.label$label <- this.label$label[!na_label]
  this.label
}
plotLollipops <- function(SNPs, feature.height, bottomHeight, baseline, 
                          type, ranges, yaxis, yaxis.gp, scoreMax, scoreMax0, scoreType,
                          LINEW, cex, ratio.yx, GAP, pin, dashline.col,
                          side=c("top", "bottom"), jitter=c("node", "label"),
                          main=TRUE){
    side <- match.arg(side)
    jitter <- match.arg(jitter)
    if(side=="top"){
        pushViewport(viewport(y=bottomHeight,
                              height=1,
                              just="bottom",
                              xscale=c(start(ranges), 
                                       end(ranges)),
                              clip="off"))
    }else{
        pushViewport(viewport(y=bottomHeight+feature.height,
                              height=1,
                              just="top",
                              xscale=c(start(ranges), 
                                       end(ranges)),
                              yscale=c(1, 0),
                              clip="off"))
    }
    
    if(type=="pie.stack" && length(SNPs$stack.factor)>0){
        stopifnot(is.vector(SNPs$stack.factor, mode="character"))
        if(length(SNPs$stack.factor.order)>0 || 
           length(SNPs$stack.factor.first)>0){
            warning("stack.factor.order and stack.factor.first are used by this function!",
                    "The values in these column will be removed.")
        }
        ## condense the SNPs
        stack.factors <- unique(as.character(SNPs$stack.factor))
        stack.factors <- sort(stack.factors)
        stack.factors.order <- seq_along(stack.factors)
        names(stack.factors.order) <- stack.factors
        SNPs <- SNPs[order(as.character(seqnames(SNPs)), start(SNPs), 
                           as.character(SNPs$stack.factor))]
        SNPs$stack.factor.order <- stack.factors.order[SNPs$stack.factor]
        SNPs$stack.factor.first <- !duplicated(SNPs)
        SNPs.condense <- SNPs
        SNPs.condense$oid <- seq_along(SNPs)
        SNPs.condense$factor <- paste(as.character(seqnames(SNPs)), start(SNPs), end(SNPs))
        SNPs.condense <- split(SNPs.condense, SNPs.condense$factor)
        SNPs.condense <- lapply(SNPs.condense, function(.ele){
            .oid <- .ele$oid
            .gr <- .ele[1]
            mcols(.gr) <- NULL
            .gr$oid <- NumericList(.oid)
            .gr
        })
        SNPs.condense <- unlist(GRangesList(SNPs.condense), use.names = FALSE)
        SNPs.condense <- sort(SNPs.condense)
        lab.pos.condense <- start(SNPs.condense)
        label.parameter.draw <- rep(TRUE, length(SNPs.condense))
        if(length(SNPs$label.parameter.draw)==length(SNPs)){
          label.parameter.draw <- 
            SNPs$label.parameter.draw[vapply(SNPs.condense$oid, `[`, i=1,
                                             FUN.VALUE = numeric(1L))]
          label.parameter.draw <- vapply(label.parameter.draw, `[`, i=1,
                                         FUN.VALUE = logical(1L))
        }
        lab.pos.condense[label.parameter.draw] <-
          jitterLables(start(SNPs.condense)[label.parameter.draw], 
                                         xscale=c(start(ranges), end(ranges)), 
                                         lineW=LINEW*cex)
        lab.pos.condense[label.parameter.draw] <-
          reAdjustLabels(lab.pos.condense[label.parameter.draw],
                         lineW=LINEW*cex)
        condense.ids <- SNPs.condense$oid
        lab.pos <- rep(lab.pos.condense, elementNROWS(condense.ids))
        lab.pos <- lab.pos[order(unlist(condense.ids))]
    }else{
        lab.pos <- start(SNPs)
        label.parameter.draw <- getDrawLabelParam(SNPs)
        lab.pos[label.parameter.draw] <-
          jitterLables(start(SNPs)[label.parameter.draw], 
                       xscale=c(start(ranges), end(ranges)), 
                       lineW=LINEW*cex)
        lab.pos[label.parameter.draw] <-
          reAdjustLabels(lab.pos[label.parameter.draw], lineW=LINEW*cex)
    }
    if(length(SNPs)>0){
        yaxisat <- NULL
        yaxisLabel <- TRUE
        diameter <- LINEW*ratio.yx
        y2 <- feature.height
        y3 <- 4*GAP*cex
        y4 <- 2.5*GAP*cex
        if(length(yaxis)>1 && is.numeric(yaxis)){
            yaxisat <- yaxis
            if(length(names(yaxis))==length(yaxis)) yaxisLabel <- names(yaxis)
            yaxis <- TRUE
        }
        if(yaxis && scoreMax>1 && !type %in% c("pie", "pie.stack")){
            if(side=="top"){
              vp <- viewport(x=.5+ifelse(main, -1, 1) *LINEW,
                             y=feature.height + y3 + y4 + scoreMax*diameter/2,
                             width=1,
                             height=(scoreMax+1)*diameter,
                             yscale=c(0, scoreMax0+scoreMax0/scoreMax))
            }else{
              vp <- viewport(x=.5+ifelse(main, -1, 1) *LINEW,
                             y=1-(feature.height + y3 + y4 +
                                    scoreMax*diameter/2),
                             width=1,
                             height=(scoreMax+1)*diameter,
                             yscale=c(scoreMax0+scoreMax0/scoreMax, 0))
            }
          grid.yaxis(at=yaxisat,
                     label=yaxisLabel,
                     main = main,
                     gp=yaxis.gp,
                     vp=vp)
        }
        
        if(length(SNPs$alpha)==length(SNPs)){
          SNPs$alpha[is.na(SNPs$alpha)] <- 0
          if(all(is.numeric(SNPs$alpha))){
            if(any(SNPs$alpha>1)){## convert to 0-1
              SNPs$alpha <- SNPs$alpha/max(SNPs$alpha)
            }
          }else{ ## not correct format.
            SNPs$alpha <- as.numeric(factor(as.character(SNPs$alpha)))
            SNPs$alpha <- (SNPs$alpha+max(SNPs$alpha))/max(SNPs$alpha)/2
          }
        }else{
          SNPs$alpha <- NULL
        }
        if(type=="circle"){
          if(length(SNPs$shape)==length(SNPs)){
            ## shape could only be "circle", "square", "diamond", "triangle_point_up", "triangle_point_down"
            if(!all(unlist(SNPs$shape) %in% c("circle", "square", "diamond", "triangle_point_up", "triangle_point_down"))){
              message('shape must be "circle", "square", "diamond", "triangle_point_up", or "triangle_point_down"')
              if(is.list(SNPs$shape)){
                stop("The shape is a list, please confirm the format.",
                     'It must be a list with elements of "circle", "square", "diamond", "triangle_point_up", or "triangle_point_down"')
              }
              SNPs$shape <- as.numeric(factor(SNPs$shape))
              SNPs$shape <- rep(c("circle", "square", "diamond", "triangle_point_up", "triangle_point_down"), 
                                max(SNPs$shape))[SNPs$shape]
            }else{
              if(is.list(SNPs$shape)){
                if(any(lengths(SNPs$shape)[SNPs$score!=0]==0)){
                  stop("The shape is a list, but zero length of shape is detected.")
                }
              }
            }
            if(scoreType){
              if(is.list(SNPs$shape)){
                if(!all(lengths(SNPs$shape)[SNPs$score!=0]==
                        SNPs$score[SNPs$score!=0])){
                  warning('Not all the length of shape equal to score.')
                }
              }
            }
          }else{
            SNPs$shape <- NULL
          }
        }
        for(m in seq_along(SNPs)){
            this.dat <- SNPs[m]
            color <- if(is.list(this.dat$color)) this.dat$color[[1]] else this.dat$color
            border <- 
                if(is.list(this.dat$border)) this.dat$border[[1]] else this.dat$border
            fill <- if(is.list(this.dat$fill)) this.dat$fill[[1]] else this.dat$fill
            alpha <- if(length(this.dat$alpha)>0) this.dat$alpha[[1]] else 1
            lwd <- if(is.list(this.dat$lwd)) this.dat$lwd[[1]] else this.dat$lwd
            shape <- if(length(this.dat$shape)>0) this.dat$shape[[1]] else "circle"
            rot <- if(length(this.dat$label.rot)>0) this.dat$label.rot else 15
            this.cex <- if(length(this.dat$cex)>0) this.dat$cex[[1]][1] else 1
            this.dashline.col <- 
              if(length(this.dat$dashline.col)>0) this.dat$dashline.col[[1]][1] else dashline.col
            ## control plot the dash line or not
            if(length(names(this.dat))<1) this.dashline.col <- NA
            if(length(SNPs$label.parameter.label)==length(SNPs) && length(SNPs)>0){
              if(this.dat$label.parameter.label[[1]]=="" ||
                 is.na(this.dat$label.parameter.label[[1]])) this.dashline.col <- NA
            }
            if(length(SNPs$label.parameter.pfm)==length(SNPs) && length(SNPs)>0){
              if(is.null(SNPs$label.parameter.pfm[[m]])){
                this.dashline.col <- NA
              }
            }
            if(length(SNPs$label.parameter.draw)==length(SNPs) && length(SNPs)>0){
              if(!(SNPs$label.parameter.draw[[m]])){
                this.dashline.col <- NA
              }
            }
            id <- 
              handleLabelParams(this.dat, cex = this.cex, prefix = "node.label.",
                                label = if(is.character(this.dat$label)) this.dat$label else
                                  if(is.character(this.dat$node.label)) this.dat$node.label else NA,
                                rot = if(length(this.dat$label.rot)>0) this.dat$label.rot else ifelse(type=="flag", 15, 0),
                                gp = gpar(cex=this.cex,
                                          col = if(length(this.dat$label.col)>0) this.dat$label.col else "black"),
                                just = "centre",
                                hjust = .5,
                                vjust = .5)
            this.dat.mcols <- mcols(this.dat)
            this.dat.mcols <- cleanDataMcols(this.dat.mcols, type)

            grid.lollipop(x1=convertX(unit(start(this.dat), "native"), "npc", 
                                      valueOnly=TRUE),  
                          y1=baseline,
                          x2=convertX(unit(ifelse(jitter=="node", 
                                                  lab.pos[m], 
                                                  start(this.dat)), 
                                           "native"), "npc", valueOnly=TRUE), 
                          y2=y2,
                          y3=y3, y4=y4, 
                          radius=LINEW/2,
                          col=color,
                          border=border,
                          percent=this.dat.mcols,
                          edges=100,
                          type=type,
                          ratio.yx=ratio.yx,
                          pin=pin,
                          scoreMax=scoreMax * LINEW,
                          scoreType=scoreType,
                          id=id,
                          cex=this.cex, lwd=lwd, dashline.col=this.dashline.col,
                          side=side, alpha=alpha, shape=shape)

        }
        this.height <- getHeight(SNPs, 
                                 ratio.yx, LINEW, GAP, cex, type,
                                 scoreMax=scoreMax,
                                 level="data")
        labels.keep <- getDrawLabelParam(SNPs)
        SNPs <- SNPs[labels.keep]
        lab.pos <- lab.pos[labels.keep]
        if(length(names(SNPs))>0){
            if(type=="pie.stack"){
                ## unique lab.pos and SNPs
                idx <- !duplicated(names(SNPs))
                lab.pos <- lab.pos[idx]
                SNPs <- SNPs[idx]
            }
            this.label <- 
              handleLabelParams(SNPs, cex = cex,
                                prefix = "label.parameter.",
                                x = lab.pos,
                                label = names(SNPs),
                                just = ifelse(side=="top", "left", "right"),
                                hjust = NULL,
                                vjust = NULL,
                                rot = 90,
                                check.overlap = FALSE,
                                default.units = "native",
                                gp = gpar(cex=cex),
                                pfm = NULL,
                                font = "Helvetica-Bold",
                                fontface = "bold",
                                ic.scale = TRUE)
            if(jitter=="label"){
              ## add guide lines
              rased.height <- 4*GAP*cex
              guide.height <- 2.5*GAP*cex
              for(i in seq_along(SNPs)){
                this.dashline.col <- 
                  if(length(SNPs[i]$dashline.col)>0) 
                    SNPs[i]$dashline.col[[1]][1] else 
                      dashline.col
                if(length(this.label$label)<1 && length(this.label$pfm)<1){
                  next
                }
                if(length(this.label$label)==length(SNPs) && length(SNPs)>0){
                  if(is.na(this.label$label[i])){
                    next
                  }
                  if(this.label$label[i]==""){
                    next
                  }
                }
                if(length(this.label$pfm)==length(SNPs) && length(SNPs)>0){
                  if(is.null(this.label$pfm[[i]])){
                    next
                  }
                }
                grid.lines(x=c(start(SNPs[i]), this.label$x[i]), 
                           y=c(this.height+feature.height-cex*LINEW, 
                               this.height+feature.height+rased.height),
                           default.units = this.label$default.units,
                           gp=gpar(col=this.dashline.col, lty=3))
                grid.lines(x=c(this.label$x[i], this.label$x[i]),
                           y=c(this.height+rased.height+feature.height,
                               this.height+rased.height+
                                 guide.height+feature.height),
                           default.units = this.label$default.units,
                           gp=gpar(col=this.dashline.col, lty=3))
              }
              ## add this height
              this.height <- this.height + rased.height + guide.height
            }
            if(length(this.label$pfm)>0){
              if(!requireNamespace("motifStack", quietly = TRUE)){
                stop("When plot motifs as labels,",
                     " the Bioconductor package 'motifStack' is required!")
              }
              for(idx in seq_along(this.label$pfm)){
                if(!is.null(this.label$pfm[[idx]])){
                  this_cex <- ifelse(length(cex)==length(this.label$pfm),
                                     cex[[idx]], cex[1])
                  this_just <- ifelse(length(this.label$just)==length(this.label$pfm),
                                      this.label$just[[idx]], this.label$just[1])
                  this_rot <- ifelse(length(this.label$rot)==length(this.label$pfm),
                                     this.label$rot[[idx]], this.label$rot[1])
                  this_font <- ifelse(length(this.label$font)==length(this.label$pfm),
                                      this.label$font[[idx]], this.label$font[1])
                  this_fontface <- ifelse(length(this.label$fontface)==length(this.label$pfm),
                                          this.label$fontface[[idx]], this.label$fontface[1])
                  this_ic.scale <- ifelse(length(this.label$ic.scale)==length(this.label$pfm),
                                          this.label$ic.scale[[idx]], this.label$ic.scale[1])
                  pushViewport(viewport(x=this.label$x[[idx]],
                                        y=this.height + feature.height,  
                                        just = this_just,
                                        width =  convertWidth(
                                          stringWidth(paste(rep("A",
                                                                ncol(this.label$pfm[[idx]]@mat)),
                                                            collapse = "")), 
                                                     unitTo="npc",
                                                     valueOnly=FALSE),
                                        height = LINEW * this_cex,
                                        angle = this_rot,
                                        default.units = this.label$default.units))
                  motifStack::plotMotifLogoA(pfm = this.label$pfm[[idx]],
                                 font=this_font,
                                 fontface = this_fontface,
                                 ic.scale = this_ic.scale)
                  popViewport()
                }
              }
            }else{
              if(any(is.na(this.label$label))){
                ## grid.text will not plot if first element is empty
                this.label <- filterThisLabel(this.label)
              }
              if(length(this.label$label)>0){
                grid.text(x=this.label$x, y=this.height + feature.height, 
                          label = this.label$label,  
                          just = this.label$just, 
                          hjust = this.label$hjust,
                          vjust = this.label$vjust,
                          rot=this.label$rot,
                          check.overlap = this.label$check.overlap,
                          default.units = this.label$default.units,
                          gp=this.label$gp)
              }
            }
        }
    }
    popViewport()
}

plotLegend <- function(legend, this.height, LINEH){
    ypos <- this.height
    pch <- 21
    if(length(legend)>0){
        if(is.list(legend)){
            thisLabels <- legend[["labels"]]
            if(is.null(legend$gp)){
              gp <- legend[!names(legend) %in% formalArgs(legendGrob)]
              class(gp) <- "gpar"
              legend$gp <- gp
            }else{
              gp <- legend$gp
            }
            if(is.null(gp$cex)) gp$cex <- 1
        }else{
          thisLabels <- names(legend)
          gp <- gpar(fill=legend, cex=1)
          legend <- list(
            labels = thisLabels,
            gp = gp
          )
        }
        if(length(thisLabels)>0){ 
          if(is.null(legend$byrow)) legend$byrow <- TRUE
          if(is.null(legend$vgap)) legend$vgap <- unit(.1*gp$cex[1], "lines")
          if(is.null(legend$hgap)) legend$hgap <- unit(.5*gp$cex[1], "lines")
          if(is.null(legend$ncol)){
            legend$ncol <- getColNum(thisLabels, cex=gp$cex)
          }
          ncol <- legend$ncol
          if(is.null(legend$pch)) legend$pch <- pch
          legend <- legend[names(legend) %in% formalArgs(legendGrob)]
          topblank <- ceiling(length(thisLabels) / ncol) * gp$cex[1]
          pushViewport(viewport(x=.5, 
                                y=ypos+(topblank+.2*gp$cex[1])*LINEH/2, 
                                width=1,
                                height=topblank*LINEH,
                                just="bottom"))
          this.height <- ypos + (topblank+.2*gp$cex[1])*LINEH
          do.call(grid.legend, legend)
          popViewport()
        }
    }
    this.height + LINEH
}