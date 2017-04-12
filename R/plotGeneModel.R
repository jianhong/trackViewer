##levels of feature "CDS"   "ncRNA" "utr3"  "utr5"
plotGeneModel <- function(track, xscale){
    unit <- 0.25
    y <- 0.5
    transcript <- track@dat
    col <- track@style@color
    strand <- as.character(strand(transcript))[1]
    unitx <- convertWidth(unit(3, "lines"), 
                          unitTo = "npc",
                          valueOnly = TRUE)
    if(unitx > 0.5) unitx <- 0.5
    unitx <- unitx * abs(xscale[2] - xscale[1]) #abs(xscale[2] - xscale[1]) / 10
    ## plot introns
    introns <- gaps(ranges(transcript))
    introns <- introns[start(introns)>min(end(transcript))
                       &end(introns)<max(start(transcript))]
    plotArrow <- function(ax0, ax1, y, col, vp){
        ay0 <- y + unit/2
        ay1 <- y
        ay2 <- y - unit/2
        grid.segments(ax0, ay0, ax1, ay1, 
                      default.units="native", 
                      gp=gpar(col=col))
        grid.segments(ax1, ay1, ax0, ay2, 
                      default.units="native", 
                      gp=gpar(col=col))
    }
    plotIntrons <- function(strand, a, b, y, xscale, col){
        if(a<xscale[1]) a <- xscale[1]
        if(b>xscale[2]) b <- xscale[2]
        w <- b - a
        times <- floor(5*w/unitx)
        for(i in seq_len(times)){
            if(b > a+i*unitx/5+unitx*3/20){
                if(strand=="+") plotArrow(a+i*unitx/5+unitx/10, 
                                          a+i*unitx/5+unitx*3/20,
                                          y, col)
                else {
                    if(strand=="-") 
                        plotArrow(a+i*unitx/5+unitx*3/20, 
                                  a+i*unitx/5+unitx/10,
                                  y, col)
                }
            }
        }
    }
    if(length(introns)>0){
        for(i in 1:length(introns)){
            if(start(introns)[i]<=xscale[2] && end(introns)[i]>=xscale[1]){
                grid.segments(start(introns)[i], y, 
                              end(introns)[i], y, default.units="native", 
                              gp=gpar(col=col))
                plotIntrons(strand, start(introns[i]), 
                            end(introns[i]),
                            y, xscale, col)
            }
        }
    }

    ## plot utr
    utr <- transcript[transcript$feature %in% c("utr3", "utr5")]
    if(length(utr)>0){
        for(i in 1:length(utr)){
            if(start(utr)[i]<=xscale[2] && end(utr)[i]>=xscale[1]){
                grid.rect(start(utr)[i], y,
                          end(utr)[i]-start(utr)[i]+1, unit, 
                          default.units="native", just=c(0, 0.5),
                          gp=gpar(col=NA, fill=col))
            }
        }
    }
    
    ## plot exons
    exons <- transcript[transcript$feature %in% c("CDS", "ncRNA")]
    if(length(exons)>0){
        for(i in 1:length(exons)){
            if(start(exons)[i]<=xscale[2] && end(exons)[i]>=xscale[1]){
                grid.rect(start(exons)[i], 
                          y, 
                          end(exons)[i]-start(exons)[i]+1, 
                          2*unit, default.units="native", just=c(0, 0.5),
                          gp=gpar(col=NA, fill=col))
            }
        }
    }

    return(invisible())
}