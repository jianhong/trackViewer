setClass("trackViewerStyle", 
         representation(
             margin="numeric",
             xlas="numeric",
             xgp="list",
             xaxis="logical",
             autolas="logical"),
         prototype(
             margin=c(.01, .05, 0, 0),
             xlas=0,
             xgp=list(),
             xaxis=FALSE,
             autolas=FALSE
             ),
         validity=function(object){
             if(!(object@xlas %in% 0:3))
                 return("xlas should be numeric in {0,1,2,3}. See ?par")
             if(any(object@margin<0 | object@margin>.8))
                 return("margin could not greater than .8 or smaller than 0")
             return(TRUE)
         }
)

trackViewerStyle <- function(...){
    new("trackViewerStyle", ...)
}

setGeneric("setTrackViewerStyleParam", function(tvs, attr, value) 
    standardGeneric("setTrackViewerStyleParam"))
setMethod("setTrackViewerStyleParam", 
          signature(tvs="trackViewerStyle", attr="character", value="ANY"),
          function(tvs, attr, value){
              if(!attr %in% c("margin", "xlas", "xgp", "xaxis", "autolas"))
                  stop("attr must be a slot name of trackViewerStyle")
              x <- tvs
              slot(x, attr, check = TRUE) <- value
              eval.parent(substitute(tvs <- x))
              return(invisible(x))
          })

setClass("pos", representation(x="numeric", y="numeric", unit="character"),
         prototype(x=0.5, y=0.5, unit="npc"),
         validity=function(object){
             if(!object@unit %in% c("npc", "cm", "inches", "mm", "points", "picas",
                                    "bigpts", "dida", "cicero", "scaledpts",
                                    "lines", "char", "native", "snpc",
                                    "strwidth", "strheight", "grobwidth",
                                    "grobheight"))
                 return("unit must be units of grid unit. See ?unit")
             return(TRUE)
         })

setClass("xscale",
         representation(from="pos",
                        to="pos",
                        label="character",
                        gp="list",
                        draw="logical"),
         prototype(draw=FALSE, gp=list())
         )

setClass("yaxisStyle",
         representation(at="numeric",
                        label="logical",
                        gp="list",
                        draw="logical",
                        main="logical"),
         prototype(draw=TRUE, label=FALSE, gp=list(), main=TRUE)
         )

setClass("trackStyle",
         representation(tracktype="character",
                        color="character",
                        height="numeric",
                        marginTop="numeric",
                        marginBottom="numeric",
                        xscale="xscale",
                        yaxis="yaxisStyle",
                        ylim="numeric",
                        ylabpos="character",
                        ylablas="numeric",
                        ylabgp="list"
                        ),
         prototype(
             marginTop=0,
             marginBottom=0.05,
             color=c("black","black"),
             tracktype="peak",
             ylabpos="left",
             ylablas=0,
             ylabgp=list()
         ),
         validity=function(object){
             if(!object@ylabpos %in% c("left", "right", "topleft", "bottomleft", "topright", "bottomright"))
                 return("ylabpos should be 'left', 'right', 'topleft', 'bottomleft', 'topright' or 'bottomright'.")
             if(!(object@ylablas %in% 0:3))
                return("ylas should be numeric in {0,1,2,3}. See ?par")
             if(!object@tracktype %in% c("peak", "cluster"))
                 return("tracktype must be on of peak or cluster")
             return(TRUE)
         }
)

setClass("track", representation(dat="GRanges",
                                 dat2="GRanges",
                                 type="character",
                                 format="character",
                                 style="trackStyle",
                                 name="character"),
         validity=function(object){
             if(!object@type %in% c("data", "gene"))
                 return("type must be 'data' or 'gene'")
             if(object@type=="data"){
                 if(!object@format %in% c("BED", "bedGraph", "WIG", "BigWig", "BAM"))
                     return("format must be one of \"BED\", 
                            \"bedGraph\", \"WIG\", \"BigWig\"")
                 if(is.null(object@dat$score))
                     return("dat should contain score metadata.")
                 if(length(object@dat2)>0){
                     if(is.null(object@dat$score))
                         return("dat2 should contain score metadata.")
                 }
                 if(object@format!="WIG"){
                     if(!inherits(object@dat$score, c("numeric", "integer")))
                         return("class of score metadata should be numeric")
                 }else{
                     if(class(object@dat$score)!="CompressedCharacterList")
                         return("Please try ?imortScore for WIG files")
                 }
             }else{
                 if(is.null(mcols(object@dat)$feature))
                     return("The metadata of dat must contain colnumn 'feature'") 
             }
             return(TRUE)
         })

setMethod("$", "track", function(x, name) slot(x, name))
setReplaceMethod("$", "track", 
                 function(x, name, value){
                     slot(x, name, check = TRUE) <- value
                     x
                 })

setGeneric("setTrackStyleParam", function(ts, attr, value) 
    standardGeneric("setTrackStyleParam"))
setMethod("setTrackStyleParam", 
          signature(ts="track", attr="character", value="ANY"),
          function(ts, attr, value){
              if(!attr %in% c("tracktype", "color", "height", "marginTop", 
                              "marginBottom", "ylim", "ylabpos", "ylablas",
                              "ylabgp"))
                  stop("attr must be a slot name (except xscale and yaxis) of trackStyle.
                       try setTrackXscaleParam for xscale slot and 
                       setTrackYaxisParam for yaxis.")
              x <- ts
              slot(x@style, attr, check = TRUE) <- value
              eval.parent(substitute(ts <- x))
              return(invisible(x))
          })

setGeneric("setTrackXscaleParam", function(ts, attr, value) 
    standardGeneric("setTrackXscaleParam"))
setMethod("setTrackXscaleParam", 
          signature(ts="track", attr="character", value="ANY"),
          function(ts, attr, value){
              if(!attr %in% c("from", "to", "label", "gp", "draw"))
                  stop("attr must be a slot name of xscale object")
              x <- ts
              slot(x@style@xscale, attr, check = TRUE) <- value
              eval.parent(substitute(ts <- x))
              return(invisible(x))
          })

setGeneric("setTrackYaxisParam", function(ts, attr, value) 
    standardGeneric("setTrackYaxisParam"))
setMethod("setTrackYaxisParam", 
          signature(ts="track", attr="character", value="ANY"),
          function(ts, attr, value){
              if(!attr %in% c("at", "label", "gp", "draw", "main"))
                  stop("attr must be a slot name of xscale object")
              x <- ts
              slot(x@style@yaxis, attr, check = TRUE) <- value
              eval.parent(substitute(ts <- x))
              return(invisible(x))
          })

setClass("trackList", contains="list", representation(names="vector"),
         validity=function(object){
             re <- sapply(object, class)
             if(any(re!="track"))
                 return("class of elements should be track")
             return(TRUE)
         })
trackList <- function(..., heightDist=NA){
    listData <- list(...)
    dots <- substitute(list(...))[-1]
    names <- as.character(sapply(dots, deparse))
    if(is.na(heightDist[1])){
        heightDist <- rep(1, length(listData))
    }else{
        if(length(heightDist)!=length(listData)){
            stop("length of heightDist should be same as length of inputs")
        }
    }
    heightDist <- heightDist * 1/sum(heightDist)
    tmpData <- list()
    recursiveList <- function(tmp, tmpData, name, heightDist){
        if(is.list(tmp)){
            if(length(tmp)==0) return(tmpData)
            for(w in 1:length(tmp)) {
                each <- tmp[[w]]
                curHeightDist <- heightDist/length(tmp)
                if(is.list(each)){
                    tmpData <- recursiveList(each, tmpData, 
                                         names(tmp)[w], curHeightDist)
                }else{
                    .name <- names(tmp)[w]
                    each@style@height <- curHeightDist
                    tmpData <- c(tmpData, each)
                    names(tmpData)[length(tmpData)] <- .name
                }
            }
        }else{
            tmp@style@height <- heightDist
            tmpData <- c(tmpData, tmp)
            names(tmpData)[length(tmpData)] <- name
        }
        tmpData
    }
    if(length(listData)>=1){
        for(i in 1:length(listData)){
            tmp <- listData[[i]]
            tmpData <- recursiveList(tmp, tmpData, names[i], heightDist[i])
            rm("tmp")
        }
    }
    listData <- tmpData
    rm("tmpData")
    
    new("trackList", listData)
}