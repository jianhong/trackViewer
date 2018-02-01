#' Class \code{"trackViewerStyle"}
#' @description An object of class \code{"trackViewerStyle"} 
#'              represents track viewer style.
#' @aliases trackViewerStyle
#' @rdname trackViewerStyle-class
#' @slot margin \code{"numeric"}, specify the bottom, left, top and right margin.
#' @slot xlas \code{"numeric"}, label direction of x-axis mark. It should 
#' be a integer 0-3. See \code{\link[graphics]{par}:las}
#' @slot xgp A \code{"list"}, object, It will convert to an object of 
#' class \code{\link[grid]{gpar}}. This is basically a list of graphical 
#' parameter settings of x-axis. For y-axis, see \code{\link{yaxisStyle}}
#' @slot xaxis \code{"logical"}, draw x-axis or not
#' @slot autolas \code{"logical"} automatic determine y label direction
#' @slot flip \code{"logical"} flip the x-axis or not, default FALSE
#' @import methods
#' @exportClass trackViewerStyle
#' @examples
#' tvs <- trackViewerStyle()
#' setTrackViewerStyleParam(tvs, "xaxis", TRUE)
#' 
setClass("trackViewerStyle", 
         representation(
             margin="numeric",
             xlas="numeric",
             xgp="list",
             xaxis="logical",
             autolas="logical",
             flip="logical"),
         prototype(
             margin=c(.01, .05, 0, 0),
             xlas=0,
             xgp=list(),
             xaxis=FALSE,
             autolas=FALSE,
             flip=FALSE
             ),
         validity=function(object){
             if(!(object@xlas %in% 0:3))
                 return("xlas should be numeric in {0,1,2,3}. See ?par")
             if(any(object@margin<0 | object@margin>.8))
                 return("margin could not greater than .8 or smaller than 0")
             return(TRUE)
         }
)

#' @rdname trackViewerStyle-class
#' @param \dots Each argument in \dots becomes an slot in the new trackViewerStyle.
#' @export

trackViewerStyle <- function(...){
    new("trackViewerStyle", ...)
}

#' @rdname trackViewerStyle-class
#' @param tvs An object of \code{trackViewerStyle}.
#' @param attr the name of slot to be changed.
#' @param value values to be assigned.
#' @exportMethod setTrackViewerStyleParam
#' @aliases setTrackViewerStyleParam
#' @aliases setTrackViewerStyleParam,trackViewerStyle,character-method
#' 
setGeneric("setTrackViewerStyleParam", function(tvs, attr, value) 
    standardGeneric("setTrackViewerStyleParam"))
#' @rdname trackViewerStyle-class
#' @aliases setTrackViewerStyleParam,trackViewerStyle,character,ANY-method
setMethod("setTrackViewerStyleParam", 
          signature(tvs="trackViewerStyle", attr="character", value="ANY"),
          function(tvs, attr, value){
              if(!attr %in% c("margin", "xlas", "xgp", "xaxis", "autolas", "flip"))
                  stop("attr must be a slot name of trackViewerStyle")
              x <- tvs
              slot(x, attr, check = TRUE) <- value
              eval.parent(substitute(tvs <- x))
              return(invisible(x))
          })

#' Class \code{"pos"}
#' @description An object of class \code{"pos"} represents a point location
#' @rdname pos-class
#' @aliases pos
#' @slot x A \code{\link{numeric}} value, indicates the x position
#' @slot y A \code{\link{numeric}} value, indicates the y position
#' @slot unit \code{"character"} apecifying the units for the corresponding
#' numeric values. See \code{\link[grid]{unit}}
#' @exportClass pos
#' 
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

#' Class \code{"xscale"}
#' @description An object of class \code{"xscale"} represents x-scale style.
#' @rdname xscale-class
#' @aliases xscale
#' @slot from A \code{\link{pos}} class, indicates the start point 
#' postion of x-scale.
#' @slot to A \code{\link{pos}} class, indicates the end point 
#' postion of x-scale.
#' @slot label \code{"character"} the label of x-scale
#' @slot gp A \code{"list"} object, It will convert to an object of 
#' class \code{\link[grid]{gpar}}. This is basically a list of graphical 
#' parameter settings of x-scale.
#' @slot draw A \code{"logical"} value indicating whether the x-scale
#' should be draw.
#' @exportClass xscale

setClass("xscale",
         representation(from="pos",
                        to="pos",
                        label="character",
                        gp="list",
                        draw="logical"),
         prototype(draw=FALSE, gp=list())
         )

#' Class \code{"yaxisStyle"}
#' @description An object of class \code{"yaxisStyle"} represents y-axis style.
#' @rdname yaxisStyle-class
#' @aliases yaxisStyle
#' @slot at \code{"numeric"} vector of y-value locations for the tick marks
#' @slot label \code{"logical"} value indicating whether to draw the 
#' labels on the tick marks.
#' @slot gp A \code{"list"} object, It will convert to an object of 
#' class \code{\link[grid]{gpar}}. This is basically a list of graphical 
#' parameter settings of y-axis.
#' @slot draw A \code{"logical"} value indicating whether the y-axis
#' should be draw.
#' @slot main A \code{"logical"} value indicating whether the y-axis
#' should be draw in left (TRUE) or right (FALSE).
#' @exportClass yaxisStyle

setClass("yaxisStyle",
         representation(at="numeric",
                        label="logical",
                        gp="list",
                        draw="logical",
                        main="logical"),
         prototype(draw=TRUE, label=FALSE, gp=list(), main=TRUE)
         )

#' Class \code{"trackStyle"}
#' @description An object of class \code{"trackStyle"} represents track style.
#' @rdname trackStyle-class
#' @aliases trackStyle
#' @slot tracktype \code{"character"} track type, could be peak or cluster. 
#' Default is "peak". "cluster" is not supported yet.
#' #' @slot color \code{"character"} track color. If the track has dat and dat2 slot,
#' it should have two values.
#' @slot height \code{"numeric"} track height. It should be a value between 0 and 1
#' @slot marginTop \code{"numeric"} track top margin
#' @slot marginBottom \code{"numeric"} track bottom margin
#' @slot xscale object of \code{\link{xscale}}, describe the details of x-scale
#' @slot yaxis object of \code{\link{yaxisStyle}}, describe the details of y-axis
#' @slot ylim \code{"numeric"} y-axis range
#' @slot ylabpos \code{"character"}, ylable postion, ylabpos should 
#' be 'left', 'right', 'topleft', 'bottomleft', 'topright' or 'bottomright'.
#' For gene type track, it also could be 'upstream' or 'downstream'
#' @slot ylablas \code{"numeric"} y lable direction. It should 
#' be a integer 0-3. See \code{\link[graphics]{par}:las}
#' @slot ylabgp A \code{"list"} object, It will convert to an object of 
#' class \code{\link[grid]{gpar}}. This is basically a list of graphical 
#' parameter settings of y-label.
#' @exportClass trackStyle
#'

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
             if(!object@ylabpos %in% c("left", "right", "topleft", "bottomleft", "topright", "bottomright", "upstream", "downstream"))
                 return("ylabpos should be 'left', 'right', 'topleft', 'bottomleft', 'topright' or 'bottomright', 'upstream', 'downstream'.")
             if(!(object@ylablas %in% 0:3))
                return("ylas should be numeric in {0,1,2,3}. See ?par")
             if(!object@tracktype %in% c("peak", "cluster"))
                 return("tracktype must be on of peak or cluster")
             return(TRUE)
         }
)

#' Class \code{"track"}
#' @description An object of class \code{"track"} represents scores of a given track.
#' @rdname trackStyle-class
#' @aliases track
#' @slot dat Object of class \code{\link[GenomicRanges:GRanges-class]{GRanges}}
#' the scores of a given track. It should contain score metadata.
#' @slot dat2 Object of class \code{\link[GenomicRanges:GRanges-class]{GRanges}}
#' the scores of a given track. It should contain score metadata. When dat2
#' and dat is paired, dat will be drawn as positive value where dat2 will be 
#' drawn as negative value (-1 * score)
#' @slot type The type of track. It could be 'data', 'gene', 'transcript' or 'lollipopData'.
#' @slot format The format of the input. It could be "BED", "bedGraph",
#' "WIG", "BigWig" or "BAM"
#' @slot style Object of class \code{\link{trackStyle}}
#' @slot name unused yet
#' @exportClass track
#' @examples 
#' extdata <- system.file("extdata", package="trackViewer",
#' mustWork=TRUE)
#' fox2 <- importScore(file.path(extdata, "fox2.bed"), format="BED")
#' setTrackStyleParam(fox2, "color", c("red","green"))
#' setTrackXscaleParam(fox2, "gp", list(cex=.5))
#' setTrackYaxisParam(fox2, "gp", list(col="blue"))
#' fox2$dat <- GRanges(score=numeric(0))
#' @seealso Please try to use \code{\link{importScore}} and \code{\link{importBam}} to 
#' generate the object.
setClass("track", representation(dat="GRanges",
                                 dat2="GRanges",
                                 type="character",
                                 format="character",
                                 style="trackStyle",
                                 name="character"),
         validity=function(object){
             if(!object@type %in% c("data", "gene", "transcript", "lollipopData"))
                 return("type must be 'data', 'transcript', 'gene', 'lollipopData'")
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
               if(object@type=="lollipopData"){
                 if(is.null(object@dat$score))
                   return("dat should contain score metadata.")
                 if(!all(width(object@dat)==1)){
                   return("Width for lollipopData must be 1")
                 }
                 if(length(object@dat2)>0){
                   if(is.null(object@dat$score))
                     return("dat2 should contain score metadata.")
                   if(!all(width(object@dat2)==1)){
                     return("Width for lollipopData must be 1")
                   }
                 }
               }else{
                 if(is.null(mcols(object@dat)$feature))
                   return("The metadata of dat must contain colnumn 'feature'") 
                 if(length(object@dat2)>0){
                   if(is.null(object@dat$score))
                     return("dat2 should contain score metadata.")
                   if(!all(width(object@dat2)==1)){
                     return("Width for lollipop data must be 1")
                   }
                 }
               }
             }
             return(TRUE)
         })

#' @rdname trackStyle-class
#' @param object an object of trackStyle.
#' @exportMethod show
#' 
#' @aliases show,track-method
setMethod("show", "track", function(object){
    cat("This is an object of track\n", "slot name:", object@name, "\n", 
        "slot type:", object@type, "\n", "slot format:", object@format, "\n")
    cat("slot dat:\n")
    show(object@dat)
    cat("slot dat2:\n")
    show(object@dat2)
    cat("slot style: try object$style to see details.\n")
})
#' Method $
#' @rdname trackStyle-class
#' @param x an object of trackStyle
#' @param name slot name of trackStyle
#' @exportMethod $
#' @aliases $,track-method
setMethod("$", "track", function(x, name) slot(x, name))
#' Method $<-
#' @rdname trackStyle-class
#' @exportMethod $<-
#' @aliases $<-,track-method
setReplaceMethod("$", "track", 
                 function(x, name, value){
                     slot(x, name, check = TRUE) <- value
                     x
                 })
#' Method setTrackStyleParam
#' @rdname trackStyle-class
#' @param ts An object of \code{track}.
#' @param attr the name of slot of \code{\link{trackStyle}} object to be changed.
#' @param value values to be assigned.
#' @exportMethod setTrackStyleParam
#' @aliases setTrackStyleParam
#' @aliases setTrackStyleParam,track,character-method
setGeneric("setTrackStyleParam", function(ts, attr, value) 
    standardGeneric("setTrackStyleParam"))
#' 
#' @rdname trackStyle-class
#' @aliases setTrackStyleParam,track,character,ANY-method
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
#' Method setTrackXscaleParam
#' @rdname trackStyle-class
#' @exportMethod setTrackXscaleParam
#' @aliases setTrackXscaleParam
#' @aliases setTrackXscaleParam,track,character-method
setGeneric("setTrackXscaleParam", function(ts, attr, value) 
    standardGeneric("setTrackXscaleParam"))
#' @rdname trackStyle-class
#' @aliases setTrackXscaleParam,track,character,ANY-method
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
#' Method setTrackYaxisParam
#' @rdname trackStyle-class
#' @exportMethod setTrackYaxisParam
#' @aliases setTrackYaxisParam
#' @aliases setTrackYaxisParam,track,character-method
setGeneric("setTrackYaxisParam", function(ts, attr, value) 
    standardGeneric("setTrackYaxisParam"))
#' @rdname trackStyle-class
#' @aliases setTrackYaxisParam,track,character,ANY-method
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

#' List of tracks
#' @description An extension of List that holds only \code{\link{track}} objects. 
#' @rdname trackList-class
#' @aliases trackList
#' @seealso \code{\link{track}}.
#' @exportClass trackList

setClass("trackList", contains="list", representation(names="vector"),
         validity=function(object){
             re <- sapply(object, class)
             if(any(re!="track"))
                 return("class of elements should be track")
             return(TRUE)
         })

#' @rdname trackList-class
#' @param \dots Each tracks in ... becomes an element in the new 
#' trackList, in the same order. This is analogous to the list constructor, except
#' every argument in ... must be derived from \code{\link{track}}.
#' @param heightDist A vector or NA to define the height of each track.
#' @export trackList
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