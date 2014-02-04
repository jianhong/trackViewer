\name{track-class}
\docType{class}
\alias{track-class}
\alias{track}
\alias{setTrackStyleParam}
\alias{setTrackStyleParam,track,character-method}
\alias{setTrackStyleParam,track,character,ANY-method}
\alias{setTrackXscaleParam}
\alias{setTrackXscaleParam,track,character-method}
\alias{setTrackXscaleParam,track,character,ANY-method}
\alias{setTrackYaxisParam}
\alias{setTrackYaxisParam,track,character-method}
\alias{setTrackYaxisParam,track,character,ANY-method}
\alias{$,track-method}
\alias{$<-,track-method}

\title{Class \code{"track"}}
\description{
  An object of class \code{"track"} represents scores of a given track.
}

\section{Objects from the Class}{
  Objects can be created by calls of the form 
  \code{new("track", dat, dat2, type, format, style, name)}.
}
\section{Slots}{
  \describe{
    \item{\code{dat}}{Object of class \code{\link[GenomicRanges]{GRanges}}
    the scores of a given track. It should contain score metadata.}
    \item{\code{dat2}}{Object of class \code{\link[GenomicRanges]{GRanges}}
    the scores of a given track. It should contain score metadata. When dat2
    and dat is paired, dat will be drawn as positive value where dat2 will be 
    drawn as negative value (-1 * score)}
    \item{\code{type}}{The type of track. It could be 'data' or 'gene'.}
    \item{\code{format}}{The format of the input. It could be "BED", "bedGraph", 
    "WIG", "BigWig" or "BAM"}
    \item{\code{style}}{Object of class \code{\link{trackStyle}}}
    \item{\code{name}}{unused yet}
  }
}
\usage{
    \S4method{setTrackStyleParam}{track,character,ANY}(ts, attr, value)
    \S4method{setTrackXscaleParam}{track,character,ANY}(ts, attr, value) 
    \S4method{setTrackYaxisParam}{track,character,ANY}(ts, attr, value)
}

\arguments{
    \item{ts}{An object of \code{track}.}
    \item{attr}{the name of slot of \code{\link{trackStyle}} object to be changed.}
    \item{value}{values to be assigned.}
}
\section{Methods}{
    \describe{
        \item{setTrackStyleParam}{change the slot values of \code{\link{trackStyle}} object for an object of \code{track}}
        \item{setTrackXscaleParam}{change the \code{\link{xscale}} slot 
        values for an object of \code{track}}
        \item{setTrackYaxisParam}{change the \code{\link{yaxisStyle}} 
        values for an object of \code{track}}
        \item{$, $<-}{Get or set the slot of \code{\link{track}}}
    }
}
\examples{
    extdata <- system.file("extdata", package="trackViewer",
                       mustWork=TRUE)
    fox2 <- importScore(file.path(extdata, "fox2.bed"), format="BED")
    setTrackStyleParam(fox2, "color", c("red","green"))
    setTrackXscaleParam(fox2, "gp", list(cex=.5))
    setTrackYaxisParam(fox2, "gp", list(col="blue"))
    fox2$dat <- GRanges(score=numeric(0))
}
\section{See Also}{
    Please try to use \code{\link{importScore}} and \code{\link{importBam}} to 
    generate the object.
}
\keyword{classes}
