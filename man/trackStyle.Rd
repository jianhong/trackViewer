\name{trackStyle-class}
\docType{class}
\alias{trackStyle-class}
\alias{trackStyle}

\title{Class \code{"trackStyle"}}
\description{
  An object of class \code{"trackStyle"} represents track style.
}

\section{Objects from the Class}{
  Objects can be created by calls of the form 
  \code{new("trackStyle", tracktype, color, height, marginTop, marginBottom,
                        xscale, yaxis, ylim, ylabpos, ylablas, ylabgp)}.
                          
}
\section{Slots}{
  \describe{
    \item{\code{tracktype}}{\code{"character"} track type, could be peak or cluster. 
    Default is "peak". "cluster" is not supported yet.}
    \item{\code{color}}{\code{"character"} track color. If the track has dat and dat2 slot,
    it should have two values.}
    \item{\code{height}}{\code{"numeric"} track height. It should be a value between 0 and 1}
    \item{\code{marginTop}}{\code{"numeric"} track top margin}
    \item{\code{marginBottom}}{\code{"numeric"} track bottom margin}
    \item{\code{xscale}}{object of \code{\link{xscale}}, describe the details of x-scale}
    \item{\code{yaxis}}{object of \code{\link{yaxisStyle}}, describe the details of y-axis}
    \item{\code{ylim}}{\code{"numeric"} y-axis range}
    \item{\code{ylabpos}}{\code{"character"}, ylable postion, ylabpos should 
    be 'left', 'right', 'topleft', 'bottomleft', 'topright' or 'bottomright'.
    For gene type track, it also could be 'upstream' or 'downstream'}
    \item{\code{ylablas}}{\code{"numeric"} y lable direction. It should 
    be a integer 0-3. See \code{\link[graphics]{par}:las}}
    \item{\code{ylabgp}}{A \code{"list"} object, It will convert to an object of 
    class \code{\link[grid]{gpar}}. This is basically a list of graphical 
    parameter settings of y-label.}
  }
}
\keyword{classes}
