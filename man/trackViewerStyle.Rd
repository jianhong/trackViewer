\name{trackViewerStyle-class}
\docType{class}
\alias{trackViewerStyle-class}
\alias{trackViewerStyle}
\alias{setTrackViewerStyleParam}
\alias{setTrackViewerStyleParam,trackViewerStyle,character-method}
\alias{setTrackViewerStyleParam,trackViewerStyle,character,ANY-method}

\title{Class \code{"trackViewerStyle"}}
\description{
  An object of class \code{"trackViewerStyle"} represents track viewer style.
}

\section{Objects from the Class}{
  Objects can be created by calls of the form 
  \code{new("trackViewerStyle", margin, xlas, 
                                xgp, xaxis, autolas)}.
                          
}
\section{constructor}{
  trackViewerStyle(\dots): Each argument in \dots becomes an slot in the new 
  trackViewerStyle.
}
\section{Slots}{
  \describe{
    \item{\code{margin}}{\code{"numeric"}, specify the bottom, left, top and right
    margin.}
    \item{\code{xlas}}{\code{"numeric"}, label direction of x-axis mark. It should 
    be a integer 0-3. See \code{\link[graphics]{par}:las}}
    \item{\code{xgp}}{A \code{"list"} object, It will convert to an object of 
    class \code{\link[grid]{gpar}}. This is basically a list of graphical 
    parameter settings of x-axis. For y-axis, see \code{\link{yaxisStyle}}}
    \item{\code{xaxis}}{\code{"logical"}, draw x-axis or not}
    \item{\code{autolas}}{\code{"logical"} automatic determine y label direction}
  }
}

\usage{
    \S4method{setTrackViewerStyleParam}{trackViewerStyle,character,ANY}(tvs, attr, value)
}

\arguments{
    \item{tvs}{An object of \code{trackViewerStyle}.}
    \item{attr}{the name of slot to be changed.}
    \item{value}{values to be assigned.}
}
\section{Methods}{
    \describe{
        \item{setTrackViewerStyleParam}{change the slot values of an object of 
            \code{trackViewerStyle}}
    }
}
\examples{
    tvs <- trackViewerStyle()
    setTrackViewerStyleParam(tvs, "xaxis", TRUE)
}
\keyword{classes}
