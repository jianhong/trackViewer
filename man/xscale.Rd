\name{xscale-class}
\docType{class}
\alias{xscale-class}
\alias{xscale}

\title{Class \code{"xscale"}}
\description{
  An object of class \code{"xscale"} represents x-scale style.
}

\section{Objects from the Class}{
  Objects can be created by calls of the form 
  \code{new("xscale", from, to, label, gp, draw)}.
                          
}
\section{Slots}{
  \describe{
    \item{\code{from}}{A \code{\link{pos}} class, indicates the start point 
    postion of x-scale.}
    \item{\code{to}}{A \code{\link{pos}} class, indicates the end point 
    postion of x-scale.}
    \item{\code{label}}{\code{"character"} the label of x-scale}
    \item{\code{gp}}{A \code{"list"} object, It will convert to an object of 
    class \code{\link[grid]{gpar}}. This is basically a list of graphical 
    parameter settings of x-scale.}
    \item{\code{draw}}{A \code{"logical"} value indicating whether the x-scale
    should be draw.}
  }
}
\keyword{classes}
