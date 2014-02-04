\name{yaxisStyle-class}
\docType{class}
\alias{yaxisStyle-class}
\alias{yaxisStyle}

\title{Class \code{"yaxisStyle"}}
\description{
  An object of class \code{"yaxisStyle"} represents y-axis style.
}

\section{Objects from the Class}{
  Objects can be created by calls of the form 
  \code{new("yaxisStyle", at, label, gp, draw, main)}.
                          
}
\section{Slots}{
  \describe{
    \item{\code{at}}{\code{"numeric"} vector of y-value locations for the tick marks}
    \item{\code{label}}{\code{"logical"} value indicating whether to draw the 
    labels on the tick marks.}
    \item{\code{gp}}{A \code{"list"} object, It will convert to an object of 
    class \code{\link[grid]{gpar}}. This is basically a list of graphical 
    parameter settings of y-axis.}
    \item{\code{draw}}{A \code{"logical"} value indicating whether the y-axis
    should be draw.}
    \item{\code{main}}{A \code{"logical"} value indicating whether the y-axis
    should be draw in left (TRUE) or right (FALSE).}
  }
}
\keyword{classes}
