\name{pos-class}
\docType{class}
\alias{pos-class}
\alias{pos}

\title{Class \code{"pos"}}
\description{
  An object of class \code{"pos"} represents a point location
}

\section{Objects from the Class}{
  Objects can be created by calls of the form 
  \code{new("pos", x, y, unit)}.
                          
}
\section{Slots}{
  \describe{
    \item{\code{x}}{A \code{\link{numeric}} value, indicates the x position}
    \item{\code{y}}{A \code{\link{numeric}} value, indicates the y position}
    \item{\code{unit}}{\code{"character"} apecifying the units for the corresponding
    numeric values. See \code{\link[grid]{unit}}}
  }
}
\keyword{classes}
