\name{trackList-class}
\docType{class}
\alias{trackList-class}
\alias{trackList}

\title{List of tracks}
\description{
  An extension of List that holds only \code{\link{track}} objects. 
}

\section{constructor}{
  trackList(\dots, heightDist=NA): Each tracks in ... becomes an element in the new 
  trackList, in the same order. This is analogous to the list constructor, except
  every argument in ... must be derived from \code{\link{track}}. The heightDist
  is vector or NA to define the height of each track.
}
\section{See Also}{
    \code{\link{track}}.
}
\keyword{classes}
