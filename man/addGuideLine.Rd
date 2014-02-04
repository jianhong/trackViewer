\name{addGuideLine}
\alias{addGuideLine}
\title{
Add guide lines to the tracks
}
\description{
A function to add lines for emphasizing the positions 
}
\usage{
addGuideLine(guideLine, col="gray", lty="dashed", lwd=1, vp=NULL)
}

\arguments{
  \item{guideLine}{The genomic coordinates to draw the lines}
  \item{col}{A vector for the line color}
  \item{lty}{A vector for the line type}
  \item{lwd}{A vector for the line width}
  \item{vp}{A Grid viewport object. It must be output of \code{\link{viewTracks}}}
}
\value{
NULL
}
\author{
Jianhong Ou
}

\seealso{
See Also as \code{\link{getCurTrackViewport}}, \code{\link{addArrowMark}}, 
\code{\link{viewTracks}}
}
\examples{
    vp <- getCurTrackViewport(trackViewerStyle(), 10000, 10200)
    addGuideLine(c(10010, 10025, 10150), vp=vp)
}
\keyword{ misc }
