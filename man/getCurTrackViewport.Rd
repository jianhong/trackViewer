\name{getCurTrackViewport}
\alias{getCurTrackViewport}
\title{
Get current track viewport
}
\description{
Get current track viewport for addGuideLine
}
\usage{
getCurTrackViewport(curViewerStyle, start, end)
}
\arguments{
  \item{curViewerStyle}{an object of \code{\link{trackViewerStyle}}}
  \item{start}{start position of current track}
  \item{end}{end position of current track}
}
\value{
an object of \code{\link[grid]{viewport}}
}
\author{
Jianhong Ou
}
\seealso{
See Also as \code{\link{addGuideLine}}
}
\examples{
    vp <- getCurTrackViewport(trackViewerStyle(), 10000, 10200)
    addGuideLine(c(10010, 10025, 10150), vp=vp)
}
\keyword{ misc }
