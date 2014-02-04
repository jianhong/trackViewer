\name{interactiveViewer}
\alias{interactiveViewer}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    plot the tracks with GUI
}
\description{
    interactively plot the tracks
}
\usage{
interactiveViewer(trackList, chromosome, start, end, strand, 
    viewerStyle = trackViewerStyle(), autoOptimizeStyle = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{trackList}{an object of \code{\link{trackList}}}
  \item{chromosome}{chromosome}
  \item{start}{start position}
  \item{end}{end position}
  \item{strand}{strand}
  \item{viewerStyle}{an object of \code{\link{trackViewerStyle}}}
  \item{autoOptimizeStyle}{should use \code{\link{optimizeStyle}} to optimize style}
}
\author{
Jianhong Ou
}
\seealso{
See Also as \code{\link{viewTracks}}
}
\examples{
if(interactive()){
    interactiveViewer()
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ plot }
