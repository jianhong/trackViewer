\name{plotGRanges}
\alias{plotGRanges}
\title{
    plot GRanges data
}
\description{
    A function to plot GRanges data for given range
}
\usage{
plotGRanges(..., range=GRanges(), 
            viewerStyle=trackViewerStyle(),
            autoOptimizeStyle=FALSE,
            newpage=TRUE)
}
\arguments{
  \item{...}{one or more objects of \code{\link[GenomicRanges]{GRanges}}}
  \item{range}{an object of \code{\link[GenomicRanges]{GRanges}}}
  \item{viewerStyle}{an object of \code{\link{trackViewerStyle}}}
  \item{autoOptimizeStyle}{should use \code{\link{optimizeStyle}} to optimize style}
  \item{newpage}{should be draw on a new page?}
}
\value{
An object of \code{\link[grid]{viewport}} for \code{\link{addGuideLine}}
}
\author{
Jianhong Ou
}
\seealso{
See Also as \code{\link{addGuideLine}}, \code{\link{addArrowMark}}
}
\examples{
    gr1 <- GRanges("chr1", IRanges(1:50, 51:100))
    gr2 <- GRanges("chr1", IRanges(seq(from=10, to=80, by=5),
                                   seq(from=20, to=90, by=5)))
    vp <- plotGRanges(gr1, gr2, range=GRanges("chr1", IRanges(1, 100)))
    addGuideLine(guideLine=c(5, 10, 50, 90), col=2:5, vp=vp)
    
    gr <- GRanges("chr1", IRanges(c(1, 11, 21, 31), width=9), 
                    score=c(5, 10, 5, 1))
    plotGRanges(gr, range=GRanges("chr1", IRanges(1, 50)))
}
\keyword{ plot }
