\name{coverageGR}
\alias{coverageGR}

\title{
calculate coverage
}
\description{
calculate coverage for \code{\link[GenomicRanges]{GRanges}}, 
\code{\link[GenomicAlignments]{GAlignments}} or \code{\link[GenomicAlignments]{GAlignmentPairs}}
}
\usage{
coverageGR(gr)
}
\arguments{
  \item{gr}{an object of RGanges, GAlignments or GAlignmentPairs}
}
\value{
an object of GRanges
}
\author{
Jianhong Ou
}

\seealso{
    See Also as \code{\link[IRanges]{coverage}}, 
    \code{\link[GenomicRanges]{coverage-methods}}
}
\examples{
    bed <- system.file("extdata", "fox2.bed", package="trackViewer",
                       mustWork=TRUE)
    fox2 <- importScore(bed)
    fox2@dat <- coverageGR(fox2@dat)
}

\keyword{ importData }
