\name{loadIdeogram}
\alias{loadIdeogram}
\title{
load ideogram from UCSC
}
\description{
Download ideogram table from UCSC 
}
\usage{
loadIdeogram(genome, chrom = NULL, ranges = NULL, ...)
}
\arguments{
  \item{genome}{
Assembly name  assigned by UCSC, such as hg38, mm10.
}
  \item{chrom}{
A character vector of chromosome names, or NULL.
}
  \item{ranges}{
A \link[IRanges]{Ranges} object with the intervals.
}
  \item{...}{
Additional arguments to pass to the \link[GenomicRanges]{GRanges} constructor.
}
}
\value{
A \link[GenomicRanges]{GRanges} object.
}
\author{
Jianhong Ou
}
\seealso{
See Also as \code{\link{ideogramPlot}}
}
\examples{
\dontrun{
 head(loadIdeogram("hg38"))
}
}

\keyword{ misc }
