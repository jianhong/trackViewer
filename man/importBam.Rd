\name{importBam}
\alias{importBam}
\title{
Reading data from a BAM file
}
\description{
Read a \code{\link{track}} object from a BAM file
}
\usage{
importBam(file, file2, ranges=GRanges(), pairs=FALSE)
}
\arguments{
  \item{file}{The path to the BAM file to read.}
  \item{file2}{The path to the second BAM file to read.}
  \item{ranges}{An object of \code{\link[GenomicRanges]{GRanges}} to indicate
  the range to be imported}
  \item{pairs}{logical object to indicate the BAM is paired or not. See 
  \code{\link[GenomicAlignments]{readGAlignmentsFromBam}}}
}
\value{
a \code{\link{track}} object 
}
\author{
Jianhong Ou
}
\seealso{
See Also as \code{\link{importScore}}, \code{\link{track}}, \code{\link{viewTracks}}
}
\examples{
bamfile <- system.file("extdata", "ex1.bam", package="Rsamtools",
                       mustWork=TRUE)
dat <- importBam(file=bamfile, ranges=GRanges("seq1", IRanges(1, 50), strand="+"))
}
\keyword{ importData }