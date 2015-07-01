\name{parseWIG}
\alias{parseWIG}
\title{
    convert WIG format track to BED format track
}
\description{
    convert WIG format track to BED format track for a given range
}
\usage{
parseWIG(trackScore, chrom, from, to)
}
\arguments{
  \item{trackScore}{
    an object of track with WIG format
}
  \item{chrom}{
    sequence name of the chromosome
}
  \item{from}{
    start coordinate
}
  \item{to}{
    end coordinate
}
}
\value{
an object of track
}
\author{
    Jianhong Ou
}
\seealso{
    \code{\link{track}}
}
\examples{
    extdata <- system.file("extdata", package="trackViewer", mustWork=TRUE)
    repA <- importScore(file.path(extdata, "cpsf160.repA_-.wig"),
                        file.path(extdata, "cpsf160.repA_+.wig"),
                        format="WIG")
    strand(repA$dat) <- "-"
    strand(repA$dat2) <- "+"
    parseWIG(repA, chrom="chr11", from=122929275, to=122930122)
}
\keyword{ misc }