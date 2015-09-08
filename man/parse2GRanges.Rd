\name{parse2GRanges}
\alias{parse2GRanges}
\title{
parse text into GRanges
}
\description{
parse text like "chr13:99,443,451-99,848,821:-" into GRanges
}
\usage{
parse2GRanges(text)
}
\arguments{
  \item{text}{
character vector like "chr13:99,443,451-99,848,821:-" or 
"chr13:99,443,451-99,848,821"
}
}
\value{
an object of \link[GenomicRanges]{GRanges}
}
\author{
Jianhong Ou
}
\examples{
    parse2GRanges("chr13:99,443,451-99,848,821:-")
}
\keyword{ misc }