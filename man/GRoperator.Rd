\name{GRoperator}
\alias{GRoperator}
\title{
    GRanges operator
}
\description{
    GRanges operations (add, aubtract, multiply, divide)
}
\usage{
GRoperator(A, B, col="score", operator = c("+", "-", "*", "/", "^", "\%\%"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
    an object of GRanges
}
  \item{B}{
    an object of GRanges
}
  \item{col}{
    colname of A and B to be calculated
}
  \item{operator}{
    operator, "+" means A + B, and so on.
}
}
\value{
    an object of GRanges
}
\author{
    Jianhong Ou
}
\examples{
    gr2 <- GRanges(seqnames=c("chr1", "chr1"),
               ranges=IRanges(c(7,13), width=3),
               strand=c("-", "-"), score=3:4)
    gr3 <- GRanges(seqnames=c("chr1", "chr1"),
               ranges=IRanges(c(1, 4), c(3, 9)),
               strand=c("-", "-"), score=c(6L, 2L))
    GRoperator(gr2, gr3, col="score", operator="+")
    GRoperator(gr2, gr3, col="score", operator="-")
    GRoperator(gr2, gr3, col="score", operator="*")
    GRoperator(gr2, gr3, col="score", operator="/")
}

\keyword{ misc }