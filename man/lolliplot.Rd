\name{lolliplot}
\alias{lolliplot}
\title{
 Lolliplots
}
\description{
 Plot variants and somatic mutations
}
\usage{
lolliplot(SNP.gr, features = NULL, ranges = NULL, 
          type = c("circle", "pie", "pin"), newpage = TRUE)
}

\arguments{
  \item{SNP.gr}{
    A object of \link[GenomicRanges]{GRanges} or \link[GenomicRanges]{GRangesList}.
    All the width of GRanges must be 1.
}
  \item{features}{
    A object of \link[GenomicRanges]{GRanges} or \link[GenomicRanges]{GRangesList}.
}
  \item{ranges}{
    A object of \link[GenomicRanges]{GRanges} or \link[GenomicRanges]{GRangesList}.
}
  \item{type}{
    Character. Could be circle, pie or pin.
}
  \item{newpage}{
  plot in the new page or not.
}
}
\details{
    In SNP.gr and features, metadata of the GRanges object will be used to control the 
    color, fill, border, height, data source of pie if the type is pie.
}
\value{
NULL
}
\author{
Jianhong Ou
}
\examples{
    SNP <- c(10, 100, 105, 108, 400, 410, 420, 600, 700, 805, 840, 1400, 1402)
    x <- sample.int(100, length(SNP))
    SNP.gr <- GRanges("chr1", IRanges(SNP, width=1, names=paste0("snp", SNP)), 
                 value1=x, value2=100-x)
    SNP.gr$color <- rep(list(c("red", 'blue')), length(SNP))
    SNP.gr$border <- sample.int(7, length(SNP), replace=TRUE)
    features <- GRanges("chr1", IRanges(c(1, 501, 1001), 
                                        width=c(120, 500, 405),
                                        names=paste0("block", 1:3)),
                        color="black",
                        fill=c("#FF8833", "#51C6E6", "#DFA32D"),
                        height=c(0.1, 0.05, 0.08))
    lolliplot(SNP.gr, features, type="pie")
}
\keyword{ misc }