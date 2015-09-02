\name{importData}
\alias{importData}
\title{
Reading data from a BED or WIG file to RleList
}
\description{
Read a \code{\link{track}} object from a BED, bedGraph, WIG or BigWig file
to RleList
}
\usage{
importData(files, format=NA, ranges=GRanges())
}
\arguments{
  \item{files}{The path to the files to read.}
  \item{format}{
    The format of import file. Could be BAM, BED, bedGraph, WIG or BigWig}
  \item{ranges}{An object of \code{\link[GenomicRanges]{GRanges}} to indicate
  the range to be imported}
}
\value{
a list of \code{\link[IRanges]{RleList}}. 
}
\author{
Jianhong Ou
}
\seealso{
See Also as \code{\link{interactiveViewer}}
}
\examples{
#import a BED file
bedfile <- system.file("tests", "test.bed", package="rtracklayer",
                       mustWork=TRUE)
dat <- importData(files=bedfile, format="BED",
                    ranges=GRanges("chr7", IRanges(127471197, 127474697)))

##import a WIG file
wigfile <- system.file("tests", "step.wig", package = "rtracklayer",
                       mustWork=TRUE)
dat <- importData(files=wigfile, format="WIG", 
                  ranges=GRanges("chr19", 
                                 IRanges(59104701, 59110920)))

##import a BigWig file
if(.Platform$OS.type!="windows"){
    ##this is because we are using rtracklayer::import
    bwfile <- system.file("tests", "test.bw", package = "rtracklayer",
                       mustWork=TRUE)
    dat <- importData(files=bwfile, format="BigWig", 
                      ranges=GRanges("chr19", IRanges(1500, 2700)))
}
}
\keyword{ importData }
