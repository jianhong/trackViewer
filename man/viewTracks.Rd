\name{viewTracks}
\alias{viewTracks}
\title{
    plot the tracks
}
\description{
    A function to plot the data for given range
}
\usage{
viewTracks(trackList, chromosome, start, end, strand, gr=GRanges(),
           ignore.strand=TRUE,
           viewerStyle=trackViewerStyle(), autoOptimizeStyle=FALSE,
           newpage=TRUE)
}
\arguments{
  \item{trackList}{an object of \code{\link{trackList}}}
  \item{chromosome}{chromosome}
  \item{start}{start position}
  \item{end}{end position}
  \item{strand}{strand}
  \item{gr}{an object of \code{\link[GenomicRanges]{GRanges}}}
  \item{ignore.strand}{ignore the strand or not when do filter. default TRUE}
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
    extdata <- system.file("extdata", package="trackViewer",
                       mustWork=TRUE)
    files <- dir(extdata, "-.wig")
    tracks <- lapply(paste(extdata, files, sep="/"), 
                        importScore, format="WIG")
    tracks <- lapply(tracks, function(.ele) {strand(.ele@dat) <- "-"; .ele})
    fox2 <- importScore(paste(extdata, "fox2.bed", sep="/"), format="BED")
    dat <- coverageGR(fox2@dat)
    fox2@dat <- dat[strand(dat)=="+"]
    fox2@dat2 <- dat[strand(dat)=="-"]
    gr <- GRanges("chr11", IRanges(122929275, 122930122), strand="-")
    viewTracks(trackList(tracks, fox2), gr=gr, autoOptimizeStyle=TRUE)
}
\keyword{ plot }
