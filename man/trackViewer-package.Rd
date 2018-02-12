\name{trackViewer-package}
\alias{trackViewer-package}
\alias{trackViewer}
\docType{package}
\title{
A R/Bioconductor package for drawing elegant interactive tracks or lollipop plot to facilitate integrated analysis of multi-omics data
}
\description{
A package that plot data and annotation information along genomic 
coordinates in an elegance style. This tool is based on Gviz but want to draw 
figures in minimal style for publication.
}
\details{
\tabular{ll}{
Package: \tab trackViewer\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-10-18\cr
License: \tab Artistic-2.0\cr
}
This package is minimal designed to plot figure for publication. 
}
\author{
Jianhong Ou, Julie Lihua Zhu

Maintainer: Jianhong Ou <jianhong.ou@umassmed.edu>
}
\keyword{ package }
\examples{
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    library(org.Hs.eg.db)
    trs <- geneModelFromTxdb(TxDb.Hsapiens.UCSC.hg19.knownGene,
                             org.Hs.eg.db,
                             chrom="chr11", 
                             start=122929275, 
                             end=122930122)
    extdata <- system.file("extdata", package="trackViewer",
                       mustWork=TRUE)
    repA <- importScore(paste(extdata, "cpsf160.repA_+.wig", sep="/"),
                        paste(extdata, "cpsf160.repA_-.wig", sep="/"),
                        format="WIG")
    strand(repA@dat) <- "+"
    strand(repA@dat2) <- "-"
    fox2 <- importScore(paste(extdata, "fox2.bed", sep="/"), format="BED")
    dat <- coverageGR(fox2@dat)
    fox2@dat <- dat[strand(dat)=="+"]
    fox2@dat2 <- dat[strand(dat)=="-"]
    gr <- GRanges("chr11", IRanges(122929275, 122930122), strand="-")
    vp <- viewTracks(trackList(repA, fox2, trs), gr=gr, autoOptimizeStyle=TRUE)
    addGuideLine(c(122929767, 122929969), vp=vp)
    addArrowMark(list(x=unit(.5, "npc"), 
                      y=unit(.39, "npc")),
                 col="blue")
}
