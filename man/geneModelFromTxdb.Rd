\name{geneModelFromTxdb}
\alias{geneModelFromTxdb}
\title{
Prepare gene model from an object of TxDb
}
\description{
Generate an object of \code{\link{track}} for \code{\link{viewTracks}} by given 
parameters.
}
\usage{
geneModelFromTxdb(txdb, chrom, start, end, strand=c("*", "+", "-"), txdump=NULL)
}
\arguments{
  \item{txdb}{An object of \code{\link[GenomicFeatures]{TxDb}}}
  \item{chrom}{chromosome name, must be a seqname of txdb}
  \item{start}{start position}
  \item{end}{end position}
  \item{strand}{strand}
  \item{txdump}{output of as.list(txdb), a list of data frames that can be used 
  to make the db again with no loss of information.}
}
\value{
An object of \code{\link{track}}
}
\author{
Jianhong Ou
}

\seealso{
See Also as \code{\link{importScore}}, \code{\link{importBam}}, \code{\link{viewTracks}}
}
\examples{
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    trs <- geneModelFromTxdb(TxDb.Hsapiens.UCSC.hg19.knownGene,
                             "chr20", 22560000, 22565000, "-")
}
\keyword{ importData }
