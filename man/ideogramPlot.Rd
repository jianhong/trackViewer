\name{ideogramPlot}
\alias{ideogramPlot}
\title{
plot ideogram with data
}
\description{
plot ideogram with data for multiple chromosomes 
}
\usage{
ideogramPlot(ideo, dataList, layout=NULL,
                         horiz=TRUE,
                         parameterList=
                           list(vp=plotViewport(margins=c(.1, 4.1, .3, .1)),
                                ideoHeight=unit(.5, "npc"), 
                                vgap=unit(.3, "lines"),
                                ylabs="auto", 
                                ylabsRot=ifelse(horiz, 0, 90), 
                                ylabsPos=unit(2.5, "lines"), 
                                xaxis=FALSE, yaxis=FALSE,
                                xlab="", 
                                types="barplot", heights=NULL, 
                                dataColumn="score", 
                                gps=gpar(col="black", fill="gray")), 
                         gp=gpar(fill=NA, lwd=2), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ideo}{
output of \link{loadIdeogram}.
}
  \item{dataList}{
 a \link[GenomicRanges]{GRangesList} of data to plot.
}
  \item{layout}{
 The layout of chromosomes. Could be a list with chromosome names
as its elements.
}
  \item{horiz}{
 a logical value. If FALSE, the ideograms are drawn vertically 
to the left. If TRUE, the ideograms are drawn horizontally at the bottom.
}
  \item{parameterList}{
parameterList a list of parameters for each dataset in the dataList. 
The elements of the parameters could be xlabs, ylabs, etc. type could be
barplot, line, point, heatmap.
}
  \item{gp}{
 parameters used for \link[grid]{grid.roundrect}.
}
  \item{...}{
 parameters not used.
}
}
\author{
Jianhong Ou
}
\seealso{
See Also as \code{\link{loadIdeogram}}
}
\examples{
\dontrun{
ideo <- loadIdeogram("hg38")
library(rtracklayer)
library(grid)
dataList <- ideo
dataList$score <- as.numeric(dataList$gieStain)
dataList <- dataList[dataList$gieStain!="gneg"]
dataList <- GRangesList(dataList)
grid.newpage()
ideogramPlot(ideo, dataList, 
             layout=list("chr1", "chr2", c("chr3", "chr22"), 
                         c("chr4", "chr21"), c("chr5", "chr20"), 
                         c("chr6", "chr19"), c("chr7", "chr18"),
                         c("chr8", "chr17"), c("chr9", "chr16"),
                         c("chr10", "chr15"), c("chr11", "chr14"),
                         c("chr12", "chr13"), c("chrX", "chrY")),
             parameterList = list(types="heatmap"))
}
}

\keyword{ plot }
