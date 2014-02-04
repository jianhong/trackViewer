\name{addArrowMark}
\alias{addArrowMark}
\title{
Add arrow mark to the figure at a given position
}
\description{
A function to add arrow mark for emphasizing peaks 
}
\usage{
addArrowMark(pos=grid.locator(), angle=15, 
             length=unit(.25, "inches"), col="red",
             type="closed")
}

\arguments{
  \item{pos}{A unit object representing the location of arrow mark to be placed 
  at current viewport. Default is the value of grid.locator, which will get the 
  location of the mouse click.}
  \item{angle}{A parameter passed into grid::arrow function. The angle of arrow 
  head in degrees (smaller numbers produce narrower, pointier arrows). Essentially 
  describes the width of the arrow head.}
  \item{length}{A parameter passed into grid::arrow function. Aunit specifying 
  the length of the arrow head.}
  \item{col}{color of the arrow}
  \item{type}{A parameter passed into grid::arrow function. One of "open" or "closed"
  indicating whether the arrow head should be a closed triangle.}
}
\value{
NULL
}
\author{
Jianhong Ou
}
\seealso{
See Also as \code{\link{addGuideLine}}, \code{\link[grid:arrow]{arrow}}
}
\examples{
    grid.newpage()
    addArrowMark(list(x=unit(.5, "npc"), 
                      y=unit(.5, "npc")), 
                 col="blue")
}
\keyword{ ~misc }
