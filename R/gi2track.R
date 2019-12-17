#' convert GInteractions to track object
#' @description Convert GInteractions object to track object
#' @param gi an object of GInteractions
#' @return an track object
#' @importFrom S4Vectors first second
#' @export
#' @examples 
#' gi <- readRDS(system.file("extdata", "nij.chr6.51120000.53200000.gi.rds", package="trackViewer"))
#' gi2track(gi)
gi2track <- function(gi){
  stopifnot(is(gi, "GInteractions"))
  a <- first(gi)
  if(length(gi$score)==length(a)){
    a$score <- gi$score
  }else{
    a$score <- 1
  }
  return(new("track", dat=a, dat2=second(gi),
             type="interactionData", format="BED"))
}