#' convert GInteractions to track object
#' @description Convert GInteractions object to track object
#' @param gi an object of GInteractions
#' @param gi2 an object of GInteractions
#' @return an track object
#' @importFrom S4Vectors first second
#' @export
#' @examples 
#' gi <- readRDS(system.file("extdata", "nij.chr6.51120000.53200000.gi.rds", package="trackViewer"))
#' gi2track(gi)
gi2track <- function(gi, gi2){
  stopifnot(is(gi, "GInteractions"))
  .oneGI2track <- function(.gi){
    a <- first(.gi)
    if(length(.gi$score)==length(a)){
      a$score <- .gi$score
    }else{
      a$score <- 1
    }
    return(new("track", dat=a, dat2=second(.gi),
               type="interactionData", format="BED"))
  }
  .gi2gr <- function(.gi){
    a <- first(.gi)
    a$target <- second(.gi)
    if(length(.gi$score)==length(a)){
      a$score <- .gi$score
    }else{
      a$score <- 1
    }
    return(a)
  }
  if(!missing(gi2)){
    if(length(gi2)>0){
      if(length(gi)==0){
        return(.oneGI2track(gi2))
      }else{
        ## put the second coordinates to target metadata
        return(new("track", dat=.gi2gr(gi),
                   dat2=.gi2gr(gi2),
                   type="interactionData",
                   format="BED"))
      }
    }
  }
  if(length(gi)==0){
    return(new("track", dat=GRanges(), dat2=GRanges(),
               type="interactionData", format="BED"))
  }
  return(.oneGI2track(gi))
}