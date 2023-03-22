#' List the available chromosome
#' @description List the chromosomes available in the file.
#' @param file character(1). File name of .hic or .cool/.mcool/.scool
#' @param format character(1). File format, "hic" or "cool".
#' @importFrom rhdf5 H5Fopen h5ls
#' @export
#' @examples
#' hicfile <- system.file("extdata", "test_chr22.hic", package="trackViewer")
#' listChromosomes(hicfile)
#' coolfile <- system.file("extdata", "test.mcool", package="trackViewer")
#' listChromosomes(coolfile, format="cool")
listChromosomes <- function(file, format=c("hic", "cool")){
  format <- match.arg(format)
  guessFormat(file, format)
  if(format=="hic"){
    .Call("_trackViewer_readHicChroms",
          file, PACKAGE = "trackViewer")
  }else{
    coolfile <- checkCoolFile(file)
    coolfileRootName <- coolfileRootName(coolfile)
    res <- h5ls(coolfile&coolfileRootName, recursive = FALSE)$name
    chroms <- lapply(res, function(.ele){
      res0 <- paste0("coolfile$", coolfileRootName, "$`", .ele, "`$chroms$name")
      eval(parse(text=res0))
    })
    sort(unique(unlist(chroms)))
  }
}
