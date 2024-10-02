#' List the available chromosome
#' @description List the chromosomes available in the file.
#' @param file character(1). File name of .hic or .cool/.mcool/.scool
#' @param format character(1). File format, "hic" or "cool".
#' @importFrom rhdf5 H5Fopen h5ls H5close
#' @importFrom strawr readHicChroms
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
    readHicChroms(fname=file)
  }else{
    coolfile <- checkCoolFile(file)
    if(isMcool(file)){
      coolfileRootName <- coolfileRootName(coolfile)
      res <- h5ls(coolfile&coolfileRootName, recursive = FALSE)$name
      res <- res[which.max(as.numeric(res))]
      chroms <- h5read(coolfile,
                       coolfilePathName(coolfileRootName, res, "chroms"))$name
    }else{
      chroms <- h5read(coolfile, "chroms")$name
    }
    chroms
  }
}
