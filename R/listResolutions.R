#' List the available resolutions
#' @description List the resolutions available in the file.
#' @param file character(1). File name of .hic or .cool/.mcool/.scool
#' @param format character(1). File format, "hic" or "cool".
#' @importFrom rhdf5 H5Fopen h5ls H5close
#' @importFrom strawr readHicBpResolutions
#' @export
#' @examples
#' hicfile <- system.file("extdata", "test_chr22.hic", package="trackViewer")
#' listResolutions(hicfile)
#' coolfile <- system.file("extdata", "test.mcool", package="trackViewer")
#' listResolutions(coolfile, format="cool")
listResolutions <- function(file, format=c("hic", "cool")){
  format <- match.arg(format)
  guessFormat(file, format)
  if(format=="hic"){
    readHicBpResolutions(fname=file)
  }else{
    if(isMcool(file)){
      coolfile <- checkCoolFile(file)
      coolfileRootName <- coolfileRootName(coolfile)
      h5ls(coolfile&coolfileRootName, recursive = FALSE)$name
    }else{
      return(NULL)
    }
  }
}
