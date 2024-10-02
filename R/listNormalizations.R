#' List the available normalizations
#' @description List the normalizations available in the file.
#' @param file character(1). File name of .hic or .cool/.mcool/.scool
#' @param format character(1). File format, "hic" or "cool".
#' @importFrom rhdf5 H5Fopen h5ls h5read H5close
#' @importFrom strawr readHicNormTypes
#' @export
#' @examples
#' hicfile <- system.file("extdata", "test_chr22.hic", package="trackViewer")
#' listNormalizations(hicfile)
#' coolfile <- system.file("extdata", "test.mcool", package="trackViewer")
#' listNormalizations(coolfile, format="cool")
listNormalizations <- function(file, format=c("hic", "cool")){
  format <- match.arg(format)
  guessFormat(file, format)
  if(format=="hic"){
    readHicNormTypes(fname=file)
  }else{
    coolfile <- checkCoolFile(file)
    if(isMcool(file)){
      coolfileRootName <- coolfileRootName(coolfile)
      resolutions <- listResolutions(file, format = format)
      available_normalization <- lapply(resolutions, function(res){
        mcool_header <- h5ls(coolfile&coolfileRootName&res&"bins",
                             recursive = FALSE)$name
        setdiff(mcool_header, c("chrom", "start", "end"))
      })
      available_normalization <- Reduce(intersect, available_normalization)
    }else{
      mcool_header <- h5ls(coolfile&"bins", recursive = FALSE)$name
      available_normalization <- setdiff(mcool_header,
                                         c("chrom", "start", "end"))
    }
    available_normalization[available_normalization=='weight'] <- 'balanced'
    c('NONE', available_normalization)
  }
}
