test_cool<-function(){
  cool <- system.file("extdata", "test.mcool", package = "trackViewer",
                      mustWork=TRUE)
  r1 <- listResolutions(cool, format='cool')
  n1 <- listNormalizations(cool, format='cool')
  c1 <- listChromosomes(cool, format='cool')
  one <- importGInteractions(file=cool, format="cool",
                      resolution = 2,
                      ranges=GRanges("chr1", IRanges(10, 28)),
                      out = "GInteractions")
  cool <- system.file("extdata", "test.cool", package = "trackViewer",
                      mustWork=TRUE)
  r2 <- listResolutions(cool, format='cool')
  n2 <- listNormalizations(cool, format='cool')
  c2 <- listChromosomes(cool, format='cool')
  two <- importGInteractions(file=cool, format="cool",
                      ranges=GRanges("chr1", IRanges(10, 28)),
                      out = "GInteractions")
  checkIdentical(one, two)
}