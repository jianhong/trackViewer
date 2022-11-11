test_cool<-function(){
  cool <- system.file("extdata", "test.mcool", package = "trackViewer",
                      mustWork=TRUE)
  one <- importGInteractions(file=cool, format="cool",
                      resolution = 2,
                      ranges=GRanges("chr1", IRanges(10, 28)),
                      out = "GInteractions")
  cool <- system.file("extdata", "test.cool", package = "trackViewer",
                      mustWork=TRUE)
  two <- importGInteractions(file=cool, format="cool",
                      ranges=GRanges("chr1", IRanges(10, 28)),
                      out = "GInteractions")
  checkIdentical(one, two)
}