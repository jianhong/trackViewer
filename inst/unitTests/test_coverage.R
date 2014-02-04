test_coverageGR<-function(){
	gr <- GRanges("chr1", IRanges(c(1:10), width=3), strand="+", score=rep(1, 10))
    gr <- coverageGR(gr)
	checkEquals(c(1, 2, 3, 2, 1), score(gr))
	checkEquals(c(1, 2, 3, 11, 12), start(gr))
	checkEquals(c(1, 2, 10, 11, 12), end(gr))
}