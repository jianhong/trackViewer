test_hic<-function(){
  if(.Platform$OS.type!="windows"){
  hicfile <- system.file("extdata", "testBp.hic", package="trackViewer")
  contactRecords <- .Call("_trackViewer_getContactRecords",
                          hicfilename=hicfile,
                                      qname1 = "0",
                                      start1 = 0,
                                      end1 = 100000000,
                                      qname2 = "0",
                                      start2 = 0,
                                      end2 = 100000000,
                                      binSize = 10000,
                                      normalization = "NONE",
                                      unit = "BP",
                          PACKAGE = "trackViewer")
	checkEquals(nrow(contactRecords), 2500)
	hicfile <- system.file("extdata", "test_chr22.hic", package="trackViewer")
	contactRecords <- .Call("_trackViewer_getContactRecords",
	                        hicfilename=hicfile,
	                                    qname1 = "chr22",
	                                    start1 = 50000000,
	                                    end1 = 100000000,
	                                    qname2 = "22",
	                                    start2 = 50000000,
	                                    end2 = 100000000,
	                                    binSize = 100000,
	                                    normalization = "KR",
	                                    unit = "BP",
	                        PACKAGE = "trackViewer")
	checkEquals(nrow(contactRecords), 70)
	
	hicfile <- system.file("extdata", "testFrag.hic", package="trackViewer")
	contactRecords <- .Call("_trackViewer_getContactRecords",
	                        hicfilename=hicfile,
	                                    qname1 = "0",
	                                    start1 = 0,
	                                    end1 = 100000000,
	                                    qname2 = "0",
	                                    start2 = 0,
	                                    end2 = 100000000,
	                                    binSize = 1,
	                                    normalization = "NONE",
	                                    unit = "FRAG",
	                        PACKAGE = "trackViewer")
  }
}