test_hic<-function(){
  if(.Platform$OS.type!="windows"){
  hicfile <- system.file("extdata", "testBp.hic", package="trackViewer")
  contactRecords <- straw(
    fname = hicfile,
    norm = "NONE",
    chr1loc = "chr0:0:100000000",
    chr2loc = "chr0:0:100000000",
    binsize = 10000,
    unit = "BP",
    matrix = "observed"
  )
  
	checkEquals(nrow(contactRecords), 2500)
	hicfile <- system.file("extdata", "test_chr22.hic", package="trackViewer")
	contactRecords <- straw(
	  fname = hicfile,
	  norm = "KR",
	  chr1loc = "22:50000000:100000000",
	  chr2loc = "22:50000000:100000000",
	  binsize = 100000,
	  unit = "BP",
	  matrix = "observed"
	)
	
	checkEquals(nrow(contactRecords), 70)
	
	hicfile <- system.file("extdata", "testFrag.hic", package="trackViewer")
	contactRecords <- straw(
	  fname = hicfile,
	  norm = "NONE",
	  chr1loc = "chr0:0:100",
	  chr2loc = "chr0:0:100",
	  binsize = 1,
	  unit = "FRAG",
	  matrix = "observed"
	)
  }
}