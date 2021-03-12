#' plot tracks for single cell RNAseq
#' @description Plot single cell RNAseq data as heatmap track for Seurat object.
#' @param object Seurat object.
#' @param files bam file to be scanned.
#' @param samplenames sample names for files.
#' @param ... parameters used by \link[GenomicAlignments:readGAlignments]{readGAlignmentsList} 
#' or \link[GenomicAlignments]{readGAlignments}
#' @param txdb TxDb object for gene model.
#' @param gene Gene name to plot. (row value)
#' @param id The id of gene used in txdb.
#' @param idents indentity class to define the groups to plot. (column value)
#' @param gr GRanges object to define the ploting region.
#' @param color vector of colors used in heatmap.  
#' @param withCoverageTrack plot coverage track or not.
#' @param flag An integer(2) vector used to filter reads based on their 
#' 'flag' entry. 
#' @importFrom GenomeInfoDb seqlevelsStyle
#' @export
#' @examples 
#' gene <- "Cd8a"
#' 
#' 

importScRNAseqScore <- function(object, files, samplenames, ...,
                         txdb, gene, id, idents, gr, color,
                         withCoverageTrack=TRUE,
                         flag=scanBamFlag(isSecondaryAlignment = FALSE,
                                          isUnmappedQuery=FALSE,
                                          isNotPassingQualityControls = FALSE,
                                          isSupplementaryAlignment = FALSE)){
  idents <- as.character(idents)
  object <- subset(object, idents = idents, features=gene)
  stopifnot(length(files)==length(samplenames))
  if(nrow(object)!=1){
    stop("No feature or more than one feature detected.")
  }
  if(ncol(object)<1){
    stop("Not valid idents.")
  }
  id <- id[gene==rownames(object)]
  gene <- rownames(object)
  what <- character(0)
  #https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/bam
  # BC: Sample index read
  # CB: Chromium cellular barcode sequence that is error-corrected and confirmed against a list of known-good barcode sequences.
  # xf: Extra alignment flags. The bit flags can be interpreted as follows: 1 - The read is confidently mapped to a feature; 2 - The read maps to a feature that the majority of other reads with this UMI did not; 8 - This read is representative for the molecule and can be treated as a UMI count. Bits 4, 16 and 32 are used internally by 10X.
  tag <- c("BC", "CB", "xf")
  tr <- geneTrack(ids = id, txdb = txdb, symbols = gene)[[1]]
  if(missing(gr)){
    gr <- range(tr$dat)
    w10 <- floor(width(gr)/10)
    start(gr) <- max(1, start(gr)-w10)
    end(gr) <- end(gr)+w10
  }
  which <- gr
  seqlens <- scanBamHeader(files[1])[[1]]$target
  seq_gr <- GRanges(names(seqlens), IRanges(rep(1, length(seqlens)), seqlens))
  seqlevelsStyle(which) <- seqlevelsStyle(seq_gr)[1]
  param <-
    ScanBamParam(what=what,
                 tag=tag,
                 which=which,
                 flag=flag)
  
  reads <- 
    lapply(files, function(bamFile) 
      readGAlignments(bamFile, ..., param=param))
  names(reads) <- samplenames
  ## filter by xf
  filterByXf <- any(vapply(reads, function(.ele){
    xf <- mcols(.ele)$xf
    length(xf)==length(.ele) && !all(xf==0)
  }, FUN.VALUE = TRUE))
  if(filterByXf){
    reads <- lapply(reads, function(.ele){
      xf <- mcols(.ele)$xf
      if(length(xf)==length(.ele)){
      .ele <- .ele[bitwAnd(xf, 8)==8]
      }else{
        .ele <- .ele[rep(FALSE, length(.ele))]
      }
      .ele
    })
  }
  
  
  ## filter by cell barcode
  cellNames <- colnames(object)
  reads <- lapply(reads, function(.ele) 
    .ele[sub("-", "_", mcols(.ele)$CB) %in% cellNames])
  
  ## split the reads by cell and get coverage
  cvg <- lapply(reads, function(.ele){
    lapply(split(.ele, mcols(.ele)$CB), coverage)
  })
  ## convert coverage to Views
  chr <- as.character(seqnames(which)[1])
  cvg <- lapply(cvg, function(.ele){
    lapply(.ele, function(.e){
      Views(.e[[chr]], ranges(which))
    })
  })
  ## convert Views to dataframe
  cvg <- lapply(cvg, function(.ele){
    .ele <- lapply(.ele, function(.e) 
      viewApply(.e, FUN=as.numeric, simplify = TRUE))
    .ele <- lapply(.ele, as.numeric)
    do.call(rbind, .ele)
  })
  ## recover the data frame for all the cells
  mat <- matrix(0, nrow = ncol(object), ncol = width(which))
  rownames(mat) <- cellNames
  cvg <- lapply(cvg, function(.ele){
    newmat <- mat
    if(!is.null(dim(.ele))){
      newmat[sub("-", "_", rownames(.ele)), ] <- .ele
    }
    newmat
  })
  ## covert it back to GRanges
  cvg <- lapply(cvg, function(.ele){
    .ele <- .ele[order(rowSums(.ele), decreasing = TRUE), ]
    .ele <- apply(.ele, 1, Rle)
    .ele <- lapply(.ele, function(.e){
      l <- runLength(.e)
      l <- cumsum(l)
      shift(IRanges(start = c(1, l[-length(l)]+1), end = l,
                    score=runValue(.e)),
            shift = start(which))
    })
    .ele <- GRanges(chr, unlist(IRangesList(.ele)))
  })
  ## split the cvg by the idents
  cvg <- lapply(cvg, function(.ele){
    .ele$ident <- object@active.ident[names(.ele)] #Idents
    split(.ele, .ele$ident)
  })
  ## switch list
  cvg <- lapply(idents, function(.ele){
    lapply(cvg, function(.e) .e[[.ele]])
  })
  names(cvg) <- idents
  cvg <- unlist(cvg)
  trs <- lapply(cvg, function(.ele){
    seqlevelsStyle(.ele) <- seqlevelsStyle(tr$dat)[1]
    tr1 <- new("track", dat=.ele, type="scRNA", format="BAM")
    if(withCoverageTrack){
      .cvg <- coverage(.ele, weight = .ele$score)
      tr2 <- new("track", dat=.cvg, type="data", format="BAM")
      list(perCell=tr1, coverage=tr2)
    }else{
      tr1
    }
  })
  trs <- unlist(trs)
  trackList <- trackList(trs, tr, heightDist = c(length(cvg), 1/2))
  return(trackList)
}
