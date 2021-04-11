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
#' @importFrom GenomeInfoDb seqlevelsStyle Seqinfo isCircular genome
#' @examples 
#' \dontrun{
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' test_file <- "https://github.com/10XGenomics/subset-bam/raw/master/test/test.bam"
#' trs <- importScSeqScore(files=test_file, 
#'                         txdb=TxDb.Hsapiens.UCSC.hg19.knownGene,
#'                         id="653635", gene = "WASH7P")
#' }
#' 

importScSeqScore <- function(object, files, samplenames, ...,
                         txdb, gene, id, idents, gr, color,
                         withCoverageTrack=TRUE,
                         flag=scanBamFlag(isSecondaryAlignment = FALSE,
                                          isUnmappedQuery=FALSE,
                                          isNotPassingQualityControls = FALSE,
                                          isSupplementaryAlignment = FALSE)){
  if(missing(txdb)){
    stop("txdb is required!")
  }
  if(missing(gr)){
    if(missing(id) || missing(gene)){
      stop("id and gene is required.")
    }
    ## get gr from id or gene
    trs <- geneTrack(ids=id, txdb = txdb, symbols = gene)
    gr <- lapply(trs, function(.ele) range(unname(.ele@dat)))
    gr <- unlist(GRangesList(unname(gr)))
  }else{
    if(missing(object)){
      if(missing(id) || missing(gene)){
        ## get id and gene from gr
        trs <- geneModelFromTxdb(txdb=txdb, gr=gr)
        ids <- unlist(lapply(trs, function(.ele) unique(.ele@dat$gene)))
        ids <- unique(ids)
        if(missing(id)) id <- ids
        if(missing(gene)) gene <- id
      }
    }
  }
  if(!missing(object)){
    idents <- as.character(idents)
    object <- subset(object, idents = idents, features=gene)
    if(nrow(object)!=1){
      stop("No feature or more than one feature detected.")
    }
    if(ncol(object)<1){
      stop("Not valid idents.")
    }
    if(missing(id) || missing(gene) || missing(txdb)){
      stop("txdb, id and gene is required.")
    }
    id <- id[gene==rownames(object)]
    gene <- rownames(object)
  }
  tr <- geneTrack(ids = id, txdb = txdb, symbols = gene)[[1]]
  if(!missing(samplenames)){
    stopifnot(length(files)==length(samplenames))
  }else{
    samplenames <- sub(".bam", "", basename(files), ignore.case = TRUE)
  }
  what <- character(0)
  #https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/bam
  # BC: Sample index read
  # CB: Chromium cellular barcode sequence that is error-corrected and confirmed against a list of known-good barcode sequences.
  # xf: Extra alignment flags. The bit flags can be interpreted as follows: 1 - The read is confidently mapped to a feature; 2 - The read maps to a feature that the majority of other reads with this UMI did not; 8 - This read is representative for the molecule and can be treated as a UMI count. Bits 4, 16 and 32 are used internally by 10X.
  tag <- c("BC", "CB", "xf")
  seqlens <- scanBamHeader(files[1])[[1]]$target
  seq_gr <- GRanges(names(seqlens), IRanges(rep(1, length(seqlens)), seqlens))
  if(missing(gr)){
    if(length(tr)){
      gr <- range(tr$dat)
      w10 <- floor(width(gr)/10)
      start(gr) <- max(1, start(gr)-w10)
      end(gr) <- end(gr)+w10
    }else{
      gr <- seq_gr
    }
  }
  which <- gr
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
    if(length(xf)){
      .ele <- .ele[!is.na(xf)]
      xf <- mcols(.ele)$xf
    }else{
      return(FALSE)
    }
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
  if(!missing(object)){
    cellNames <- colnames(object)
    reads <- lapply(reads, function(.ele) 
      .ele[sub("-", "_", mcols(.ele)$CB) %in% cellNames])
  }
  
  ## split the reads by cell and get coverage
  cvg <- lapply(reads, function(.ele){
    .e <- .ele
    sep <- "__SEQN__"
    CB <- unique(mcols(.ele)$CB)
    seqn <- paste(rep(CB, each=length(seqlevels(.ele))),
                  rep(seqlevels(.ele), length(CB)),
                  sep=sep)
    seqinf <- Seqinfo(seqnames = c(seqlevels(.ele), seqn),
                      seqlengths = unname(rep(seqlengths(.ele), length(CB)+1)),
                      isCircular = unname(rep(isCircular(.ele), length(CB)+1)),
                      genome= unname(rep(genome(.ele), length(CB)+1)))
    seqname <- paste(mcols(.e)$CB, as.character(seqnames(.e)),
                     sep=sep)
    seqinfo(.e) <- seqinf[c(seqlevels(.ele), seqn[seqn %in% seqname])]
    seqnames(.e) <- factor(seqname, levels = seqlevels(.e))
    .cvg <- coverage(.e)
    .cvg <- split(.cvg, 
                  sub(paste0("^(.*?)", sep, ".*$"), "\\1", names(.cvg)))
    .cvg <- .cvg[names(.cvg) %in% CB]
    .cvg <- lapply(.cvg, function(.e){
      names(.e) <- sub(paste0("^.*?", sep), "", names(.e))
      .e
    })
    .cvg
  })
  if(!missing(object)){
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
  }else{
    ## split the cvg by the idents
    idents <- unique(unlist(lapply(cvg, function(.ele) names(.ele))))
    cvg <- lapply(cvg, function(.ele){
      .ele <- lapply(.ele[idents], function(.e){
        .e <- as(.e, "GRanges")
        .e[.e$score>0]
        })
    })
  }
  
  ## switch list
  cvg <- lapply(idents, function(.ele){
    lapply(cvg, function(.e) .e[[.ele]])
  })
  names(cvg) <- idents
  cvg <- unlist(cvg)
  trs <- lapply(cvg, function(.ele){
    seqlevelsStyle(.ele) <- seqlevelsStyle(tr$dat)[1]
    tr1 <- new("track", dat=.ele, type="scSeq", format="BAM")
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
