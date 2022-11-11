# codes for cooler
guessFormat <- function(file, format){
  if(is.character(file)){
    if(!grepl(format, file)){
      warning("The input file is:", file, 
              " But the input format is:", format)
    }
  }
}
checkCoolFile <- function(coolfile){
  stopifnot(length(coolfile)==1)
  stopifnot(is.character(coolfile) || is(coolfile, "H5IdComponent"))
  if(is.character(coolfile)){
    coolfile <- H5Fopen(coolfile, flags="H5F_ACC_RDONLY")
  }
  return(coolfile)
}
isMcool <- function(coolfile){
  coolRootName <- coolfileRootName(coolfile)
  return(!all(c("bins", "chroms", "indexes", "pixels") %in%
           coolRootName))
}

coolfileRootName <- function(coolfile){
  coolfile <- checkCoolFile(coolfile)
  h5ls(coolfile, recursive = FALSE)$name
}

coolfileGetByPath <- function(coolfile, resolution, subPath){
  coolfile <- checkCoolFile(coolfile)
  if(!isMcool(coolfile)){
    res0 <- "coolfile"
  }else{
    available_res <- listResolutions(coolfile, "cool")
    if(!resolution %in% available_res){
      stop(paste("The given resolution is not available in file.",
                 "Available resolutions: ",
                 paste(available_res, collapse = ","), "."
      ))
    }
    coolfileRootName <- coolfileRootName(coolfile)
    res0 <- paste0("coolfile$", coolfileRootName, "$`", resolution, "`")
  }
  if(!missing(subPath)){
    subPath <- subPath[1]
    if(subPath!="" && !is.na(subPath)){
      res0 <- paste0(res0, "$", subPath)
    }
  }
  eval(parse(text=res0))
}

cooler_indexes <- function(coolfile, resolution){
  coolfileGetByPath(coolfile, resolution, "indexes")
}

coolfilePathName <- function(rootName, resolution, subPath){
  stopifnot(is.character(subPath) && length(subPath)==1)
  ifelse(subPath %in% rootName, subPath, 
         paste(rootName, resolution, subPath, sep="/"))
}

cooler_bins <- function(coolfile, resolution, seqname){
  seqs <- coolfileGetByPath(coolfile, resolution, "chroms")
  if(!all(seqname %in% seqs$name)){
    stop(paste("The given seqname is not available in file.",
               "Available seqnames: ",
               paste(seqs, collapse = ","), "."
    ))
  }
  rootName <- coolfileRootName(coolfile)
  indexes <- cooler_indexes(coolfile, resolution)
  ir1 <- IRanges(indexes$chrom_offset[-length(indexes$chrom_offset)]+1,
                 indexes$chrom_offset[-1], 
                 names = seqs$name)
  ir0 <- ir1[seqname]
  name <- coolfilePathName(rootName, resolution, 'bins')
  bins <- list()
  on.exit(h5closeAll())
  for(i in seq_along(ir0)){
    h5closeAll()
    bins[[names(ir0)[i]]] <- 
      as.data.frame(h5read(coolfile, 
                           name=name, 
                           start=start(ir0[i]),
                           count = width(ir0[i])))
    bins[[names(ir0)[i]]]$bin <- seq(start(ir0[i]), end(ir0[i]))-1
  }
  bins <- do.call(rbind, bins)
  gr <- GRanges(seqnames = bins$chrom, 
                IRanges(bins$start+1, bins$end, names = bins$bin))
  if(length(bins$weight)){
    gr$weight <- bins$weight
  }else{
    gr$weight <- 1
  }
  gr
}

cooler_pixels <- function(coolfile, resolution, gr=GRanges()){
  stopifnot(is(gr, "GRanges"))
  stopifnot(length(gr)<3)
  if(length(gr)==0) return(NULL)
  seqname <- as.character(seqnames(gr))
  bins <- cooler_bins(coolfile, resolution, seqname)
  bins1 <- subsetByOverlaps(bins, gr[1])
  if(length(gr)==1){
    bins2 <- bins1
  }else{
    bins2 <- subsetByOverlaps(bins, gr[2])
  }
  if(length(bins)==0) return(NULL)
  indexes <- cooler_indexes(coolfile, resolution)
  bins_index <- IRanges(indexes$bin1_offset[-length(indexes$bin1_offset)]+1,
                        indexes$bin1_offset[-1],
                        names = seq.int(length(indexes$bin1_offset)-1)-1)
  ## read bin1
  bins_index_1 <- bins_index[names(bins1)]
  bins_index_1 <- bins_index_1[width(bins_index_1)>0]
  ir0 <- reduce(bins_index_1)
  rootName <- coolfileRootName(coolfile)
  name <- coolfilePathName(rootName, resolution, "pixels")
  pixels <- list()
  on.exit(h5closeAll())
  for(i in seq_along(ir0)){
    h5closeAll()
    curr_pixels <- 
      as.data.frame(h5read(coolfile, 
                           name=name, 
                           start=start(ir0[i]),
                           count = width(ir0[i])))
    curr_pixels <- 
      curr_pixels[curr_pixels$bin2_id %in% names(bins2),,drop=FALSE]
    pixels[[i]] <- curr_pixels
  }
  pixels <- do.call(rbind, pixels)
  gi <- GInteractions(anchor1 = bins[as.character(pixels$bin1_id)],
                      anchor2 = bins[as.character(pixels$bin2_id)],
                      count = pixels$count)
  gi$balanced <- gi$count/(gi$anchor1.weight*gi$anchor2.weight)
  gi$anchor1.weight <- NULL
  gi$anchor2.weight <- NULL
  gi
}



