getLabelLeghtRate <- function(SNPs, ratio.yx){
  if(length(SNPs$label.parameter.rot)>0) {
    labels.rot <- SNPs$label.parameter.rot
  }else{
    labels.rot <- 90
  }
  labels.cex <- 1
  if(length(SNPs$label.parameter.gp)>0){
    if(length(SNPs$label.parameter.gp$cex)>0)
      labels.cex <- SNPs$label.parameter.gp$cex[[1]][1]
  }
  labels.length.rate <- labels.cex * 
    max(cospi((labels.rot-90)/180), 0) * ratio.yx
  return(labels.length.rate)
}
getLabelsMaxHeight <- function(SNPs, ratio.yx){
  labels.length.rate <- getLabelLeghtRate(SNPs, ratio.yx)
  stringH <- as.numeric(convertY(stringHeight("W"), "npc"))
  if(length(names(SNPs))>0){
    maxStrHeight <- 
      max(as.numeric(
        convertX(stringWidth(names(SNPs)), "npc")
      ))+stringH
  }else{
    maxStrHeight <- 0
  }
  maxStrHeight <- maxStrHeight * labels.length.rate
}