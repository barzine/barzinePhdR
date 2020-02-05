#' Exclude the genes (rows) from a data.frame based on a given list
#'
#' @param DF data.frame to be filtered
#' @param genelist a vector of ID (e.g. Ensembl ID for mitochondria genes) to be excluded from DF
#'
#' @return data.frame without the filtered out genes from the specified list
#' @export
#'
excludeGenes<-function(DF,genelist){
  DF[!(rownames(DF) %in% genelist),]
}


#' Convert R/FPKM to TPM
#'
#' @param DF numeric data.frame
#'
#' @return numeric data.frame
#' @export
#'
rpkm2tpm<-function(DF){
  vColSums<-colSums(DF)
  DF1<-(DF/vColSums)*10^6
  return(DF1)
}

#' Select specifically the genes that are expressed above a threshold
#' for a given number of times in a data.frame
#'
#' @param DF numeric (expression) data.frame
#' @param threshold numeric; default: 0. Minimal expression to be considered
#' @param expBreadth integer; default: 1. Exact number of times to be expressed to be considered
#' @param verbose boolean; default: TRUE. Print the number of rows of the returned data.frame
#'
#' @return a data.frame that contains genes (ie rows) that
#'        are expressed exactly a number of times (expBreadth) above the given threshold
#' @export
#'

selectSpecific<-function(DF,threshold=0,expBreadth=1,verbose=TRUE){
  if(threshold){
    DF.nb <- rowSums( DF  %>=%  threshold)
  }else{
    DF.nb <- rowSums( DF >  0)
  }
  DF.spe<-DF[DF.nb==expBreadth,]
  if(verbose) print(nrow(DF.spe))
  return(DF.spe)
}





