#' Wrapper around colorRampalette with two prechosen colours.
#'
#' @param Colvec vector of two string of characters.
#'               Default: c('aliceblue','darkcyan')
#'
#' @return a function
#' @export
#'
colCbc<-function(Colvec=c('aliceblue','darkcyan')){
  colorRampPalette(Colvec)
}

