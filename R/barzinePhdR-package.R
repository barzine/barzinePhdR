#' barzinePhdR
#'
#' This package contains all the require R functions supporting the analyses
#' in [https://github.com/barzine/phd-analyses](https://github.com/barzine/phd-analyses)
#' which are presented in [https://github.com/barzine/thesis](https://github.com/barzine/thesis)
#'
#'
#' @importFrom reshape2 melt dcast
#' @importFrom stats t.test phyper setNames na.omit hclust relevel as.dist as.formula binom.test pbinom sd aggregate
#' @importFrom utils combn read.table
#' @importFrom grid textGrob gpar gTree grid.text pushViewport popViewport
#' @importFrom grDevices rgb colorRampPalette cairo_pdf dev.off
#' @importFrom gplots heatmap.2
#' @importFrom gridExtra grid.arrange arrangeGrob
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggthemes geom_rangeframe
#' @importFrom WGCNA cor
#' @export
#'
initialise<-function(){
  sharedEnv<-new.env(parent = emptyenv())
  assign("epsilon",1e-10,envir = sharedEnv)
  assign("debug",FALSE,envir = sharedEnv)
  #assign('colCbc',colCbc,envir=sharedEnv)
  suppressWarnings(utils::data("List"))
  #assign("sharedEnv",sharedEnv,envir = globalenv())
}

#utils::globalVariables("sharedEnv", package="barzinePhdR")
utils::globalVariables("sharedEnv")

#fix-me: use rlang::.data or something less ugly than globalVariables
utils::globalVariables(c("..count..", "..density..")) # for ggplot2
utils::globalVariables(c("Group", "Tissue", "Comparison", "nb.tissues",
                         "variable","value","Value","Type","Correlation")) # for ggplot2
utils::globalVariables(c("Bins", "Cumulp","Comparison")) # for Cumulative_spe
utils::globalVariables("Label") # for  bibarplotsDiversityCond
utils::globalVariables("Genes") # for cross.sharedDistrib_firstLast
utils::globalVariables(c("Rank","TSpercent")) # for cumulSpeSimulAggregated
utils::globalVariables(c('Ranked.cor','Cor')) # for sortedGenesCorr
utils::globalVariables("Cut") # for evolCorrCumul

