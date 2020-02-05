.onLoad <- function(libname, pkgname) {
  initialise()
  WGCNA::allowWGCNAThreads()
}
