.onLoad <- function(libname, pkgname) {
  initialise()
}

.onAttach<-function(libname, pkgname){
  packageStartupMessage("initialise WGCNA::allowWGCNAThreads() to allow for multi-threading")
}
