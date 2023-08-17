#' @importFrom grDevices dev.off pdf
#' @importFrom graphics plot.new
#' @importFrom methods is
#' @importFrom stats quantile runif
#' @importFrom utils str
#' @import methods

.onLoad <- function(libname, pkgname) {
  options(list("Matrix.quiet" = TRUE))
}


.onUnload <- function(libpath) {
  library.dynam.unload("fmesher", libpath)
}
