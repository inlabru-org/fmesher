#' @importFrom R.methodsS3 setMethodS3
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics plot.new
#' @importFrom methods is
#' @importFrom stats quantile runif
#' @importFrom utils str
#' @import Matrix

.onLoad <- function(libname, pkgname) {
    options(list("Matrix.quiet"=TRUE))
}
