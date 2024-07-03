#' Sparse partial inverse
#'
#' Compute sparse partial matrix inverse. Requires special build of the `fmesher`
#' package.
#'
#' @param A A symmetric positive definite matrix
#'
#' @keywords internal

fm_qinv <- function(A) {
  A_C <- fm_as_dgCMatrix(A)
  stopifnot(nrow(A_C) == ncol(A_C))
  if (!identical(A_C, Matrix::t(A_C))) {
    warning("Asymmetric matrix A detected, but only lower left triangle will be used.")
  }
  fmesher_qinv(A_C)
}
