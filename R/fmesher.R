#' Sparse partial inverse
#'
#' @param A A symmetric matrix
#'
#' @export

fm_qinv <- function(A) {
  A_C <- as(as(A, "CsparseMatrix"), "dgCMatrix")
  stopifnot(nrow(A_C) == ncol(A_C))
  if (!identical(A_C,t(A_C))) {
    warning("Asymmetric matrix A detected, but only lower left triangle will be used.")
  }
  C_qinv(A_C)
}
