#' @importFrom methods as
# Explicit import of something from Matrix to appease automated checks:
#' @importFrom Matrix as.matrix

fm_as_dgCMatrix <- function(x) {
  if (inherits(x, "dgCMatrix")) {
    x
  } else {
    as(as(as(x, "dMatrix"), "generalMatrix"), "CsparseMatrix")
  }
}

fm_as_dgTMatrix <- function(x, unique = TRUE) {
  if (unique) {
    as(fm_as_dgCMatrix(x), "TsparseMatrix")
  } else {
    if (inherits(x, "dgTMatrix")) {
      x
    } else {
      as(as(as(x, "dMatrix"), "generalMatrix"), "TsparseMatrix")
    }
  }
}

#' Sparse partial inverse
#'
#' @param A A symmetric matrix
#'
#' @export

fm_qinv <- function(A) {
  A_C <- fm_as_dgCMatrix(A)
  stopifnot(nrow(A_C) == ncol(A_C))
  if (!identical(A_C, Matrix::t(A_C))) {
    warning("Asymmetric matrix A detected, but only lower left triangle will be used.")
  }
  C_qinv(A_C)
}
