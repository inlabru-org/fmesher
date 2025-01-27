#' Sparse partial inverse
#'
#' Compute sparse partial matrix inverse. Slow R implementation of the Takahashi
#' recursion method, unless a special build of the `fmesher` package is used.
#'
#' @param A A symmetric positive definite matrix
#'
#' @keywords internal
#' @export
#' @examples
#' A <- Matrix::Matrix(
#'   c(2, -1, 0, 0, -1, 2, -1, 0, 0, -1, 2, -1, 0, 0, -1, 2),
#'   4,
#'   4
#' )
#' (S <- fm_qinv(A))
#' (S2 <- solve(A))
#' c(sum((S - Matrix::t(S))^2), sum((S2 - Matrix::t(S2))^2))
#' sum((S - S2)[S != 0]^2)
fm_qinv <- function(A) {
  A_C <- fm_as_dgCMatrix(A)
  stopifnot(nrow(A_C) == ncol(A_C))
  if (!identical(A_C, Matrix::t(A_C))) {
    warning(
      "Asymmetric matrix A detected, ",
      "but only lower left triangle will be used."
    )
  }
  fmesher_qinv_R(A_C)
  # fmesher_qinv(A_C)
}


fmesher_qinv_R <- function(A) {
  C <- Matrix::Cholesky(A)
  LP <- Matrix::expand2(C)

  n <- nrow(A)

  S <- A
  for (i in rev(seq_len(n))) {
    if (i < n) {
      jj <- which(LP$L1[(i + 1L):n, i] != 0)
      if (length(jj)) {
        jj <- jj + i
      }
    } else {
      jj <- integer(0)
    }
    for (j in c(rev(jj), i)) {
      if (i == j) {
        S[i, i] <- 1 / LP$D[i, i]
        if (i < n) {
          S[i, i] <- S[i, i] - sum(LP$L1[jj, i] * S[jj, i])
        }
      } else {
        if (i < n) {
          S[i, j] <- -sum(LP$L1[jj, i] * S[jj, j])
          S[j, i] <- S[i, j]
        }
      }
    }
  }

  # A = P'LDL'P
  # S = A^-1 = P' (LDL')^-1 P
  S <- LP$P1. %*% S %*% LP$P1
  S
}
