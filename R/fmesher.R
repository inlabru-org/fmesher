#' Sparse partial inverse
#'
#' Compute sparse partial matrix inverse. As of `0.2.0.9010`, an R
#' implementation of the Takahashi recursion method, unless a special build of
#' the `fmesher` package is used.
#'
#' @param A A sparse symmetric positive definite matrix
#' @returns A sparse symmetric matrix, with the elements of the inverse of `A`
#' for the non-zero pattern of `A` plus potential Cholesky in-fill locations.
#'
#' @export
#' @examples
#' A <- Matrix::Matrix(
#'   c(2, -1, 0, 0, -1, 2, -1, 0, 0, -1, 2, -1, 0, 0, -1, 2),
#'   4,
#'   4
#' )
#' # Partial inverse:
#' (S <- fm_qinv(A))
#' # Full inverse (not guaranteed to be symmetric):
#' (S2 <- solve(A))
#' # Matrix symmetry:
#' c(sum((S - Matrix::t(S))^2), sum((S2 - Matrix::t(S2))^2))
#' # Accuracy (not that S2 is non-symmetric, and S may be more accurate):
#' sum((S - S2)[S != 0]^2)
fm_qinv <- function(A) {
  A_C <- fm_as_dgCMatrix(A)
  stopifnot(nrow(A_C) == ncol(A_C))
  # if (sum(abs(A_C - Matrix::t(A_C))) > 0.0) {
  #   warning("Asymmetric matrix A detected, ",
  #           "but only lower left triangle will be used.")
  # }
  fmesher_qinv_R(A_C)
  # fmesher_qinv(A_C)
}


fmesher_qinv_R <- function(A) {
  C <- Matrix::Cholesky(A)
  LP <- Matrix::expand2(C)

  n <- nrow(A)

  S <- A
  for (i in rev(seq_len(n))) {
    if (i == n) {
      S[i, i] <- 1 / LP$D[i, i]
    } else {
      jj <- sort(unique(
        c(
          which(LP$L1[(i + 1L):n, i] != 0),
          which(S[i, (i + 1L):n] != 0)
        )
      ))
      if (length(jj) > 0) {
        jj <- jj + i
        Lvals <- LP$L1[jj, i]
        result <- -as.vector(Lvals %*% S[jj, jj])
        S[i, jj] <- result
        S[jj, i] <- result
        S[i, i] <- 1 / LP$D[i, i] - as.vector(Lvals %*% S[jj, i])
      } else {
        S[i, i] <- 1 / LP$D[i, i]
      }
    }
  }

  # A = P'LDL'P
  # S = A^-1 = P' (LDL')^-1 P
  S <- LP$P1. %*% S %*% LP$P1
  S
}
