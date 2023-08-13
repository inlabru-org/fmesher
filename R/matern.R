#' @title Mat\'ern process precision
#' @description
#' `r lifecycle::badge("experimental")`
#' Construct the (sparse) precision matrix for the basis weights for
#' Whittle-Mat\'ern SPDE models.
#'
#' @param x A mesh object
#' @param alpha The SPDE operator order. The smoothness index is `alpha - dim/2`.
#' @param rho The range parameter.
#' @param sigma The nominal standard deviation.
#'
#' @export
#' @examples
#' library(Matrix)
#' mesh <- fm_mesh_1d((0:10)^2, degree = 2)
#' Q <- fm_matern_precision(mesh, alpha = 2, rho = 15, sigma = 1)
#' x <- seq(0, 100, length.out = 301)
#' A <- fm_basis(mesh, x)
#' plot(x, as.vector(Matrix::diag(A %*% solve(Q, t(A)))), type = "l")
#'
#' plot(x, fm_evaluate(mesh, loc = x, field = fm_gmrf_sample(1, Q)[, 1]), type = "l")
fm_matern_precision <- function(x, alpha, rho, sigma) {
  mesh <- x
  d <- fm_manifold_dim(mesh)
  nu <- alpha - d / 2
  kappa <- sqrt(8 * nu) / rho
  scaling <- gamma(nu) / gamma(alpha) / (4 * pi)^(d / 2) / kappa^(2 * nu)

  fem <- fm_fem(mesh, order = ceiling(alpha))

  if (inherits(mesh, "fm_mesh_1d") && (mesh$degree <- 2)) {
    C <- fem$c1
  } else {
    C <- fem$c0
  }
  if (alpha == 2) {
    Q <- (C * kappa^4 + 2 * kappa^2 * fem$g1 + fem$g2) / sigma^2 * scaling
  } else if (alpha == 1) {
    Q <- (C * kappa^2 + fem$g1) / sigma^2 * scaling
  } else {
    stop("This version only supports alpha = 1 and 2.")
  }

  Q
}


#' @describeIn fm_matern_precision
#' Generate `n` samples based on a sparse precision matrix `Q`
#' @param n The number of samples to generate
#' @param Q A precision matrix
#' @export
fm_gmrf_sample <- function(n, Q) {
  # Find P and L such that P Q P' = L L',
  # i.e. Q = P' L L' P and Q^-1 = P' solve(L L') P
  fact <- Matrix::expand(Matrix::Cholesky(Q, perm = TRUE))
  # L' P x = w gives L' P S_x P' L = I, S_x = P' (L L')^-1 P
  # so we need to solve L' x0 = w and then compute x = P' x0
  x0 <- Matrix::solve(
    fact$L,
    Matrix::Matrix(stats::rnorm(n * nrow(Q)), nrow(Q), n)
  )
  as.matrix(Matrix::solve(fact$P * 1, x0))
}
