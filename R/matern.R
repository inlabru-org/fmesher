#' @title SPDE, GMRF, and Matérn process methods
#' @description
#' `r lifecycle::badge("experimental")`
#' Methods for SPDEs and GMRFs.
#' @name fm_gmrf
NULL

#' @describeIn fm_gmrf
#' Construct the (sparse) precision matrix for the basis weights for
#' Whittle-Matérn SPDE models.  The boundary behaviour is determined by the
#' provided mesh function space.
#'
#' @param x A mesh object, e.g. from `fm_mesh_1d()` or `fm_mesh_2d()`.
#' @param alpha The SPDE operator order. The resulting smoothness index
#' is `nu = alpha - dim / 2`.
#' @param rho The Matérn range parameter
#' (scale parameter `kappa = sqrt(8 * nu) / rho`)
#' @param sigma The nominal Matérn std.dev. parameter
#'
#' @export
#' @examples
#' library(Matrix)
#' mesh <- fm_mesh_1d(-20:120, degree = 2)
#' Q <- fm_matern_precision(mesh, alpha = 2, rho = 15, sigma = 1)
#' x <- seq(0, 100, length.out = 601)
#' A <- fm_basis(mesh, x)
#' plot(x,
#'   as.vector(Matrix::diag(fm_covariance(Q, A))),
#'   type = "l",
#'   ylab = "marginal variances"
#' )
#'
#' plot(x,
#'   fm_evaluate(mesh, loc = x, field = fm_sample(1, Q)[, 1]),
#'   type = "l",
#'   ylab = "process sample"
#' )
#'
fm_matern_precision <- function(x, alpha, rho, sigma) {
  mesh <- x
  d <- fm_manifold_dim(mesh)
  nu <- alpha - d / 2
  kappa <- sqrt(8 * nu) / rho
  scaling <- gamma(nu) / gamma(alpha) / (4 * pi)^(d / 2) / kappa^(2 * nu)

  fem <- fm_fem(mesh, order = ceiling(alpha))

  if (inherits(mesh, "fm_mesh_1d") && (mesh$degree == 2)) {
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

#' @describeIn fm_gmrf
#' Construct the (sparse) precision matrix for the basis weights of anisotropic
#' SPDE models.
#'
#' @param x A mesh object, e.g. from `fm_mesh_2d()`.
#' @param aniso List `[kappa,vec]` where `kappa` controls the (inverse) correlation range
#' and (the half angle version of) `vec` controls the main directions of the anisotropy
#' @param log_sima The logarithmm of the marginal variance sigma (is a constant)
#' @export
#' @examples
#' .....
fm_aniso_precision <- function(x, aniso, log_sigma = 0) {
  sigma <- exp(log_sigma)
  scaling <- 1 / (4 * pi * sigma^2) #Calculates scaling so that Q_fem * scaling has variance sigma
  fem <- fm_fem_aniso(x, aniso)
  Q <- (fem$c0 + 2 * fem$g1 + fem$g2) * scaling
  Q
}



#' @describeIn fm_gmrf
#' Simulate a Matérn field given a mesh and
#' covariance function parameters, and optionally evaluate at given locations.
#'
#' @param loc locations to evaluate the random field, compatible with
#' `fm_evaluate(x, loc = loc, field = ...)`
#'
#' @return `fm_matern_sample()` returns a matrix, where each column is a sampled
#' field. If `loc` is `NULL`, the `fm_dof(mesh)` basis weights are given.
#' Otherwise, the evaluated field at the `nrow(loc)` locations `loc` are given.
#' @export

fm_matern_sample <- function(x, alpha = 2, rho, sigma, n = 1, loc = NULL) {
  Q <- fm_matern_precision(x, alpha = alpha, rho = rho, sigma = sigma)
  x <- fm_sample(n = n, Q = Q)
  if (!is.null(loc)) {
    x <- fm_evaluate(x, loc = loc, field = x)
  }
  x
}

#' @describeIn fm_gmrf
#' Simulate an anisotropic field given a mesh and
#' anisotropic parameters, and optionally evaluate at given locations.
#'
#' @param loc locations to evaluate the random field, compatible with
#' `fm_evaluate(x, loc = loc, field = ...)`
#'
#' @return `fm_aniso_sample()` returns a matrix, where each column is a sampled
#' field. If `loc` is `NULL`, the `fm_dof(mesh)` basis weights are given.
#' Otherwise, the evaluated field at the `nrow(loc)` locations `loc` are given.
#' @export

fm_aniso_sample <- function(x, aniso, n = 1, loc = NULL) {
  Q <- fm_aniso_precision(x, aniso)
  x <- fm_sample(n = n, Q = Q)
  if (!is.null(loc)) {
    x <- fm_evaluate(x, loc = loc, field = x)
  }
  x
}

#' @describeIn fm_gmrf
#' Simulates the basis weights u_i in u(x) = sum(u_j phi_j(x))
#'
#' @param x A 2d mesh object.
#' @param aniso List `[kappa,vec]` where `kappa` controls the (inverse) correlation range
#' and (the half angle version of) `vec` controls the main directions of the anisotropy
#' @param sigma The std.dev. parameter of anisotropic u
#'
#' @return `fm_aniso_sample()` returns a vector whose j_th component is a sample of the
#' weight u_j of the jth basis vector.
#' @export

fm_aniso_basis_weights_sample <- function(x, aniso, n = 1, log_sigma = 0) {
  Q <- fm_aniso_precision(x, aniso, log_sigma) # Calculate the precision
  L_solve <- function(fact, b) {
    Matrix::solve(fact, Matrix::solve(fact, b, system = "P"), system = "L")
  }
  Lt_solve <- function(fact, b) {
    Matrix::solve(fact, Matrix::solve(fact, b, system = "Lt"), system = "Pt")
  }
  # Find P and L such that Q = P' L L' P
  fact <- Matrix::Cholesky(Q, perm = TRUE, LDL = FALSE)
  # To find the value at the weights we need to solve PLu = w
  u <- Lt_solve(
    fact,
    Matrix::Matrix(stats::rnorm(n * nrow(Q)), nrow(Q), n)
  )
  return(u)
}



#' @describeIn fm_gmrf Compute the covariance between "A1 x" and "A2 x", when
#' x is a basis vector with precision matrix `Q`.
#' @param A1,A2 Matrices, typically obtained from [fm_basis()] and/or [fm_block()].
#' @param partial `r lifecycle::badge("experimental")` If `TRUE`, compute the
#' partial inverse of `Q`, i.e. the elements of the inverse corresponding to
#' the non-zero pattern of `Q`. (Note: This can be done efficiently with
#' the Takahashi recursion method, but to avoid an RcppEigen dependency this
#' is currently disabled, and a slower method is used until the efficient method
#' is reimplemented.)
#' @export
fm_covariance <- function(Q, A1 = NULL, A2 = NULL, partial = FALSE) {
  if (is.null(A2)) {
    A2 <- A1
  }
  if (partial) {
    if (!is.null(A1)) {
      warning("'partial=TRUE', but `A1` is not NULL; ignoring `A1` and `A2`.")
    }
    fact <- Matrix::Cholesky(as(Q, "sparseMatrix"), perm = TRUE)
    Q <- fm_as_dgTMatrix(Q)
    Q_idx <- data.frame(i = Q@i, j = Q@j)
    block_size <- 50
    n_blocks <- (nrow(Q) %/% block_size) + (nrow(Q) %% block_size > 0L)
    output <- list(i = integer(0), j = integer(0), x = numeric(0))
    for (k in seq_len(n_blocks)) {
      i_offset <- (k - 1L) * block_size
      i_len <- min(i_offset + block_size, nrow(Q)) - i_offset
      e <- Matrix::sparseMatrix(
        i = seq_len(i_len) + i_offset,
        j = seq_len(i_len),
        x = rep(1, i_len),
        dims = c(nrow(Q), i_len)
      )
      Q_idx_block <-
        Q_idx[(Q_idx$j + 1L > i_offset) &
          (Q_idx$j + 1L <= i_offset + block_size), , drop = FALSE]
      result <- fm_as_dgTMatrix(Matrix::solve(fact, e))
      ok <- (result@i %in% Q_idx_block$i)
      ok[ok] <- ((result@j[ok] + i_offset) %in% Q_idx_block$j)
      idx <-
        base::merge(
          Q_idx_block,
          data.frame(
            i = result@i[ok],
            j = result@j[ok] + i_offset,
            x = result@x[ok]
          ),
          sort = FALSE
        )
      output$i <- c(output$i, idx$i + 1L)
      output$j <- c(output$j, idx$j + 1L)
      output$x <- c(output$x, idx$x)
    }
    return(Matrix::sparseMatrix(
      i = output$i,
      j = output$j,
      x = output$x,
      dims = c(nrow(Q), nrow(Q))
    ))
  }

  if (is.null(A1) || is.null(A2)) {
    return(Matrix::solve(Q))
  }
  A1 %*% Matrix::solve(Q, Matrix::t(A2))
}


#' @describeIn fm_gmrf
#' Generate `n` samples based on a sparse precision matrix `Q`
#' @param n The number of samples to generate
#' @param Q A precision matrix
#' @param mu Optional mean vector
#' @param constr Optional list of constraint information, with elements
#' `A` and `e`. Should only be used for a small number of exact constraints.
#' @export
fm_sample <- function(n, Q, mu = 0, constr = NULL) {
  L_solve <- function(fact, b) {
    Matrix::solve(fact, Matrix::solve(fact, b, system = "P"), system = "L")
  }
  Lt_solve <- function(fact, b) {
    Matrix::solve(fact, Matrix::solve(fact, b, system = "Lt"), system = "Pt")
  }

  # TODO: Work out how to use LDLt factorisations; How to access D^(0.5) ?
  #
  # Find P and L such that P Q P' = L L',
  # i.e. Q = P' L L' P and Q^-1 = P' solve(L L') P = P' L^-T L^-1 P
  fact <- Matrix::Cholesky(Q, perm = TRUE, LDL = FALSE)
  # Not currently needed: fact_exp <- Matrix::expand(fact)
  # L' P x = w gives L' P S_x P' L = I, S_x = P' (L L')^-1 P
  # so we need to solve L' x0 = w and then compute x = P' x0
  x <- Lt_solve(
    fact,
    Matrix::Matrix(stats::rnorm(n * nrow(Q)), nrow(Q), n)
  )
  result <- mu + x
  if (!is.null(constr) && !is.null(constr[["A"]])) {
    if (is.null(constr[["e"]])) {
      constr$e <- matrix(0, nrow(constr$A), n)
    } else {
      constr$e <- as.matrix(constr$e)
      if (ncol(constr$e) == 1) {
        constr$e <- matrix(as.vector(constr$e), nrow(constr$A), n)
      }
    }
    # See gmrf.pdf section 4.2
    A_tilde_T <- L_solve(fact, Matrix::t(constr$A))
    result <-
      result - Lt_solve(
        fact,
        qr.solve(
          Matrix::t(A_tilde_T),
          constr$A %*% result - constr$e
        )
      )
    # Keeping this version for when soft constraints are added, as it is closer
    # to what's needed for that:
    # result <- result - Lt_solve(
    #   fact,
    #   A_tilde_T %*%
    #     Matrix::solve(
    #       Matrix::t(A_tilde_T) %*% A_tilde_T,
    #       constr$A %*% result - constr$e
    #     )
    # )
  }
  return(as.matrix(result))
}
