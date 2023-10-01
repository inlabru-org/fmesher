#' @title PC Prior Calculation
#' @description Calculate the PC prior based on given hyperparameters and vectors.
#'
#' @param lambda A hyperparameter controlling the size of kappa.
#' @param lambda1 A hyperparameter controlling the size of |v|.
#' @param kappa Inverse correlation range
#' @param v Vector that controls anisotropy
#'
#' @return The calculated PC prior value.
#' @export
#' @examples
#' lambda <- 1
#' lambda1 <- 1
#' kappa <- 0.5
#' v <- c(1, 2)
#' pc_prior(lambda, lambda1, kappa, v)
pc_prior <- function(lambda, lambda1, kappa, v) {
  f <- function(r) {
    sqrt((1 / (48 * pi)) * (3 * cosh(2 * r) + 1))
  }

  f_prime <- function(r) {
    sinh(2 * r) / (4 * sqrt(pi * cosh(2 * r) + pi / 3))
  }

  v_norm <- sqrt(sum(v^2))
  f_val <- f(v_norm)
  f_prime_val <- f_prime(v_norm)

  term1 <- (lambda * lambda1 * abs(f_prime_val) * f_val) / (2 * pi * v_norm)
  term2 <- exp(-lambda1 * (f_val - f(0)) - lambda * f_val * kappa)

  return(term1 * term2)
}

# Example usage
lambda <- 1
lambda1 <- 1
kappa <- 0.5
v <- c(1, 2)

result <- pc_prior(lambda, lambda1, kappa, v)
print(result)

#' @title Log PC Prior Calculation
#' @description Calculate the log of the PC prior based on given hyperparameters and vectors.
#'
#' @param lambda A hyperparameter controlling the scale.
#' @param lambda1 A hyperparameter controlling the shape.
#' @param kappa A parameter related to the base model.
#' @param v A 2D vector.
#'
#' @return The calculated log PC prior value.
#' @export
#' @examples
#' lambda <- 1
#' lambda1 <- 1
#' kappa <- 0.5
#' v <- c(1, 2)
#' log_pc_prior(lambda, lambda1, kappa, v)
log_pc_prior <- function(lambda, lambda1, kappa, v) {

  # Calculate the log PC prior
  log_pc_prior_value <- log(pc_prior(lambda, lambda1, kappa, v))

  return(log_pc_prior_value)
}

#' @title Sparse Matrix Determinant using Cholesky Factorization
#' @description Calculate the determinant of a sparse matrix using Cholesky factorization.
#'
#' @param Q A sparse matrix.
#'
#' @return The determinant of the sparse matrix.
#' @export
#' @examples
#' library(Matrix)
#' Q <- Matrix(c(4, -1, -1, 4), nrow = 2, sparse = TRUE)
#' sparse_determinant_chol(Q)
sparse_determinant_chol <- function(Q) {
  # Perform Cholesky factorization
  chol_fact <- Matrix::Cholesky(Q, perm = TRUE, LDL = FALSE)

  # Extract diagonal elements of the Cholesky factor L
  diag_L <- Matrix::diag(chol_fact)

  # Calculate determinant using Cholesky factor
  det_val <- prod(diag_L)

  return(det_val)
}



#' @title Calculate the norm ||x||_Q:= <Qx,x>
#'
#' @description
#' Calculate the term <Qx, x> appearing in Gaussian density
#'
#' @param Q A sparse matrix representing the precision matrix
#' @param x A numeric vector representing the point x
#'
#' @return The calculated norm term
#' @export
norm_Q <- function(Q, x) {
  norm <- t(x) %*% Q %*% x
  return(norm)
}

#' @title Calculate the log Gaussian density
#'
#' @description
#' Calculate the log Gaussian density using the precision matrix Q and mean mu
#'
#' @param Q A sparse matrix representing the precision matrix
#' @param x A numeric vector representing the point x
#' @param mu A numeric vector representing the mean mu
#'
#' @return The calculated log Gaussian density
#' @export
logGdensity <- function(Q, x, mu) {
  # Calculate determinant using Cholesky factorization
  det_val <- sparse_determinant_chol(Q)

  # Calculate the norm term
  norm_term <- norm_Q(Q, x - mu)

  # Calculate the log Gaussian density
  logGdty <- -0.5 * (log(det_val) + norm_term)

  return(logGdty)
}

