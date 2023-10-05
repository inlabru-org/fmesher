#' @importFrom Matrix t solve
#' @title PC Prior Calculation
#' @description Calculate the PC prior based on given hyperparameters and vectors.
#'
#' @param lambda A hyperparameter controlling the size of kappa.
#' @param lambda1 A hyperparameter controlling the size of |v|.
#' @param kappa Inverse correlation range
#' @param v Vector that controls anisotropy
#'
#' @return The calculated PC prior density .
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
log_pc_prior <- function(kappa, v, lambda, lambda1, ) {
  # Calculate the log PC prior
  log_pc_prior_value <- log(pc_prior(lambda = lambda, lambda1 =  lambda1 , kappa = kappa, v = v))

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
sparse_log_determinant_chol <- function(Q) {
  # Perform Cholesky factorization
  chol_fact <- Matrix::Cholesky(Q, perm = TRUE, LDL = FALSE)

  # Extract square diagonal elements of the Cholesky factor L
  diag_L <- Matrix::diag(chol_fact)

  # Calculate determinant using Cholesky factor
  det_val <- sum(log(diag_L))

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
  # x <- as.vector(x)
  norm <- sum(x * (Q %*% x))
  return(norm)
}

#' @title Calculate the log Gaussian density
#'
#' @description
#' Calculate the log Gaussian density using the precision matrix Q and mean mu
#'
#' @param x A numeric vector representing the point x
#' @param mu A numeric vector representing the mean mu
#' @param Q A sparse matrix representing the precision matrix
#'
#' @return The calculated log Gaussian density
#' @export
logGdensity <- function(x, mu, Q) {
  # Calculate determinant using Cholesky factorization
  log_det_val <- sparse_log_determinant_chol(Q)

  # Calculate the norm term
  norm_term <- norm_Q(Q, x - mu)

  # Calculate the log Gaussian density
  logGdty <- 0.5 * (log_det_val - norm_term) - nrow(Q) / 2 * log(2 * pi)

  return(logGdty)
}

#' @title Calculate the log-posterior for anisotropic parameters (kappa,v) with PC prior.
#'
#' @description
#' Calculate the log-posterior based on the prior and likelihood. Only stationary parameters accepted.
#'
#' @param mesh The mesh
#' @param kappa Inverse correlation range
#' @param v Vector that controls anisotropy
#' @param lambda A hyperparameter controlling the size of kappa.
#' @param lambda1 A hyperparameter controlling the size of |v|.
#' @param y A vector with length equal to the number of basis elements n representing the observed data.
#' @param A Matrix of size nxn representing the transformation A
#' @param Q_epsilon A sparse matrix of size nxn representing the noise precision matrix
#' @param m_u A vector with length n representing the prior mean m_u
#'
#' @return The calculated log-posterior
#' @export


log_posterior <- function(mesh, kappa, v, lambda, lambda1, y, A, Q_epsilon, m_u) {
  # Calculate log-prior
  log_pc_value <- log_pc_prior(kappa, v, lambda, lambda1)

  # Calculate anisotropy
  n <- nrow(mesh$loc)
  kappa_values <- rep(kappa, n)
  vec_values <- matrix(v, n, 2, byrow = TRUE)
  aniso <- list(kappa = kappa_values, vec = vec_values)

  # Calculate log-density of the distribution of u knowing (kappa, v)
  Q_u <- fm_aniso_precision(mesh, aniso)
  if (length(m_u) == 1) {
    m_u <- rep(m_u, n)
  }
  u <- m_u
  logGdty_prior <- logGdensity(u, m_u, Q_u)

  # Calculate Q_{u|y,theta} and m_{u|y,theta}
  Q_uy_theta <- Q_u + t(A) %*% Q_epsilon %*% A
  m_uy_theta <- solve(Q_uy_theta, Q_u %*% m_u + t(A) %*% Q_epsilon %*% y)

  # Calculate log-density of the posterior distribution of u given y and theta
  logGdty_posterior <- logGdensity(u, m_uy_theta, Q_uy_theta)

  # Calculate log-density of the observation of y|u, theta
  logGdty_observation <- logGdensity(y, A %*% u, Q_epsilon)

  # Calculate log-posterior
  log_posterior_val <- log_pc_value + logGdty_prior + logGdty_observation - logGdty_posterior

  return(log_posterior_val)
}

#' @title Calculate the MAP estimate for linear noisy observation of field.
#'
#' @description
#' Calculated by maximizing log posterior using optim. Only stationary parameters accepted.
#'
#' @param mesh The mesh
#' @param lambda A hyperparameter controlling the size of kappa.
#' @param lambda1 A hyperparameter controlling the size of |v|.
#' @param y A vector with length equal to the number of basis elements n representing the observed data.
#' @param A Matrix of size nxn representing the transformation A
#' @param Q_epsilon A sparse matrix of size nxn representing the noise precision matrix
#' @param m_u A vector with length n representing the prior mean m_u
#'
#' @return The calculated log-posterior
#' @export


MAP <- function(mesh, lambda, lambda1, y, A, Q_epsilon, m_u) {
  # Writes the log-posterior as a function of (kappa,v)
  log_post <- function(theta) {
    log_kappa <- theta[1]
    v <- theta[2:3]
    # epsilon =
    return(log_posterior(mesh = mesh, kappa = exp(log_kappa), v = v, lambda = lambda, lambda1 = lambda1, y = y, A = A, Q_epsilon = Q_epsilon, m_u = m_u))
  }
  aniso_0 <- c(0, 0, 0)
  # To do: calculate the gradient of log posterior
  # gradient= grad_log_posterior(mesh, kappa, v, lambda, lambda1, y, A, Q_epsilon, m_u)
  return(aniso_0, optim(aniso_0, log_post))
}

#' @title Calculate the gradient of the log posterior of a linear observation y = A u + noise
#'
#' @description
#' Calculate the gradient of the log posterior of a linear observation y = A u + noise
#'
#' @param mesh The mesh
#' @param kappa Inverse correlation range
#' @param v Vector that controls anisotropy
#' @param lambda A hyperparameter controlling the size of kappa.
#' @param lambda1 A hyperparameter controlling the size of |v|.
#' @param y A vector with length equal to the number of basis elements n representing the observed data.
#' @param A Matrix of size nxn representing the transformation A
#' @param Q_epsilon A sparse matrix of size nxn representing the noise precision matrix
#' @param m_u A vector with length n representing the prior mean m_u
#'
#' @return The calculated log-posterior
#' @export


grad_log_posterior <- function(mesh, kappa, v, lambda, lambda1, y, A, Q_epsilon, m_u) {
  # Calculation
  grad_log_prior <- 0
  # d|Q|= |Q| Tr(Q^{-1}dQ)
  # Need to differentiate Q = C_kappa + 2 G_v + G_v C_kp^{-1} G_v
  # Equivalent to differentiating C, and G component by component.
  # Should we add the calculation of gradient over each triangle to the calcCaniso, calcGaniso, calcQaniso functions
  # dC^{-1} = - dC * M * dC   (don't forget lumped mass)
  # Should I calDo we need to differentitate manifold case?
  grad_logGdty_prior <- 0
  grad_logGdty_posterior <- 0
  return(grad_log_prior + grad_logGdty_prior + grad_logGdty_posterior)
}
