#' @importFrom Matrix t solve
#' @title PC Prior Calculation
#' @description Calculates  the PC prior based on given hyperparameters and vectors.
#'
#' @param lambda A hyperparameter controlling the size of kappa.
#' @param lambda1 A hyperparameter controlling the size of |v|.
#' @param log_kappa Logarithm of inverse correlation range.
#' @param v Vector that controls anisotropy.
#'
#' @return The calculated PC prior density .
#' @export
#' @examples
#' lambda <- 1
#' lambda1 <- 1
#' log_kappa <- 0.5
#' v <- c(1, 2)
#' pc_prior(lambda = lambda, lambda1 = lambda1, log_kappa = log_kappa, v = v)
pc_prior <- function(lambda, lambda1, log_kappa, v) {
  kappa <- exp(log_kappa)
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

#' @title Log PC Prior Calculation for anisotropy parameters.
#' @description Calculates  the log of the PC prior for (kappa, v) based on given hyperparameters and vectors.
#'
#' @param lambda A hyperparameter controlling the size of kappa.
#' @param lambda1 A hyperparameter controlling the size of |v|.
#' @param log_kappa Logarithm of inverse correlation range.
#' @param v Vector that controls anisotropy.
#'
#' @return The calculated log PC prior value.
#' @export
#' @examples
#' lambda <- 1
#' lambda1 <- 1
#' log_kappa <- 0.5
#' v <- c(1, 2)
#' log_pc_prior_aniso(lambda = lambda, lambda1 = lambda1, log_kappa = log_kappa, v = v)
log_pc_prior_aniso <- function(lambda, lambda1, log_kappa, v) {
  # Calculates  the log PC prior
  log_pc_prior_aniso_value <- log(pc_prior(lambda = lambda, lambda1 =  lambda1 , log_kappa = log_kappa, v = v))

  return(log_pc_prior_aniso_value)
}

#' @title Log PC Prior calcualtion for variance of noise
#' @description Calculates  the log of the PC prior for sigma_epsilon based on given hyperparameters and vectors.
#'
#' @param lambda_epsilon A hyperparameter controlling the size of epsilon.
#' @param log_sigma_epsilon The logarithm of the variance of the additive noise.
#'
#' @return The calculated log PC prior value.
#' @export
#' @examples
#' lambda_epsilon <- 1
#' log_pc_prior_noise_variance(lambda_epsilon = lambda_epsilon, log_sigma_epsilon = log_sigma_epsilon)
log_pc_prior_noise_variance <- function(lambda_epsilon, log_sigma_epsilon) {
  sigma_epsilon = exp(log_sigma_epsilon)
  # Calculates  the logarithm of exponential density.
  return(log(lambda_epsilon) - lambda_epsilon * sigma_epsilon)
}

#' @title Sparse Matrix Determinant using Cholesky Factorization
#' @description Calculates  the determinant of a sparse matrix using Cholesky factorization.
#'
#' @param Q A sparse matrix.
#'
#' @return The determinant of the sparse matrix.
#' @export
#' @examples
#' library(Matrix)
#' Q <- Matrix(c(4, -1, -1, 4), nrow = 2, sparse = TRUE)
#' sparse_log_determinant_chol(Q)
sparse_log_determinant_chol <- function(Q) {
  # Perform Cholesky factorization
  chol_fact <- Matrix::Cholesky(Q, perm = TRUE, LDL = FALSE)

  # Extract square diagonal elements of the Cholesky factor L
  diag_L <- Matrix::diag(chol_fact)

  # Calculates  determinant using Cholesky factor
  det_val <- sum(log(diag_L))

  return(det_val)
}



#' @title Calculates  the norm ||x||_Q:= <Qx,x>
#'
#' @description
#' Calculates  the term <Qx, x> appearing in Gaussian density
#'
#' @param Q A sparse matrix representing the precision matrix
#' @param x A numeric vector representing the point x
#'
#' @return The calculated norm term
#' @export
norm_Q <- function(Q, x) {
  norm <- sum(x * (Q %*% x))
  return(norm)
}

#' @title Calculates  the log Gaussian density
#'
#' @description
#' Calculates  the log Gaussian density using the precision matrix Q and mean mu
#'
#' @param x A numeric vector representing the point x
#' @param mu A numeric vector representing the mean mu
#' @param Q A sparse matrix representing the precision matrix
#'
#' @return The calculated log Gaussian density
#' @export
logGdensity <- function(x, mu, Q) {
  # Calculates  determinant using Cholesky factorization
  log_det_val <- sparse_log_determinant_chol(Q)

  # Calculates  the norm term
  norm_term <- norm_Q(Q, x - mu)

  # Calculates  the log Gaussian density
  logGdty <- 0.5 * (log_det_val - norm_term) - nrow(Q) / 2 * log(2 * pi)

  return(logGdty)
}

#' @title Calculates  the log-posterior density for parameters (log_kappa,v,log_sigma_epsilon) with PC prior.
#'
#' @description
#' Calculates  the log-posterior density based on the prior density and the likelihood. Only stationary parameters are accepted.
#'
#' @param mesh The mesh
#' @param log_kappa Logarithm of inverse correlation range
#' @param v 2D vector that controls anisotropy
#' @param log_sigma_epsilon Variance of noise
#' @param lambda A hyperparameter controlling the size of kappa.
#' @param lambda1 A hyperparameter controlling the size of |v|.
#' @param lambda_epsilon A hyperparameter controlling the size of sigma_epsilon.
#' @param y A vector with length equal to the number of basis elements n representing the observed data.
#' @param A Matrix of size nxn representing the transformation A
#' @param Q_epsilon A sparse matrix of size nxn representing the noise precision matrix
#' @param m_u A vector with length n representing the prior mean m_u. If a number is given is transformed into (m_u, m_u,..., m_u)
#'
#' @return The calculated log-posterior
#' @export


log_posterior <- function(mesh, log_kappa, log_sigma_epsilon, v, lambda, lambda1, lambda_epsilon, y, A, m_u) {
  kappa <-exp(log_kappa)
  sigma_epsilon <- exp(log_sigma_epsilon)
  # Calculates  log-prior
  log_pc_aniso_value <- log_pc_prior_aniso(lambda = lambda, lambda1 = lambda1, log_kappa = log_kappa, v = v)
  log_pc_noise_value <- log_pc_prior_noise_variance(lambda_epsilon =lambda_epsilon, log_sigma_epsilon = log_sigma_epsilon)
  log_pc_value <- log_pc_aniso_value + log_pc_noise_value

  # Calculates  anisotropy
  n <- nrow(mesh$loc)
  kappa_values <- rep(kappa, n)
  vec_values <- matrix(v, n, 2, byrow = TRUE)
  aniso <- list(kappa = kappa_values, vec = vec_values)

  # Calculates  log-density of the distribution of u knowing (kappa, v)
  Q_u <- fm_aniso_precision(mesh, aniso)
  if (length(m_u) == 1) {
    m_u <- rep(m_u, n)
  }
  u <- m_u
  logGdty_prior <- logGdensity(x = u, mu = m_u, Q = Q_u)

  # Calculates Q_epsilon,  Q_{u|y,theta} and m_{u|y,theta}
  Q_epsilon <-Matrix::Diagonal(n, sigma_epsilon^2)
  Q_uy_theta <- Q_u + t(A) %*% Q_epsilon %*% A
  m_uy_theta <- solve(Q_uy_theta, Q_u %*% m_u + t(A) %*% Q_epsilon %*% y)

  # Calculates  log-density of the posterior distribution of u given y and theta
  logGdty_posterior <- logGdensity(x = u, mu = m_uy_theta, Q = Q_uy_theta)

  # Calculates  log-density of the observation of y|u, theta
  logGdty_observation <- logGdensity(x = y, mu = A %*% u, Q = Q_epsilon)

  # Calculates  log-posterior
  log_posterior_val <- log_pc_value + logGdty_prior + logGdty_observation - logGdty_posterior

  return(log_posterior_val)
}

#' @title Calculates  the MAP estimate for linear noisy observation of field.
#'
#' @description
#' Calculated by maximizing log posterior using optim. Only stationary parameters accepted.
#'
#' @param mesh The mesh
#' @param lambda A hyperparameter controlling the size of kappa.
#' @param lambda1 A hyperparameter controlling the size of |v|.
#' @param lambda_epsilon A hyperparameter controlling the size of sigma_epsilon.
#' @param y A vector with length equal to the number of basis elements n representing the observed data.
#' @param A Matrix of size nxn representing the transformation A
#' @param m_u A vector with length n representing the prior mean m_u
#'
#' @return The calculated log-posterior
#' @export


MAP <- function(mesh, lambda, lambda1, lambda_epsilon, y, A, m_u) {
  # Writes the log-posterior as a function of (log_kappa, v, sigma)
  log_post <- function(theta) {
    log_kappa <- theta[1]
    v <- theta[2:3]
    log_sigma_epsilon <- theta[4]
    return(log_posterior(mesh = mesh, log_kappa = log_kappa, v = v, log_sigma_epsilon = log_sigma_epsilon, lambda = lambda, lambda1 = lambda1, lambda_epsilon = lambda_epsilon, y = y, A = A, m_u = m_u))
  }
  aniso_0 <- c(log(0.5), c(1,2), 1)
  # To do: calculate the gradient of log posterior
  # gradient= grad_log_posterior(mesh, kappa, v, lambda, lambda1, y, A, Q_epsilon, m_u)
  return(optim(par = aniso_0, fn = log_post, control= list(fnscale = -1, maxit=300)))
}

#' @title Calculates  the gradient of the log posterior of a linear observation y = A u + noise
#'
#' @description
#' Calculates  the gradient of the log posterior of a linear observation y = A u + noise
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
