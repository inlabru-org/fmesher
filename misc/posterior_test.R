library(devtools)
library(ggplot2)
library(Matrix)
library(sp)
library(INLA)

#Hyperparameters for PC priors
lambda <- 1; lambda1 <- 1; lambda_epsirlon <- 1; lambda_u <- 1

#Anisotropy parameters
kappa <- 3; log_kappa <- log(kappa);
v <- c(1,2)

#Correlation range calculation
rho <- sqrt(8)/kappa/sqrt(exp(sqrt(v[1]^2 + v[2]^2)))

#Noise parameters
sigma_u <- 0.01; log_sigma_u <- log(sigma_u)
sigma_epsilon <- 0.01 ; log_sigma_epsilon <- log(sigma_epsilon)


#Testing PC priors
log_pc_aniso <- log_pc_prior_aniso(lambda, lambda1, kappa, v)
log_pc_noise <- log_pc_prior_noise_variance(lambda_epsilon = lambda_epsilon, log_sigma_epsilon = log_sigma_epsilon)
log_pc_sigma_u <- log_pc_prior_noise_variance(lambda_epsilon = lambda_u, log_sigma_epsilon = log_sigma_u)
log_pc_value <- log_pc_aniso + log_pc_noise + log_pc_sigma_u
#Mesh definition
library(sf)
boundary_sf = st_sfc(st_polygon(list(rbind(c(0, 0), c(10, 0), c(10, 10), c(0, 10),c(0,0)))))
boundary = fm_as_segm(boundary_sf)
mesh <- fm_mesh_2d_inla(  boundary = boundary,  max.edge = c(0.5, 0.5))
nodes <- mesh$loc
n <- mesh$n
plot(mesh)

#Defining anisotropy
kappa_values <- rep(kappa, n)
vec_values <- matrix(v, n, 2, byrow = TRUE)
aniso <- list(kappa = kappa_values, vec = vec_values)

#Calculating precision and testing log determinant
Q <-fm_aniso_precision(mesh, aniso, log_sigma = log_sigma_u)
log_det <- sparse_log_determinant_chol(Q)
#Checking value
log_det_2 <- determinant(Q,logarithm = TRUE)
log_det_2_v <- as.numeric(log_det_2$modulus)
log_det - log_det_2_v

#Calculating average marginal variance, should give sigma_u
simulated_variance <- mean(diag(INLA::inla.qinv(Q)))

#Calculating |x|_Q:= x^T * Q * x
x <- rep(2,n)
x_norm <- norm_Q(Q,x)

#Calculating log Gaussian density
m_u = as.vector(rep(0,n))
logdtyQ <- logGdensity(x,m_u,Q) #This should be small if x != m_u and large if x = m_u

#Sampling from noisy data
x = fm_aniso_basis_weights_sample(x = mesh, aniso = aniso, log_sigma = log_sigma_u)
A = Matrix::Diagonal(n, 1)
y = A %*% x + exp(log_sigma_epsilon) * stats::rnorm(nrow(Q))

#Calculate log posterior and map
log_posterior_true <- log_posterior(mesh = mesh,
                                    log_kappa = log_kappa,
                                    v = c(1,2),
                                    log_sigma_epsilon = log_sigma_epsilon,
                                    log_sigma_u = log_sigma_u,
                                    lambda =lambda,
                                    lambda1 = lambda1,
                                    lambda_epsilon = lambda_epsilon,
                                    lambda_u = lambda_u,
                                    y= y,
                                    A = A,
                                    m_u = m_u
                                    )

map <- MAP(mesh = mesh,
    lambda =lambda, lambda1 = lambda1, lambda_epsilon = lambda_epsilon, lambda_u = lambda_u,
    y= y, A = A, m_u =m_u, maxiterations = 300, log_sigma_epsilon = log_sigma_epsilon )
print(map)
cov2cor(solve(-map$hessian))
par <- map$par
real_par <- c(log_kappa,v,log_sigma_epsilon)
print(par-real_par)

##Trying to see what doesn't work
sigma_u <- 0.01; log_sigma_u <- log(sigma_u)
sigma_epsilon <- 0.01 ; log_sigma_epsilon <- log(sigma_epsilon)


# Calculates log-prior
log_pc_aniso_value <- log_pc_prior_aniso(lambda = lambda, lambda1 = lambda1, log_kappa = log_kappa, v = v)
log_pc_noise_value <- log_pc_prior_noise_variance(lambda_epsilon = lambda_epsilon, log_sigma_epsilon = log_sigma_epsilon)
log_pc_sigma_u_value <- log_pc_prior_noise_variance(lambda_epsilon = lambda_u, log_sigma_epsilon = log_sigma_u)
log_pc_value <- log_pc_aniso_value + log_pc_noise_value + log_pc_sigma_u_value

# Calculates anisotropy
n <- nrow(mesh$loc)
kappa_values <- rep(kappa, n)
vec_values <- matrix(v, n, 2, byrow = TRUE)
aniso <- list(kappa = kappa_values, vec = vec_values)

# Calculates log-density of the distribution of u at m_u knowing (kappa, v)
Q_u <- fm_aniso_precision(mesh, aniso, log_sigma = log_sigma_u)
if (length(m_u) == 1) {
  m_u <- rep(m_u, n)
}
u <- m_u
logGdty_prior <- logGdensity(x = u, mu = m_u, Q = Q_u)

# Calculates Q_epsilon,  Q_{u|y,theta} and m_{u|y,theta}
Q_epsilon <- Matrix::Diagonal(n, sigma_epsilon^2)
Q_uy_theta <- Q_u + t(A) %*% Q_epsilon %*% A
m_uy_theta <- solve(Q_uy_theta, Q_u %*% m_u + t(A) %*% Q_epsilon %*% y)

# Calculates  log-density of the posterior distribution of u given y and theta
logGdty_posterior <- logGdensity(x = u, mu = m_uy_theta, Q = Q_uy_theta)

# Calculates  log-density of the observation of y given u, theta
logGdty_observation <- logGdensity(x = y, mu = A %*% u, Q = Matrix::Diagonal(n, 1))

# Calculates  log-posterior
log_posterior_val <- log_pc_value + logGdty_prior + logGdty_observation - logGdty_posterior



ggplot()+ gg(data=mesh,color = x)
ggplot()+ gg(data=mesh,color = sqrt(as.vector(diag(INLA::inla.qinv(Q)))))
ggplot()+ gg(data=mesh,color = as.vector(diag(INLA::inla.qinv(Q))))

