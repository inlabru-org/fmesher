library(devtools)
library(ggplot2)
library(Matrix)
library(sp)
library(INLA)

#Hyperparameters for PC priors
lambda <- 1; lambda1 <- 1; lambda_epsilon <- 1; lambda_u <- 1

#Anisotropy parameters
kappa <- 1; logkappa <- log(kappa);
v <- c(1,2)

#Correlation range calculation
rho <- sqrt(8)/kappa/sqrt(exp(sqrt(v[1]^2 + v[2]^2)))

#Noise parameters
sigma_u <- 0.01; log_sigma_u <- log(sigma_u)
sigma_epsilon <- 0.01; log_sigma_epsilon <- log(sigma_epsilon)


#Testing PC priors
log_pc_aniso <- log_pc_prior_aniso(lambda, lambda1, kappa, v)
log_pc_noise <- log_pc_prior_noise_variance(lambda_epsilon = lambda_epsilon, log_sigma_epsilon = log_sigma_epsilon)
log_pc_u <- log_pc_prior_noise_variance(lambda_epsilon = lambda_u, log_sigma_epsilon = log_sigma_u)

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
                                    log_kappa = logkappa,
                                    v = v,
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
    y= y, A = A, m_u =m_u, maxiterations = 20, log_sigma_epsilon = log_sigma_epsilon )
print(map)
cov2cor(solve(-map$hessian))
par <- map$par


ggplot()+ gg(data=mesh,color = x)
ggplot()+ gg(data=mesh,color = sqrt(as.vector(diag(INLA::inla.qinv(Q)))))
ggplot()+ gg(data=mesh,color = as.vector(diag(INLA::inla.qinv(Q))))

