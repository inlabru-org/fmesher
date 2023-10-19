library(devtools)
library(ggplot2)
library(Matrix)
library(INLA)
document()
load_all()

#Defining parameters
lambda <- 1
lambda1 <- 1
lambda_epsilon <- 1
lambda_u <- 1
kp <- 3
sigma_u <- 0.01
logkp <- log(kp)
log_sigma_u <- log(sigma_u)
log_sigma_epsilon <- log(0.01)
v <- c(1,2)
kappa <- function(x) {
  return(kp)
}

vec <- function(x) {
  return(v)
}

#Testing PC priors
log_pc_value <- log_pc_prior_aniso(lambda, lambda1, kp, v)
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
kappa_values <- apply(nodes, 1, kappa)
vec_values <- t(apply(nodes, 1, vec))
aniso=list(kappa_values,vec_values)

#Calculating log determinant
scaling <- 1/(4 * pi * sigma_u ^2)
Q <- scaling * fm_aniso_precision(mesh,aniso)
log_det <- sparse_log_determinant_chol(Q)
#Checking value
log_det_2 <- determinant(Q,logarithm = TRUE)
log_det_2_v <- as.numeric(log_det_2$modulus)
log_det - log_det_2_v

#Calculating |x|_Q
x <- rep(2,n)
x_norm <- norm_Q(Q,x)

#Calculating logGdty
m_u = as.vector(rep(0,n))
logdtyQ <- logGdensity(x,m_u,Q) #This should be small

#Sampling from noisy data
x = fm_aniso_basis_weights_sample(x = mesh, aniso = aniso, log_sigma = log_sigma_u)
A = Matrix::Diagonal(n, 1)
y = A %*% x + exp(log_sigma_epsilon)*stats::rnorm(nrow(Q))

#Calculate log posterior and map
log_posterior(mesh = mesh, log_kappa = logkp,
              log_sigma_epsilon = log_sigma_epsilon, log_sigma_u = log_sigma_u,
              v = v, lambda =lambda, lambda1 = lambda1, lambda_epsilon = lambda_epsilon,
              lambda_u = lambda_u,
              y= y, A = A, m_u =m_u )


map <- MAP(mesh = mesh,
    lambda =lambda, lambda1 = lambda1, lambda_epsilon = lambda_epsilon, lambda_u = lambda_u,
    y= y, A = A, m_u =m_u, maxiiterations = 3 )
print(map)
cov2cor(solve(-map$hessian))
map$par

ggplot()+ gg(data=mesh,color = x)
ggplot()+ gg(data=mesh,color = sqrt(as.vector(diag(INLA::inla.qinv(Q)))))
ggplot()+ gg(data=mesh,color = as.vector(diag(INLA::inla.qinv(Q))))
mean(diag(INLA::inla.qinv(Q)))

