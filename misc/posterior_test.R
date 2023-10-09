library(devtools)
library(ggplot2)
library(Matrix)
document()
load_all()
# Use document() and devtools::load_all() to load files
# Example usage
lambda <- 1
lambda1 <- 1
lambda_epsilon <- 1
kp <- 0.5
logkp <- log(kp)
v <- c(1,2)
kappa <- function(x) {
  return(kp)
}

vec <- function(x) {
  return(v)
}

#Testing PC priors
log_pc_value <- log_pc_prior_aniso(lambda, lambda1, kp, v)
log_pc_noise <- log_pc_prior_noise_variance(lambda_epsilon = lambda_epsilon, log_sigma_epsilon = 2)

#Mesh definition
library(sf)
boundary_sf = st_sfc(st_polygon(list(rbind(c(0, 0), c(10, 0), c(10, 10), c(0, 10),c(0,0)))))
boundary = fm_as_segm(boundary_sf)

mesh <- fm_mesh_2d_inla(  boundary = boundary,  max.edge = c(4, 4))
nodes <- mesh$loc
n <- mesh$n
plot(mesh)

#Defining anisotropy
kappa_values <- apply(nodes, 1, kappa)
vec_values <- t(apply(nodes, 1, vec))
aniso=list(kappa_values,vec_values)

#Calculating log determinant
Q <- fm_aniso_precision(mesh,aniso)
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
logGdensity(x,m_u,Q) #This should be small

#Calculate of log posterior
log_sigma_epsilon <- 1
y = rep(0,n)
A = Matrix::Diagonal(n, 1)
log_posterior(mesh = mesh, log_kappa = logkp,log_sigma_epsilon = log_sigma_epsilon, v = v,
              lambda =lambda, lambda1 = lambda1, lambda_epsilon = lambda_epsilon,
              y= y, A = A, m_u =m_u )


map <- MAP(mesh = mesh,
    lambda =lambda, lambda1 = lambda1, lambda_epsilon = lambda_epsilon,
    y= y, A = A, m_u =m_u )
print(map)




