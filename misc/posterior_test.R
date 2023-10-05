library(devtools)
library(ggplot2)
library(Matrix)
# load_all()
# Example usage
lambda <- 1
lambda1 <- 1
kp <- 0.5
v <- c(1,2)
kappa <- function(x) {
  return(kp)
}

vec <- function(x) {
  return(v)
}

#Testing PC priors
pc_value <- pc_prior(lambda, lambda1, kp, v)
log_pc_value <- log_pc_prior(lambda, lambda1, kp, v)
print(pc_value)
print(log(pc_value)==log_pc_value)

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
#Which version to use?
log_det_2 <- determinant(Q,logarithm = TRUE)
log_det_2_v <- as.numeric(log_det_2$modulus)
log_det - log_det_2_v

#Calculating |x|_Q
x <- rep(2,n)
norm_Q(Q,x)

#Calculating logGdty
m_u = as.vector(rep(0,n))
logGdensity(x,mu,Q) #This should be small
#Calculate log posterior
y = rep(1,n)
A = Matrix::Diagonal(n, 3)
Q_epsilon = Matrix::Diagonal(n, 1)

#For some reason this doesn't execute, had similar problem with norm function what causes it?
log_posterior(mesh, kp, v, lambda, lambda1, y, A, Q_epsilon, m_u )

#Calculate log-density of the distribution of u knowing (kappa, v)
Q_u <- fm_aniso_precision(mesh, aniso)
u <- m_u  #Could it be anything?
logGdty_prior <- logGdensity(Q_u, u, m_u)

# Calculate Q_{u|y,theta} and m_{u|y,theta}
Q_uy_theta <- Q_u + t(A) %*% Q_epsilon %*% A
m_uy_theta <- solve(Q_uy_theta, Q_u %*% m_u + t(A) %*% Q_epsilon %*% y)

# Calculate log-density of the posterior distribution of u given y and theta
logGdty_posterior <- logGdensity(Q_uy_theta, u, m_uy_theta)

# Calculate log-posterior
log_posterior_val <- log_pc_value + logGdty_prior + logGdty_posterior

