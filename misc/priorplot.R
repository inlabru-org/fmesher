# Load required packages
library(Matrix)
library(pracma)


# Define Matern covariance function
sigma <- function(kappa, alpha, d) {
  sqrt(gamma(alpha - d / 2) / (gamma(alpha) * (4 * pi)^(d / 2) * kappa^(2 * (alpha - d / 2))))
}

K <- function(x, kappa, alpha, d) {
  nu = alpha - d / 2
  sigma(kappa, alpha, d)^2 / (2^(nu - 1) * gamma(nu)) * (kappa * norm(x,type = "2"))^nu * besselK(kappa * norm(x,type = "2"), nu)
}

# Define deformed metric
H <- function(v1, v2) {
  matrix(c(cosh(sqrt(v1^2 + v2^2)) + (v1 * sinh(sqrt(v1^2 + v2^2))) / sqrt(v1^2 + v2^2),
           (v2 * sinh(sqrt(v1^2 + v2^2))) / sqrt(v1^2 + v2^2),
           (v2 * sinh(sqrt(v1^2 + v2^2))) / sqrt(v1^2 + v2^2),
           cosh(sqrt(v1^2 + v2^2)) - (v1 * sinh(sqrt(v1^2 + v2^2))) / sqrt(v1^2 + v2^2)), nrow = 2)
}

deformednorm <- function(kappa, v, x) {
  Hmat <- H(v[1], v[2])
  Fmatinv <- sqrtm(Hmat)$Binv
  kappa * norm(Fmatinv%*% x, type = "2")
}

# Define covariance function of sulution to SPDE with stationary diffusion H_v and inverse correlation kappa
covarianceH <- function(kappa, v, x) {
  nu = 1
  sigma(kappa, 2, 2)^2 / (2^(nu - 1) * gamma(nu)) * (kappa * norm(x, type = "2"))^nu * besselK(kappa * deformednorm(kappa, v, x), nu)
}

# Create a density plot
kappa <- 1
v <- c(-1,0)
l <- 4
pxl <- expand.grid(x = seq(-l, l, length.out = 100), y = seq(-l, l, length.out = 100))
pxl$uaniso <- mapply(function(x, y) covarianceH(kappa, v, c(x, y)), pxl$x, pxl$y)

# Plot
library(ggplot2)
ggplot(pxl, aes(x = x, y = y, fill = uaniso)) +
  geom_tile() +
  scale_fill_gradientn(colors = c("blue", "white", "red")) +
  labs(x = "X Coordinate", y = "Y Coordinate", title = "Covariance u stationary anisotropic")

