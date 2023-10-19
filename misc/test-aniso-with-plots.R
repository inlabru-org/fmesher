library(devtools)
setwd("C:/Users/liaml/OneDrive/Documentos")
load_all("fmesher")
library(ggplot2)

# Parameter values stationary, anisotropic field
kp <- 1
v1 <- 0
v2 <- 5
lambda <- exp(sqrt(v1^2 + v2^2))
stretch <- sqrt(lambda)

# Parameter values Matérn field
nu <- 2 - 2 / 2
rh <- sqrt(8 * nu) / kp


# Kappa and vector field
kappa <- function(x) {
  return(kp)
}

vec <- function(x) {
  return(c(0, -1 ))
}

# Square mesh for field
l <- 4
library(sf)
boundary_sf <- st_sfc(st_polygon(list(rbind(c(-l, -l), c(l, -l), c(l, l), c(-l, l), c(-l, -l)))))
boundary <- fm_as_segm(boundary_sf)

mesh <- fm_mesh_2d_inla(
  boundary = boundary,
  max.edge = c(0.2, 1)
)
nodes1 <- mesh$loc
#plot(mesh)

# Defining anisotropy
kappa_values <- apply(nodes1, 1, kappa)
vec_values <- t(apply(nodes1, 1, vec))
aniso <- list(kappa_values, vec_values)

# Sample of Matérn field u' with correlation range kp
sample_matern <- fm_matern_sample(mesh, alpha = 2, rho = rh, sigma = 1)

# Sample of anisotropic field u, should be equal to  u'(H^{-1/2}x)
sample_aniso <- fm_aniso_sample(mesh, aniso)
sample_matern <- fm_matern_sample(mesh, alpha = 2, rho = rh, sigma = 1)

# Data for plotting
field_matern <- data.frame(
  x = nodes1[, 1],
  y = nodes1[, 2],
  u = sample_matern
)
field_aniso <- data.frame(
  x = nodes1[, 1],
  y = nodes1[, 2],
  u = sample_aniso
)

# Definining pixels for plotting
pxl <- fm_pixels(mesh, mask = boundary_sf)
pxl$uisotropic <- fm_evaluate(mesh,
  loc = pxl,
  field = field_matern$u
)
pxl$uaniso <- fm_evaluate(mesh,
  loc = pxl,
  field = field_aniso$u
)

# Plotting Matérn field
ggplot(pxl) +
  geom_tile(aes(geometry = geometry, fill = uisotropic),
    stat = "sf_coordinates", alpha = 1
  ) +
  scale_fill_gradientn(colours = c("#0000FFFF", "#FFFFFFFF", "#FF0000FF"), limits = c(min(field_matern$u), max(field_matern$u))) +
  coord_equal() +
  xlab("X Coordinate") +
  ylab("Y Coordinate") +
  labs(fill = "Field u isotropic")

# Plotting anisotropic field
ggplot(pxl) +
  geom_tile(aes(geometry = geometry, fill = uaniso),
    stat = "sf_coordinates", alpha = 1
  ) +
  scale_fill_gradientn(colours = c("#0000FFFF", "#FFFFFFFF", "#FF0000FF"), limits = c(min(field_aniso$u), max(field_aniso$u))) +
  coord_equal() +
  xlab("X Coordinate") +
  ylab("Y Coordinate") +
  labs(fill = "Field u")
