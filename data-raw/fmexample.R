## code to prepare `fmexample` dataset goes here

library(fmesher)

set.seed(1234L)
loc <- matrix(rnorm(20), 10, 2)
boundary <- fm_extensions(loc, c(1, 3))

mesh <- fm_mesh_2d_inla(boundary = boundary, max.edge = c(0.5, 2))
plot(mesh)

fmexample <- list(
  loc = loc,
  boundary = boundary,
  boundary_sp = lapply(boundary, sf::as_Spatial),
  mesh = mesh
)

usethis::use_data(fmexample, overwrite = TRUE)
