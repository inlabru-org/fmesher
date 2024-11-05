## code to prepare `fmexample` dataset goes here
## Also see the fmexample_sp() function

library(fmesher)

set.seed(1234L)
loc <- matrix(rnorm(20), 10, 2)
loc_sf <- sf::st_geometry(sf::st_as_sf(as.data.frame(loc), coords = 1:2))
boundary_sf <- fm_extensions(loc, c(1, 3))

mesh <- fm_mesh_2d_inla(boundary = boundary_sf, max.edge = c(0.515, 2))
plot(mesh)

fmexample <- list(
  loc = loc,
  loc_sf = loc_sf,
  boundary_fm = fm_as_segm_list(boundary_sf),
  boundary_sf = boundary_sf,
  mesh = mesh
)

# usethis::use_data(fmexample, overwrite = TRUE)
