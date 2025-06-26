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

mesh_loc <- rbind(
  c(0, 0, 0),
  c(0, 1, 0),
  c(0, 1, 1),
  c(0, 0, 1),
  c(1, 0, 1),
  c(1, 0, 0),
  c(1, 1, 0),
  c(1, 1, 1)
)
mesh_tv <- rbind(
  c(1, 2, 8, 3),
  c(1, 3, 8, 4),
  c(1, 4, 8, 5),
  c(1, 5, 8, 6),
  c(1, 6, 8, 7),
  c(1, 7, 8, 2)
)
mesh <- fm_mesh_3d(loc = mesh_loc, tv = mesh_tv)

set.seed(1234L)
loc <- matrix(runif(30), 10, 3)

b <- fm_bary(mesh, loc)
fm_bary_loc(mesh, b) - loc

fmexample3d <- list(
  loc = loc,
  mesh = mesh
)

# usethis::use_data(fmexample3d, overwrite = TRUE)
