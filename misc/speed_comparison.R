suppressPackageStartupMessages(library(inlabru))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(INLA))
suppressPackageStartupMessages(library(fmesher))

#### Mesh creation ####

# Some seed points

set.seed(12345L)
loc <- matrix(rnorm(500), 250, 2)

# Basic case

extend <- list(offset = 1)
refine <- list(max.edge = 1)

fmesher_INLA <- inla.mesh.create(loc = loc, extend = extend, refine = refine)
fmesher_Rcpp <- fm_rcdt_2d_inla(loc = loc, extend = extend, refine = refine)
(ggplot() +
  gg(fmesher_INLA) +
  ggtitle(paste0("fmesher/INLA: nV =  ", fmesher_INLA$n)) |
  ggplot() +
    gg(fmesher_Rcpp) +
    ggtitle(paste0("fmesher/Rcpp: nV =  ", fmesher_Rcpp$n)))

bench::mark(
  INLA = inla.mesh.create(loc = loc, extend = extend, refine = refine),
  Rcpp = fm_rcdt_2d_inla(loc = loc, extend = extend, refine = refine),
  check = FALSE
)

# Basic case with more triangles

extend <- list(offset = 1)
refine <- list(max.edge = 0.2)

fmesher_INLA <- inla.mesh.create(loc = loc, extend = extend, refine = refine)
fmesher_Rcpp <- fm_rcdt_2d_inla(loc = loc, extend = extend, refine = refine)
(ggplot() +
  gg(fmesher_INLA) +
  ggtitle(paste0("fmesher/INLA: nV =  ", fmesher_INLA$n)) |
  ggplot() +
    gg(fmesher_Rcpp) +
    ggtitle(paste0("fmesher/Rcpp: nV =  ", fmesher_Rcpp$n)))

bench::mark(
  INLA = inla.mesh.create(loc = loc, extend = extend, refine = refine),
  Rcpp = fm_rcdt_2d_inla(loc = loc, extend = extend, refine = refine),
  check = FALSE
)

# Two-stage construction:

offset <- c(1, 2)
max.edge <- c(0.2, 1)

fmesher_INLA <- inla.mesh.2d(loc = loc, offset = offset, max.edge = max.edge)
fmesher_Rcpp <- fm_mesh_2d_inla(loc = loc, offset = offset, max.edge = max.edge)
(ggplot() +
  gg(fmesher_INLA) +
  ggtitle(paste0("fmesher/INLA: nV =  ", fmesher_INLA$n)) |
  ggplot() +
    gg(fmesher_Rcpp) +
    ggtitle(paste0("fmesher/Rcpp: nV =  ", fmesher_Rcpp$n)))

bench::mark(
  INLA = inla.mesh.2d(loc = loc, offset = offset, max.edge = max.edge),
  Rcpp = fm_mesh_2d_inla(loc = loc, offset = offset, max.edge = max.edge),
  check = FALSE
)


#### Basis evaluation ####

offset <- c(1, 2)
max.edge <- c(0.2, 1)

fmesher_INLA <- inla.mesh.2d(loc = loc, offset = offset, max.edge = max.edge)
fmesher_Rcpp <- fm_mesh_2d_inla(loc = loc, offset = offset, max.edge = max.edge)

set.seed(12345L)
loc_find <- list(
  matrix(rnorm(50), 25, 2),
  matrix(rnorm(500), 250, 2),
  matrix(rnorm(5000), 2500, 2),
  matrix(rnorm(50000), 25000, 2),
  matrix(rnorm(500000), 250000, 2)
)

for (lf in loc_find) {
  message(paste0("nrow(locations) = ", nrow(lf)))
  print(bench::mark(
    INLA = inla.mesh.projector(fmesher_INLA, loc = lf)$proj,
    Rcpp = fm_evaluator(fmesher_Rcpp, loc = lf)$proj,
    check = FALSE
  ))
}


#### profvis ####

lf <- loc_find[[3]]
system.time({
  for (k in seq_len(100)) {
    INLA <- inla.mesh.projector(fmesher_INLA, loc = lf)$proj
  }
})
system.time({
  for (k in seq_len(100)) {
    Rcpp <- fm_evaluator(fmesher_Rcpp, loc = lf)$proj
  }
})

if (FALSE) {
  profvis::profvis({
    for (k in seq_len(100)) {
      INLA <- inla.mesh.projector(fmesher_INLA, loc = lf)$proj
      Rcpp <- fm_evaluator(fmesher_Rcpp, loc = lf)$proj
    }
  })

  devtools::load_all() # For source
  lf <- loc_find[[3]]
  profvis::profvis({
    for (k in seq_len(100)) {
      Rcpp <- fm_evaluator(fmesher_Rcpp, loc = lf)$proj
    }
  })
}
