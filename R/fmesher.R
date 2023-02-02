#' @importFrom methods as
# Explicit import of something from Matrix to appease automated checks:
#' @importFrom Matrix as.matrix

fm_as_dgCMatrix <- function(x) {
  if (inherits(x, "dgCMatrix")) {
    x
  } else {
    as(as(as(x, "dMatrix"), "generalMatrix"), "CsparseMatrix")
  }
}

fm_as_dgTMatrix <- function(x, unique = TRUE) {
  if (unique) {
    as(fm_as_dgCMatrix(x), "TsparseMatrix")
  } else {
    if (inherits(x, "dgTMatrix")) {
      x
    } else {
      as(as(as(x, "dMatrix"), "generalMatrix"), "TsparseMatrix")
    }
  }
}

#' Sparse partial inverse
#'
#' @param A A symmetric matrix
#'
#' @export

fm_qinv <- function(A) {
  A_C <- fm_as_dgCMatrix(A)
  stopifnot(nrow(A_C) == ncol(A_C))
  if (!identical(A_C, Matrix::t(A_C))) {
    warning("Asymmetric matrix A detected, but only lower left triangle will be used.")
  }
  C_qinv(A_C)
}






# Convert loc information to raw matrix coordinates for the mesh
fm_get_native_loc <- function(mesh, loc) {
  # Support INLA <= 22.11.27 by converting globes to spheres
  # TODO: Handle the > 22.11.27 more efficiently
  if (!fm_crs_is_null(mesh$crs)) {
    if (fm_crs_is_geocent(mesh$crs)) {
      crs.sphere <- fm_CRS("sphere")
      if (!fm_identical_CRS(mesh$crs, crs.sphere)) {
        ## Convert the mesh to a perfect sphere.
        mesh <- fm_transform(mesh, crs = crs.sphere)
      }
      if (!is.matrix(loc)) {
        if (!fm_identical_CRS(crs, crs.sphere)) {
          loc <- fm_transform(loc, crs = crs.sphere, crs0 = crs)
        }
      } else {
        if (fm_crs_is_null(crs)) {
          loc <- loc / rowSums(loc^2)^0.5
        } else {
          loc <- fm_transform(loc, crs = crs.sphere, crs0 = crs)
        }
      }
    } else if (!fm_identical_CRS(crs, mesh$crs)) {
      mesh <- fm_transform(mesh, crs = mesh$crs, crs0 = crs, passthrough = TRUE)
    }
  } else if (identical(mesh$manifold, "S2")) {
    mesh$loc <- mesh$loc / rowSums(mesh$loc^2)^0.5
    loc <- loc / rowSums(loc^2)^0.5
  }

  if (inherits(loc, c("SpatialPoints", "SpatialPointsDataFrame"))) {
    loc <- sp::coordinates(loc)
  } else if (inherits(loc, c("sf", "sfc", "sfg"))) {
    loc <- sf::st_coordinates(loc)
    c_names <- colnames(loc)
    c_names <- intersect(c_names, c("X", "Y", "Z"))
    loc <- loc[, c_names, drop = FALSE]
  }

  loc
}



#' @title Compute barycentric coordinates
#'
#' @description Identify triangles and compute barycentric coordinates
#'
#' @param mesh `inla.mesh` object
#' @param loc Points for which to identify the containing triangle, and
#' corresponding barycentric coordinates. May be raw matrix coordinates, `sf`, or `sp`
#' point information.
#'
#' @return A list with elements `t` (vector of triangle indices) and `bary`
#' (3-column matrix of barycentric coordinates). Points that were not found
#' give `NA` entries in `t` and `bary`.
#'
#' @export
fm_bary <- function(mesh, loc) {
  loc <- fm_get_native_loc(mesh, loc)
  jj <-
    which(rowSums(matrix(is.na(as.vector(loc)),
                         nrow = nrow(loc),
                         ncol = ncol(loc)
    )) == 0)
  result <- fmesher_bary(loc = loc[jj, , drop = FALSE],
                         mesh_loc = mesh$loc,
                         mesh_tv = mesh$graph$tv,
                         options = list())
  tri <- rep(NA_integer_, nrow(loc))
  bary <- matrix(NA_real_, nrow(loc), 3)
  tri[jj[result$t >= 0]] <- result$t[result$t >= 0] + 1L
  bary[jj[result$t > 0], ] <- result$bary[result$t >= 0, ]
  list(t = tri, bary = bary)
}
