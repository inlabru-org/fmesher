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
fm_onto_mesh <- function(mesh, loc, crs = NULL) {
  if (!is.matrix(loc) && !fm_crs_is_null(crs)) {
    warning("loc is non-matrix but crs specified; will be ignored")
  }
  if (inherits(loc, c("SpatialPoints", "SpatialPointsDataFrame",
                      "sf", "sfc", "sfg"))) {
    crs <- fm_crs(loc)
  }
  mesh_crs <- fm_crs(mesh)

  loc_needs_normalisation <- FALSE
  if (!fm_crs_is_null(crs) && !fm_crs_is_null(mesh_crs)) {
    if (fm_crs_is_geocent(mesh_crs)) {
      if (!is.matrix(loc)) {
        if (!fm_identical_CRS(crs, mesh_crs)) {
          loc <- fm_transform(loc, crs = mesh_crs, crs0 = crs)
        }
      }
    } else if (!fm_identical_CRS(crs, mesh_crs)) {
      loc <- fm_transform(loc, crs = mesh_crs, crs0 = crs, passthrough = TRUE)
    }
  } else if (identical(mesh$manifold, "S2")) {
    loc_needs_normalisation <- TRUE
  }

  if (inherits(loc, c("SpatialPoints", "SpatialPointsDataFrame"))) {
    loc <- sp::coordinates(loc)
  } else if (inherits(loc, c("sf", "sfc", "sfg"))) {
    loc <- sf::st_coordinates(loc)
    c_names <- colnames(loc)
    c_names <- intersect(c_names, c("X", "Y", "Z"))
    loc <- loc[, c_names, drop = FALSE]
  } else if (!is.matrix(loc)) {
    warning(
      paste0(
        "Unclear if the 'loc' class ('",
        paste0(class(loc), collapse = "', '"),
        "') is of a type we know how to handle."
      ),
      immediate. = TRUE
    )
  }
  if (loc_needs_normalisation) {
    loc <- loc / rowSums(loc^2)^0.5 * mean(rowSums(mesh$loc^2)^0.5)
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
#' @param crs Optional crs information for `loc`
#'
#' @return A list with elements `t` (vector of triangle indices) and `bary`
#' (3-column matrix of barycentric coordinates). Points that were not found
#' give `NA` entries in `t` and `bary`.
#'
#' @export
fm_bary <- function(mesh, loc, crs = NULL) {
  # Support INLA <= 22.11.27 by converting globes to spheres
  # TODO: Handle the > 22.11.27 more efficiently
  if (fm_crs_is_geocent(fm_crs(mesh))) {
    crs.sphere <- fm_crs("sphere")
    if (!fm_identical_CRS(fm_crs(mesh), crs.sphere)) {
      ## Convert the mesh to a perfect sphere.
      mesh <- fm_transform(mesh, crs = crs.sphere)
    }
  }
  loc <- fm_onto_mesh(mesh, loc, crs = crs)

  # Avoid sphere accuracy issues by scaling to unit sphere
  scale <- 1
  if (identical(mesh$manifold, "S2")) {
    scale <- 1 / mean(rowSums(mesh$loc^2)^0.5)
    loc <- loc / rowSums(loc^2)^0.5
  }

  pre_ok_idx <-
    which(rowSums(matrix(is.na(as.vector(loc)),
                         nrow = nrow(loc),
                         ncol = ncol(loc)
    )) == 0)
  result <- fmesher_bary(loc = loc[pre_ok_idx, , drop = FALSE],
                         mesh_loc = mesh$loc * scale,
                         mesh_tv = mesh$graph$tv - 1L,
                         options = list())
  tri <- rep(NA_integer_, nrow(loc))
  bary <- matrix(NA_real_, nrow(loc), 3)
  ok <- result$t >= 0
  tri[pre_ok_idx[ok]] <- result$t[ok] + 1L
  bary[pre_ok_idx[ok], ] <- result$bary[ok, ]
  list(t = tri, bary = bary)
}
