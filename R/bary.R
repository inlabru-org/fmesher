#' @include deprecated.R

# fm_bary ####

#' @title Compute barycentric coordinates
#'
#' @description Identify knot intervals or triangles and compute barycentric coordinates
#'
#' @param mesh `fm_mesh_1d` or `fm_mesh_2d` object
#' @param loc Points for which to identify the containing triangle, and
#' corresponding barycentric coordinates. May be a vector (for 1d) or
#' raw matrix coordinates, `sf`, or `sp` point information (for 2d).
#' @param \dots Arguments forwarded to sub-methods.
#'
#' @export
fm_bary <- function(mesh, loc, ...) {
  UseMethod("fm_bary")
}


## Binary split method, returning the index of the left knot for the
## interval containing each location. Points to the left are assigned index 1,
## and points to the right are assigned index length(knots)-1.
do.the.split <- function(knots, loc) {
  n <- length(knots)
  if (n <= 2L) {
    return(rep(1L, length(loc)))
  }
  split <- 1L + (n - 1L) %/% 2L ## Split point
  upper <- (loc >= knots[split])
  idx <- rep(0, length(loc))
  idx[!upper] <- do.the.split(knots[1:split], loc[!upper])
  idx[upper] <- split - 1L + do.the.split(knots[split:n], loc[upper])
  return(idx)
}



#' @describeIn fm_bary Return a list with elements
#' `t` (start and endpoint knot indices) and `bary` (barycentric coordinates), both
#' 2-column matrices. For backwards compatibility with old inla code, a copy `index=t`
#' is also included in the list.
#'
#' For `method = "nearest"`, `t[,1]` contains the index of the nearest mesh knot,
#' and each row of `bary` contains `c(1, 0)`.
#' @param method character; method for defining the barycentric coordinates,
#' "linear" (default) or "nearest"
#' @param restricted logical, used for `method="linear"`.
#' If `FALSE` (default), points outside the mesh interval will be given
#' barycentric weights less than 0 and greater than 1, according to linear
#' extrapolation. If `TRUE`, the barycentric weights are clamped to the (0, 1)
#' interval.
#' @export
fm_bary.fm_mesh_1d <- function(mesh,
                               loc,
                               method = c("linear", "nearest"),
                               restricted = FALSE, ...) {
  method <- match.arg(method)

  if (mesh$cyclic) {
    knots <- c(mesh$loc - mesh$loc[1], diff(mesh$interval))
    loc <- (loc - mesh$loc[1]) %% diff(mesh$interval)
  } else {
    knots <- mesh$loc - mesh$loc[1]
    loc <- loc - mesh$loc[1]
  }

  idx <- do.the.split(knots, loc)
  u <- (loc - knots[idx]) / (knots[idx + 1L] - knots[idx])

  if (method == "nearest") {
    if (mesh$cyclic) {
      idx <- idx + (u > 0.5)
      u <- numeric(length(loc))
      idx <- (idx - 1L) %% mesh$n + 1L
      idx_next <- idx %% mesh$n + 1L
    } else { # !cyclic
      idx <- idx + (u > 0.5)
      idx_next <- idx + 1L
      u <- numeric(length(loc))
      found <- (idx == mesh$n)
      idx_next[found] <- mesh$n - 1L
      u[found] <- 0.0
    }
  } else { ## (method=="linear") {
    if (mesh$cyclic) {
      idx_next <- idx %% mesh$n + 1L
    } else { # !cyclic
      idx_next <- idx + 1L
      if (restricted) {
        u[u < 0.0] <- 0.0
        u[u > 1.0] <- 1.0
      }
    }
  }

  index <- cbind(idx, idx_next)
  bary <- cbind(1 - u, u)

  return(list(t = index, bary = bary, index = index))
}


#' @param crs Optional crs information for `loc`
#'
#' @describeIn fm_bary A list with elements `t` (vector of triangle indices) and `bary`
#' (3-column matrix of barycentric coordinates). Points that were not found
#' give `NA` entries in `t` and `bary`.
#'
#' @export
fm_bary.fm_mesh_2d <- function(mesh, loc, crs = NULL, ...) {
  loc <- fm_onto_mesh(mesh, loc, crs = crs)

  # Avoid sphere accuracy issues by scaling to unit sphere
  scale <- 1
  if (fm_manifold(mesh, "S2")) {
    scale <- 1 / mean(rowSums(mesh$loc^2)^0.5)
    loc <- loc / rowSums(loc^2)^0.5
  }

  pre_ok_idx <-
    which(rowSums(matrix(
      is.na(as.vector(loc)),
      nrow = nrow(loc),
      ncol = ncol(loc)
    )) == 0)
  result <- fmesher_bary(
    mesh_loc = mesh$loc * scale,
    mesh_tv = mesh$graph$tv - 1L,
    loc = loc[pre_ok_idx, , drop = FALSE],
    options = list()
  )
  tri <- rep(NA_integer_, nrow(loc))
  bary <- matrix(NA_real_, nrow(loc), 3)
  ok <- result$t >= 0
  tri[pre_ok_idx[ok]] <- result$t[ok] + 1L
  bary[pre_ok_idx[ok], ] <- result$bary[ok, ]
  list(t = tri, bary = bary)
}



#' @rdname fm_bary
#' @export
#' @method fm_bary inla.mesh
fm_bary.inla.mesh <- function(mesh, ...) {
  fm_bary.fm_mesh_2d(fm_as_mesh_2d(mesh), ...)
}

#' @rdname fm_bary
#' @export
#' @method fm_bary inla.mesh.1d
fm_bary.inla.mesh.1d <- function(mesh, ...) {
  fm_bary.fm_mesh_1d(fm_as_mesh_1d(mesh), ...)
}
