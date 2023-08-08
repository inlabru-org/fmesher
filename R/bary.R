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


#' @describeIn fm_bary Return a list with elements
#' `t` (start and endpoint knot indices) and `bary` (barycentric coordinates), both
#' 2-column matrices. For backwards compatibility with old inla code, a copy `index=t`
#' is also included in the list.
#' @param method character; method for defining the barycentric coordinates
#' @export
fm_bary.fm_mesh_1d <- function(mesh, loc, method = c("linear", "nearest"), ...) {
  method <- match.arg(method)

  if (method == "linear") {
    if (mesh$cyclic) {
      mloc <- c(mesh$loc - mesh$loc[1], diff(mesh$interval))
      loc <- (loc - mesh$loc[1]) %% diff(mesh$interval)
    } else {
      mloc <- c(mesh$loc - mesh$loc[1], diff(mesh$interval))
      loc <- pmax(0, pmin(diff(mesh$interval), loc - mesh$loc[1]))
    }
  } else {
    if (mesh$cyclic) {
      mloc <-
        c(
          mesh$loc[mesh$n] - diff(mesh$interval),
          mesh$loc,
          diff(mesh$interval)
        )
      mloc <- (mloc[-(mesh$n + 2)] + mloc[-1]) / 2
      loc <- (loc - mloc[1]) %% diff(mesh$interval)
      mloc <- mloc - mloc[1]
    } else {
      mloc <-
        c(
          0,
          (mesh$loc[1:(mesh$n - 1L)] +
            mesh$loc[2:mesh$n]) / 2 - mesh$loc[1],
          diff(mesh$interval)
        )
      loc <- pmax(0, pmin(diff(mesh$interval), loc - mesh$loc[1]))
    }
  }

  ## Binary split method:
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

  idx <- do.the.split(mloc, loc)

  if (method == "nearest") {
    u <- rep(0, length(loc))
    if (mesh$cyclic) {
      found <- which(idx == (mesh$n + 1L))
      idx[found] <- 1L
    }
  } else { ## (method=="linear") {
    u <- pmax(0, pmin(1, (loc - mloc[idx]) / (mloc[idx + 1L] - mloc[idx])))
    if (!mesh$cyclic) {
      found <- which(idx == mesh$n)
      idx[found] <- mesh$n - 1L
      u[found] <- 1
    }
  }

  if (mesh$cyclic) {
    index <- matrix(c(idx, (idx %% mesh$n) + 1L), length(idx), 2)
    bary <- matrix(c(1 - u, u), length(idx), 2)
  } else {
    index <- matrix(c(idx, idx + 1L), length(idx), 2)
    bary <- matrix(c(1 - u, u), length(idx), 2)
  }

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
