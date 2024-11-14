#' @include deprecated.R

# fm_bary ####

#' @title Compute barycentric coordinates
#'
#' @description Identify knot intervals or triangles and compute barycentric
#'   coordinates
#'
#' @param mesh `fm_mesh_1d` or `fm_mesh_2d` object
#' @param loc Points for which to identify the containing interval/triangle, and
#'   corresponding barycentric coordinates. May be a vector (for 1d) or a matrix
#'   of raw coordinates, `sf`, or `sp` point information (for 2d).
#' @param \dots Arguments forwarded to sub-methods.
#' @returns A `fm_bary` object, a `tibble` with columns `index`; either
#' \itemize{
#' \item{vector of triangle indices (triangle meshes),}
#' \item{matrix of interval knot indices (1D meshes), or}
#' \item{matrix of lower left box indices (2D lattices),}
#' }
#' and `where`, a matrix of barycentric coordinates.
#'
#' @seealso [fm_bary_simplex()]
#' @export
#' @examples
#' str(fm_bary(fmexample$mesh, fmexample$loc_sf))
#' str(fm_bary(fm_mesh_1d(1:4), seq(0, 5, by = 0.5)))
fm_bary <- function(...) {
  UseMethod("fm_bary")
}

#' @describeIn fm_bary Returns the `bary` input unchanged
#' @param bary An `fm_bary` object, or an object that can be converted to
#' `fm_bary`.
#' @export
fm_bary.fm_bary <- function(bary, ...) {
  bary
}

#' @describeIn fm_bary Converts a `list` `bary` to `fm_bary`.
#' In the list elements are unnamed, the names `index` and `where` are assumed.
#' @export
fm_bary.list <- function(bary, ...) {
  if (is.null(names(bary))) {
    names(bary) <- c("index", "where")
  }
  bary <- tibble::tibble(
    index = bary[["index"]],
    where = bary[["where"]]
  )
  storage.mode(bary[["index"]]) <- "integer"
  structure(
    bary,
    class = c("fm_bary", class(bary))
  )
}

#' @describeIn fm_bary Converts a [tibble::tibble()] `bary` to `fm_bary`
#' @export
fm_bary.tbl_df <- function(bary, ...) {
  stopifnot(
    all(c("index", "where") %in% names(bary))
  )
  storage.mode(bary[["index"]]) <- "integer"
  structure(
    bary,
    class = c("fm_bary", class(bary))
  )
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



#' @describeIn fm_bary Return an `fm_bary` object with elements `index`
#'   (starting knot indices indices) and `where` (barycentric coordinates), both
#'   2-column matrices. Use [fm_bary_simplex()] to obtain the corresponding
#'   endpoint knot indices.
#'
#'   For `method = "nearest"`, `index` contains the index of the nearest mesh
#'   knot, and `where` is a single-column all-ones matrix.
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
      idx <- (idx - 1L) %% mesh$n + 1L
    } else { # !cyclic
      idx <- idx + (u > 0.5)
      idx_next <- idx + 1L
    }
    bary <- matrix(1.0, length(loc), 1)
  } else { ## (method=="linear") {
    if (!mesh$cyclic && restricted) {
      u[u < 0.0] <- 0.0
      u[u > 1.0] <- 1.0
    }
    bary <- cbind(1 - u, u, deparse.level = 0)
  }

  fm_bary(
    tibble::tibble(
      index = idx,
      where = bary
    )
  )
}


#' @describeIn fm_bary An `fm_bary` object with columns `index` (vector of
#'   triangle indices) and `where` (3-column matrix of barycentric coordinates).
#'   Points that were not found give `NA` entries in `index` and `where`.
#' @param crs Optional crs information for `loc`
#' @param max_batch_size integer; maximum number of points to process in a
#'   single batch. This speeds up calculations by avoiding repeated large
#'   internal memory allocations and data copies. The default, `NULL`, uses
#'   `max_batch_size = 2e5L`, chosen based on empirical time measurements to
#'   give an approximately optimal runtime.
#'
#' @export
fm_bary.fm_mesh_2d <- function(mesh,
                               loc,
                               crs = NULL,
                               ...,
                               max_batch_size = NULL) {
  if (is.null(max_batch_size)) {
    max_batch_size <- 2e5L
  }

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
  if (length(pre_ok_idx) <= max_batch_size) {
    result <- fmesher_bary(
      mesh_loc = mesh$loc * scale,
      mesh_tv = mesh$graph$tv - 1L,
      loc = loc[pre_ok_idx, , drop = FALSE],
      options = list()
    )
    tri <- rep(NA_integer_, nrow(loc))
    where <- matrix(NA_real_, nrow(loc), 3)
    ok <- result$index >= 0
    tri[pre_ok_idx[ok]] <- result$index[ok] + 1L
    where[pre_ok_idx[ok], ] <- result$where[ok, ]
  } else {
    tri <- rep(NA_integer_, nrow(loc))
    where <- matrix(NA_real_, nrow(loc), 3)
    n_batches <- ceiling(length(pre_ok_idx) / max_batch_size)
    batch_idx <- round(seq(0, length(pre_ok_idx), length.out = n_batches + 1))
    subindex <- split(pre_ok_idx, rep(seq_len(n_batches), diff(batch_idx)))
    for (k in seq_along(subindex)) {
      result <- fmesher_bary(
        mesh_loc = mesh$loc * scale,
        mesh_tv = mesh$graph$tv - 1L,
        loc = loc[subindex[[k]], , drop = FALSE],
        options = list()
      )
      ok <- result$index >= 0
      tri[subindex[[k]][ok]] <- result$index[ok] + 1L
      where[subindex[[k]][ok], ] <- result$where[ok, ]
    }
  }

  fm_bary(
    tibble::tibble(
      index = tri,
      where = where
    )
  )
}


#' @title Extract Simplex information for Barycentric coordinates
#'
#' @description
#' Extract the simplex vertex information for a combination of a mesh
#' and `fm_bary` coordinates.
#'
#' @param mesh A mesh object, e.g. [fm_mesh_2d] or [fm_mesh_1d].
#' @param bary An `fm_bary` object. If NULL, return the full simplex
#' information for the mesh.
#' @param \dots Further arguments potentially used by sub-methods.
#' @returns A matrix of vertex indices, one row per point in `bary`.
#' @seealso [fm_bary()]
#' @export
fm_bary_simplex <- function(mesh, bary = NULL, ...) {
  UseMethod("fm_bary_simplex")
}

#' @describeIn fm_bary_simplex Extract the triangle vertex indices for a 2D mesh
#' @export
#'
#' @examples
#' bary <- fm_bary(fmexample$mesh, fmexample$loc_sf)
#' fm_bary_simplex(fmexample$mesh, bary)
fm_bary_simplex.fm_mesh_2d <- function(mesh, bary = NULL, ...) {
  if (is.null(bary)) {
    return(mesh$graph$tv)
  }
  if (NROW(bary) == 0L) {
    return(matrix(integer(1), 0L, 3L))
  }
  mesh$graph$tv[bary$index, , drop = FALSE]
}

#' @describeIn fm_bary_simplex Extract the edge vertex indices for a 1D mesh
#'
#' @export
#' @examples
#' mesh1 <- fm_mesh_1d(1:4)
#' (bary1 <- fm_bary(mesh1, seq(0, 5, by = 0.5)))
#' (bary1 <- fm_bary(mesh1, seq(0, 5, by = 0.5), restricted = TRUE))
#' fm_bary_simplex(mesh1, bary1)
fm_bary_simplex.fm_mesh_1d <- function(mesh, bary = NULL, ...) {
  if (is.null(bary)) {
    if (mesh$cyclic) {
      return(cbind(seq_len(mesh$n), seq_len(mesh$n) %% mesh$n + 1L))
    }
    return(cbind(seq_len(mesh$n - 1L), seq_len(mesh$n - 1L) + 1L))
  }
  if (NCOL(bary$where) == 1L) {
    return(matrix(bary$index, NROW(bary), 1L))
  }
  if (NROW(bary) == 0L) {
    return(matrix(integer(1), 0L, 2L))
  }
  if (mesh$cyclic) {
    idx_next <- bary$index %% mesh$n + 1L
  } else { # !cyclic
    idx_next <- bary$index + 1L
  }
  cbind(bary$index, idx_next, deparse.level = 0)
}
