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
#' @seealso [fm_bary_simplex()], [fm_bary_loc()]
#'
#' @export
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
    result <- rep(1L, length(loc))
    result[is.na(loc)] <- NA_integer_
    return(result)
  }
  ok <- !is.na(loc)
  split <- 1L + (n - 1L) %/% 2L ## Split point
  upper <- (loc[ok] >= knots[split])
  idx <- rep(0, length(loc))
  idx[ok][!upper] <- do.the.split(knots[1:split], loc[ok][!upper])
  idx[ok][upper] <- split - 1L + do.the.split(knots[split:n], loc[ok][upper])
  idx[!ok] <- NA_integer_
  return(idx)
}



#' @describeIn fm_bary Return an `fm_bary` object with elements `index`
#'   (edge index vector pointing to the first knot of each edge) and
#'   `where` (barycentric coordinates,
#'   2-column matrices). Use [fm_bary_simplex()] to obtain the corresponding
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
#' @examples
#' bary <- fm_bary(fm_mesh_1d(1:4), seq(0, 5, by = 0.5))
#' bary
fm_bary.fm_mesh_1d <- function(mesh,
                               loc,
                               method = c("linear", "nearest"),
                               restricted = FALSE, ...) {
  method <- match.arg(method)

  if (inherits(loc, "fm_bary")) {
    if (method == "nearest") {
      if (ncol(loc$where) == 1L) {
        return(loc)
      }
      simplex <- fm_bary_simplex(mesh, loc)
      ok <- !is.na(loc$index)
      is_second <- loc$where[, 1] < loc$where[, 2]
      idx <- rep(NA_integer_, length(loc$index))
      idx[ok & !is_second] <- simplex[ok & !is_second, 1L]
      idx[ok & is_second] <- simplex[ok & is_second, 2L]
      return(fm_bary(list(index = idx, where = matrix(1.0, nrow(loc), 1L))))
    }
    if (ncol(loc$where) == 2L) {
      return(loc)
    }
    loc_ <- fm_bary_loc(mesh, loc)
    return(fm_bary(mesh, loc_, method = "linear"))
  }

  if (mesh$cyclic) {
    knots <- c(mesh$loc - mesh$loc[1], diff(mesh$interval))
    loc <- (loc - mesh$loc[1]) %% diff(mesh$interval)
  } else {
    knots <- mesh$loc - mesh$loc[1]
    loc <- loc - mesh$loc[1]
  }

  idx <- do.the.split(knots, loc)
  ok <- !is.na(idx)

  u <- numeric(length(loc))
  u[ok] <- (loc[ok] - knots[idx[ok]]) / (knots[idx[ok] + 1L] - knots[idx[ok]])

  if (method == "nearest") {
    idx[ok] <- idx[ok] + (u[ok] > 0.5)
    if (mesh$cyclic) {
      idx[ok] <- (idx[ok] - 1L) %% mesh$n + 1L
    }
    bary <- matrix(1.0, length(loc), 1)
    bary[!ok, 1L] <- NA_real_
  } else { ## (method=="linear") {
    if (!mesh$cyclic && restricted) {
      u[ok][u[ok] < 0.0] <- 0.0
      u[ok][u[ok] > 1.0] <- 1.0
    }
    u[!ok] <- NA_real_
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
#' @examples
#' str(fm_bary(fmexample$mesh, fmexample$loc_sf))
fm_bary.fm_mesh_2d <- function(mesh,
                               loc,
                               crs = NULL,
                               ...,
                               max_batch_size = NULL) {
  if (inherits(loc, "fm_bary")) {
    return(loc)
  }

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


# Simplex extraction ####

#' @title Extract Simplex information for Barycentric coordinates
#'
#' @description
#' Extract the simplex vertex information for a combination of a mesh
#' and [fm_bary] coordinates.
#'
#' @param mesh A mesh object, e.g. [fm_mesh_2d] or [fm_mesh_1d].
#' @param bary An [fm_bary] object. If NULL, return the full simplex
#' information for the mesh.
#' @param \dots Further arguments potentially used by sub-methods.
#' @returns A matrix of vertex indices, one row per point in `bary`.
#' @seealso [fm_bary()], [fm_bary_loc()]
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


# Location extraction ####

#' @title Extract Euclidean Sgeometry from Barycentric coordinates
#'
#' @description
#' Extract the Euclidean coordinates for location identified by an [fm_bary]
#' object. This acts as the inverse of `fm_bary()`.
#'
#' @param mesh A mesh object, e.g. [fm_mesh_2d] or [fm_mesh_1d].
#' @param bary An `fm_bary` object. If `NULL`, return the mesh nodes is the mesh
#' class supports it, otherwise gives an error.
#' @param \dots Further arguments potentially used by sub-methods.
#' @param format Optional format for the output. If `NULL`, the output format
#' is determined by the default for the mesh object.
#' @returns Output format depends on the mesh `class`.
#' @seealso [fm_bary()], [fm_bary_simplex()]
#' @export
fm_bary_loc <- function(mesh, bary = NULL, ..., format = NULL) {
  UseMethod("fm_bary_loc")
}

#' @describeIn fm_bary_loc Extract points on a triangle mesh. Implemented
#' formats are `"matrix"` (default) and `"sf"`.
#' @export
#'
#' @examples
#' head(fm_bary_loc(fmexample$mesh))
#' bary <- fm_bary(fmexample$mesh, fmexample$loc_sf)
#' fm_bary_loc(fmexample$mesh, bary, format = "matrix")
#' fm_bary_loc(fmexample$mesh, bary, format = "sf")
fm_bary_loc.fm_mesh_2d <- function(mesh, bary = NULL, ..., format = NULL) {
  format = match.arg(format, c("matrix", "sf"))
  if (is.null(bary)) {
    loc <- mesh$loc
  } else if (NROW(bary) == 0L) {
    loc <- matrix(0.0, 0L, ncol(mesh$loc))
  } else {
    loc <- matrix(NA_real_, NROW(bary), ncol(mesh$loc))
    ok <- !is.na(bary$index)
    simplex <- fm_bary_simplex(mesh, bary = bary[ok, ])
    loc[ok, ] <- (mesh$loc[simplex[, 1L], , drop = FALSE] * bary$where[ok, 1] +
              mesh$loc[simplex[, 2L], , drop = FALSE] * bary$where[ok, 2] +
              mesh$loc[simplex[, 3L], , drop = FALSE] * bary$where[ok, 3])
    if (fm_manifold(mesh, "S2")) {
      loc[ok, ] <- loc[ok, ] / rowSums(loc[ok, ]^2)^0.5 *
        mean(rowSums(mesh$loc^2)^0.5)
    }
  }
  if (format == "sf") {
    loc <- sf::st_as_sf(as.data.frame(loc),
                        coords = seq_len(ncol(loc)),
                        crs = fm_crs(loc))
  }
  loc
}

#' @describeIn fm_bary_loc Extract points on a 1D mesh. Implemented
#' formats are `"numeric"` (default).
#'
#' @export
#' @examples
#' mesh1 <- fm_mesh_1d(1:4)
#' fm_bary_loc(mesh1)
#' (bary1 <- fm_bary(mesh1, seq(0, 5, by = 0.5)))
#' fm_bary_loc(mesh1, bary1)
#' (bary1 <- fm_bary(mesh1, seq(0, 5, by = 0.5), restricted = TRUE))
#' fm_bary_loc(mesh1, bary1)
#' fm_basis(mesh1, bary1)
#' (bary1 <- fm_bary(mesh1, bary1, method = "nearest"))
#' fm_bary_loc(mesh1, bary1)
#' fm_basis(mesh1, bary1)
#' (bary1 <- fm_bary(mesh1, bary1, method = "linear"))
#' fm_bary_loc(mesh1, bary1)
#' fm_basis(mesh1, bary1)
fm_bary_loc.fm_mesh_1d <- function(mesh, bary = NULL, ..., format = NULL) {
  format = match.arg(format, c("numeric"))
  if (is.null(bary)) {
    loc <- mesh$loc
  } else if (NROW(bary) == 0L) {
    loc <- numeric(0L)
  } else {
    loc <- rep(NA_real_, NROW(bary))
    ok <- !is.na(bary$index)
    if (ncol(bary$where) == 1L) {
      loc[ok] <- mesh$loc[bary$index[ok]]
    } else {
      simplex <- fm_bary_simplex(mesh, bary = bary[ok, ])
      loc[ok] <- (mesh$loc[simplex[, 1L]] * bary$where[ok, 1] +
                    mesh$loc[simplex[, 2L]] * bary$where[ok, 2])
    }
  }
  loc
}
