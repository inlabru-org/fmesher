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
#' \item{vector of knot indices (1D meshes, either for edges or individual
#' knots), or}
#' \item{vector of lower left box indices (2D lattices),}
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
#' @param extra_class character; If non-`NULL` and not already in the class
#'   vector of `bary`, add it to the front of the class vector.
#' @export
fm_bary.fm_bary <- function(bary, ..., extra_class = NULL) {
  if (!is.null(extra_class) && !(extra_class %in% class(bary))) {
    class(bary) <- c(extra_class, class(bary))
  }
  bary
}

#' @describeIn fm_bary Converts a `list` `bary` to `fm_bary`.
#' In the list elements are unnamed, the names `index` and `where` are assumed.
#' @export
fm_bary.list <- function(bary, ..., extra_class = NULL) {
  if (is.null(names(bary))) {
    names(bary) <- c("index", "where")
  }
  bary <- tibble::tibble(
    index = bary[["index"]],
    where = bary[["where"]]
  )
  storage.mode(bary[["index"]]) <- "integer"
  fm_bary(
    structure(
      bary,
      class = c("fm_bary", class(bary))
    ),
    extra_class = extra_class
  )
}

#' @describeIn fm_bary Converts a [tibble::tibble()] `bary` to `fm_bary`
#' @export
fm_bary.tbl_df <- function(bary, ..., extra_class = NULL) {
  stopifnot(
    all(c("index", "where") %in% names(bary))
  )
  storage.mode(bary[["index"]]) <- "integer"
  fm_bary(
    structure(
      bary,
      class = c("fm_bary", class(bary))
    ),
    extra_class = extra_class
  )
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

  idx <- findInterval(loc, knots, all.inside = TRUE)
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
    if (any(ok)) {
      tri[pre_ok_idx[ok]] <- result$index[ok] + 1L
      where_ok <- result$where[ok, , drop = FALSE]
      where_ok <- matrix(pmax(0.0, where_ok), nrow(where_ok), 3)
      where_ok <- where_ok / rowSums(where_ok)
      where[pre_ok_idx[ok], ] <- where_ok
    }
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
      if (any(ok)) {
        tri[subindex[[k]][ok]] <- result$index[ok] + 1L
        where_ok <- result$where[ok, ]
        where_ok <- matrix(pmax(0.0, where_ok), nrow(where_ok), 3)
        where_ok <- where_ok / rowSums(where_ok)
        where[subindex[[k]][ok], ] <- where_ok
      }
    }
  }

  fm_bary(
    tibble::tibble(
      index = tri,
      where = where
    )
  )
}

#' @describeIn fm_bary An `fm_bary` object with columns `index` (vector of
#'   triangle indices) and `where` (4-column matrix of barycentric coordinates).
#'   Points that were not found give `NA` entries in `index` and `where`.
#' @param max_batch_size integer; maximum number of points to process in a
#'   single batch. This speeds up calculations by avoiding repeated large
#'   internal memory allocations and data copies. The default, `NULL`, uses
#'   `max_batch_size = 2e5L`, chosen based on empirical time measurements to
#'   give an approximately optimal runtime.
#'
#' @export
#' @examples
#' m <- fm_mesh_3d(
#'   rbind(
#'     c(1, 0, 0),
#'     c(0, 1, 0),
#'     c(0, 0, 1),
#'     c(0, 0, 0)
#'   ),
#'   matrix(c(1, 2, 3, 4), 1, 4)
#' )
#' b <- fm_bary(m, matrix(c(1, 1, 1) / 4, 1, 3))
fm_bary.fm_mesh_3d <- function(mesh,
                               loc,
                               ...,
                               max_batch_size = NULL) {
  if (inherits(loc, "fm_bary")) {
    return(loc)
  }

  if (is.null(max_batch_size)) {
    max_batch_size <- 2e5L
  }

  pre_ok_idx <-
    which(rowSums(matrix(
      is.na(as.vector(loc)),
      nrow = nrow(loc),
      ncol = ncol(loc)
    )) == 0)
  if (length(pre_ok_idx) <= max_batch_size) {
    result <- fmesher_bary3d(
      mesh_loc = mesh$loc,
      mesh_tv = mesh$graph$tv - 1L,
      loc = loc[pre_ok_idx, , drop = FALSE],
      options = list()
    )
    tet <- rep(NA_integer_, nrow(loc))
    where <- matrix(NA_real_, nrow(loc), 4)
    ok <- result$index >= 0
    tet[pre_ok_idx[ok]] <- result$index[ok] + 1L
    where[pre_ok_idx[ok], ] <- result$where[ok, ]
  } else {
    tet <- rep(NA_integer_, nrow(loc))
    where <- matrix(NA_real_, nrow(loc), 4)
    n_batches <- ceiling(length(pre_ok_idx) / max_batch_size)
    batch_idx <- round(seq(0, length(pre_ok_idx), length.out = n_batches + 1))
    subindex <- split(pre_ok_idx, rep(seq_len(n_batches), diff(batch_idx)))
    for (k in seq_along(subindex)) {
      result <- fmesher_bary3d(
        mesh_loc = mesh$loc,
        mesh_tv = mesh$graph$tv - 1L,
        loc = loc[subindex[[k]], , drop = FALSE],
        options = list()
      )
      ok <- result$index >= 0
      tet[subindex[[k]][ok]] <- result$index[ok] + 1L
      where[subindex[[k]][ok], ] <- result$where[ok, ]
    }
  }

  fm_bary(
    tibble::tibble(
      index = tet,
      where = where
    )
  )
}



#' @describeIn fm_bary An `fm_bary` object with columns `index` (vector of
#'   lattice cell indices) and `where` (4-column matrix of barycentric
#'   coordinates). Points that are outside the lattice are given `NA` entries in
#'   `index` and `where`.
#' @param crs Optional crs information for `loc`
#'
#' @export
#' @examples
#' str(fm_bary(fmexample$mesh, fmexample$loc_sf))
fm_bary.fm_lattice_2d <- function(mesh,
                                  loc,
                                  crs = NULL,
                                  ...) {
  if (inherits(loc, "fm_bary")) {
    if ((nrow(loc) > 0) && (
      min(loc[["index"]]) < 1L ||
        max(loc[["index"]]) > (length(mesh$x) - 1L) * (length(mesh$y) - 1L))) {
      warning("Some 'index' information is outside the lattice.")
    }
    if (ncol(loc[["where"]]) != 4L) {
      stop("Invalid 'where' matrix; should have 4 columns.")
    }
    return(loc)
  }

  loc <- fm_transform(loc, crs = mesh$crs0, crs0 = crs, passthrough = TRUE)
  if (inherits(loc, "sf")) {
    loc <- sf::st_coordinates(loc)[, c("X", "Y"), drop = FALSE]
  } else if (inherits(loc, "Spatial")) {
    stopifnot(fm_safe_sp())
    loc <- sp::coordinates(loc)
  }

  # # Avoid sphere accuracy issues by scaling to unit sphere
  # scale <- 1
  # if (fm_manifold(mesh, "S2")) {
  #   scale <- 1 / mean(rowSums(mesh$loc^2)^0.5)
  #   loc <- loc / rowSums(loc^2)^0.5
  # }

  pre_ok <-
    which(rowSums(matrix(
      is.na(as.vector(loc)),
      nrow = nrow(loc),
      ncol = ncol(loc)
    )) == 0)

  loc <- loc[pre_ok, , drop = FALSE]
  x_idx <- findInterval(loc[, 1L], mesh$x, rightmost.closed = TRUE)
  y_idx <- findInterval(loc[, 2L], mesh$y, rightmost.closed = TRUE)
  ok <- which(x_idx > 0 &
    y_idx > 0 &
    x_idx < length(mesh$x) &
    y_idx < length(mesh$y))
  x_loc <- (loc[ok, 1] - mesh$x[x_idx[ok]]) / diff(mesh$x)[x_idx[ok]]
  y_loc <- (loc[ok, 2] - mesh$y[y_idx[ok]]) / diff(mesh$y)[y_idx[ok]]
  simplex_idx <- x_idx + (y_idx - 1L) * (length(mesh$x) - 1L)

  bary <- fm_bary(
    tibble::tibble(
      index = simplex_idx,
      where = cbind(
        (1 - x_loc) * (1 - y_loc),
        x_loc * (1 - y_loc),
        x_loc * y_loc,
        (1 - x_loc) * y_loc
      )
    )
  )

  bary
}

#' @describeIn fm_bary An `fm_bary` object with columns `index` (vector of
#'   lattice cell indices) and `where` `2^d`-column matrix of barycentric
#'   coordinates). Points that are outside the lattice are given `NA` entries in
#'   `index` and `where`.
#'
#' @export
fm_bary.fm_lattice_Nd <- function(mesh,
                                  loc,
                                  ...) {
  d <- length(mesh$dims)
  d_bary <- 2^d
  if (inherits(loc, "fm_bary")) {
    if ((nrow(loc) > 0) && (
      min(loc[["index"]]) < 1L ||
        max(loc[["index"]]) > prod(mesh$dims - 1L))) {
      warning("Some 'index' information is outside the lattice.")
    }
    if (ncol(loc[["where"]]) != d_bary) {
      stop("Invalid 'where' matrix; should have ", d_bary, " columns.")
    }
    return(loc)
  }

  pre_ok <-
    which(rowSums(matrix(
      is.na(as.vector(loc)),
      nrow = nrow(loc),
      ncol = ncol(loc)
    )) == 0)

  loc <- loc[pre_ok, , drop = FALSE]
  x_idx <- do.call(
    cbind,
    lapply(
      seq_len(d),
      function(k) {
        findInterval(loc[, k],
          mesh$values[[k]],
          rightmost.closed = TRUE
        )
      }
    )
  )
  ok <- rowSums(do.call(
    cbind,
    lapply(
      seq_len(d),
      function(k) {
        x_idx[, k] > 0 &
          x_idx[, k] < mesh$dims[k]
      }
    )
  )) == d
  x_loc <- do.call(
    cbind,
    lapply(
      seq_len(d),
      function(k) {
        (loc[ok, k] - mesh$values[[k]][x_idx[ok, k]]) /
          diff(mesh$values[[k]])[x_idx[ok, k]]
      }
    )
  )
  simplex_idx <- x_idx[ok, 1]
  for (k in seq_len(d - 1) + 1) {
    simplex_idx <- simplex_idx + (x_idx[ok, k] - 1L) *
      prod(mesh$dims[seq_len(k - 1)] - 1L)
  }

  # Vertex order
  # 2d lattice convention: 00 10 11 01
  # Nd lattice convention: 000 100 010 110 001 101 011 111
  # (this gives the same ordering as a local lattice_Nd)
  B <- matrix(0.0, nrow(x_loc), d_bary)
  for (k in seq_len(d_bary)) {
    idx <- (k - 1) %/% 2^(seq_len(d) - 1) %% 2
    xx <- 1.0
    for (j in seq_along(idx)) {
      if (idx[j] == 1) {
        xx <- xx * x_loc[, j]
      } else {
        xx <- xx * (1 - x_loc[, j])
      }
    }
    B[, k] <- xx
  }

  index <- rep(NA_integer_, nrow(loc))
  where <- matrix(NA_real_, nrow(loc), d_bary)
  index[pre_ok[ok]] <- simplex_idx
  where[pre_ok[ok], ] <- B

  bary <- fm_bary(
    tibble::tibble(
      index = index,
      where = where
    )
  )

  bary
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

#' @describeIn fm_bary_simplex Extract the tetrahedron vertex indices for a 3D
#'   mesh
#' @export
#'
#' @examples
#' (m <- fm_mesh_3d(
#'   matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0), 4, 3, byrow = TRUE),
#'   matrix(c(1, 2, 3, 4), 1, 4, byrow = TRUE)
#' ))
#' (bary <- fm_bary(m, rbind(
#'   cbind(0.1, 0.2, 0.3),
#'   cbind(-0.1, 0.2, 0.3)
#' )))
#' fm_bary_simplex(m, bary)
fm_bary_simplex.fm_mesh_3d <- function(mesh, bary = NULL, ...) {
  if (is.null(bary)) {
    return(mesh$graph$tv)
  }
  if (NROW(bary) == 0L) {
    return(matrix(integer(1), 0L, 4L))
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

#' @describeIn fm_bary_simplex Extract the cell vertex indices for a 2D lattice
#' @export
#'
#' @examples
#' m <- fm_lattice_2d(x = 1:3, y = 1:4)
#' bary <- fm_bary(m, cbind(1.5, 3.2))
#' fm_bary_simplex(m, bary)
fm_bary_simplex.fm_lattice_2d <- function(mesh, bary = NULL, ...) {
  simplex <- matrix(0L,
    nrow = (length(mesh$x) - 1L) * (length(mesh$y) - 1L),
    ncol = 4L
  )
  simplex[, 1L] <-
    rep(seq_len(length(mesh$x) - 1L),
      times = length(mesh$y) - 1L
    ) +
    rep((seq_len(length(mesh$y) - 1L) - 1L) * length(mesh$x),
      each = length(mesh$x) - 1L
    )
  simplex[, 2L] <-
    rep(seq_len(length(mesh$x) - 1L) + 1L,
      times = length(mesh$y) - 1L
    ) +
    rep((seq_len(length(mesh$y) - 1L) - 1L) * length(mesh$x),
      each = length(mesh$x) - 1L
    )
  simplex[, 3L] <-
    rep(seq_len(length(mesh$x) - 1L) + 1L,
      times = length(mesh$y) - 1L
    ) +
    rep((seq_len(length(mesh$y) - 1L) - 1L + 1L) * length(mesh$x),
      each = length(mesh$x) - 1L
    )
  simplex[, 4L] <-
    rep(seq_len(length(mesh$x) - 1L),
      times = length(mesh$y) - 1L
    ) +
    rep((seq_len(length(mesh$y) - 1L) - 1L + 1L) * length(mesh$x),
      each = length(mesh$x) - 1L
    )
  if (is.null(bary)) {
    return(simplex)
  }
  if (NROW(bary) == 0L) {
    return(matrix(integer(1), 0L, 4L))
  }
  simplex[bary$index, , drop = FALSE]
}


#' @describeIn fm_bary_simplex Extract the cell vertex indices for a ND lattice
#' @export
#'
#' @examples
#' m <- fm_lattice_Nd(list(x = 1:3, y = 1:4, z = 1:2))
#' (bary <- fm_bary(m, cbind(1.5, 3.2, 1.5)))
#' (fm_bary_simplex(m, bary))
#' fm_bary_loc(m, bary)
fm_bary_simplex.fm_lattice_Nd <- function(mesh, bary = NULL, ...) {
  d <- length(mesh$dims)
  d_bary <- 2^d
  simplex <- matrix(0L,
    nrow = prod(mesh$dims - 1L),
    ncol = d_bary
  )

  simplex <- matrix(0L, prod(mesh$dims - 1L), d_bary)

  # Vertex order
  # 2d lattice convention: 00 10 11 01
  # Nd lattice convention: 000 100 010 110 001 101 011 111
  # (this gives the same ordering as a local lattice_Nd)
  # Local simplex
  local_simplex <- matrix(0L, d_bary, d)
  for (k in seq_len(d_bary)) {
    local_simplex[k, ] <- (k - 1) %/% 2^(seq_len(d) - 1) %% 2
  }
  simplex_root <- matrix(0L, prod(mesh$dims - 1L), ncol = d)
  for (k in seq_len(d)) {
    if (k == 1) {
      simplex_root[, 1] <- rep(seq_len(mesh$dims[1] - 1L),
        times = prod(mesh$dims[-1] - 1L)
      )
    } else {
      simplex_root[, k] <- rep(
        rep(seq_len(mesh$dims[k] - 1L),
          each = prod(mesh$dims[seq_len(k - 1)] - 1L)
        ),
        times = prod(mesh$dims[-seq_len(k)] - 1L)
      )
    }
  }
  for (k in seq_len(d_bary)) {
    simplex[, k] <- simplex_root[, 1] + local_simplex[k, 1]
    for (j in seq_len(d - 1) + 1) {
      simplex[, k] <-
        simplex[, k] + (simplex_root[, j] + local_simplex[k, j] - 1L) *
          prod(mesh$dims[seq_len(j - 1)])
    }
  }

  if (is.null(bary)) {
    return(simplex)
  }
  if (NROW(bary) == 0L) {
    return(matrix(integer(1), 0L, d_bary))
  }
  simplex[bary$index, , drop = FALSE]
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
  format <- match.arg(format, c("matrix", "sf"))
  if (is.null(bary)) {
    loc <- mesh$loc
  } else if (NROW(bary) == 0L) {
    loc <- matrix(0.0, 0L, ncol(mesh$loc))
  } else {
    loc <- matrix(NA_real_, NROW(bary), ncol(mesh$loc))
    ok <- !is.na(bary$index)
    simplex <- fm_bary_simplex(mesh, bary = bary[ok, ])
    loc[ok, ] <- (
      mesh$loc[simplex[, 1L], , drop = FALSE] * bary$where[ok, 1] +
        mesh$loc[simplex[, 2L], , drop = FALSE] * bary$where[ok, 2] +
        mesh$loc[simplex[, 3L], , drop = FALSE] * bary$where[ok, 3]
    )
    if (fm_manifold(mesh, "S2")) {
      loc[ok, ] <- loc[ok, ] / rowSums(loc[ok, ]^2)^0.5 *
        mean(rowSums(mesh$loc^2)^0.5)
    }
  }
  if (format == "sf") {
    loc <- sf::st_as_sf(
      as.data.frame(loc),
      coords = seq_len(ncol(loc)),
      crs = fm_crs(loc)
    )
  }
  loc
}

#' @describeIn fm_bary_loc Extract points on a tetrahedron mesh. Implemented
#' format is `"matrix"` (default).
#' @export
#'
#' @examples
#' (m <- fm_mesh_3d(
#'   matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0), 4, 3, byrow = TRUE),
#'   matrix(c(1, 2, 3, 4), 1, 4, byrow = TRUE)
#' ))
#' (bary <- fm_bary(m, rbind(
#'   cbind(0.1, 0.2, 0.3),
#'   cbind(-0.1, 0.2, 0.3)
#' )))
#' fm_bary_loc(m, bary)
fm_bary_loc.fm_mesh_3d <- function(mesh, bary = NULL, ..., format = NULL) {
  format <- match.arg(format, c("matrix"))
  if (is.null(bary)) {
    loc <- mesh$loc
  } else if (NROW(bary) == 0L) {
    loc <- matrix(0.0, 0L, ncol(mesh$loc))
  } else {
    loc <- matrix(NA_real_, NROW(bary), ncol(mesh$loc))
    ok <- !is.na(bary$index)
    simplex <- fm_bary_simplex(mesh, bary = bary[ok, ])
    loc[ok, ] <- (
      mesh$loc[simplex[, 1L], , drop = FALSE] * bary$where[ok, 1] +
        mesh$loc[simplex[, 2L], , drop = FALSE] * bary$where[ok, 2] +
        mesh$loc[simplex[, 3L], , drop = FALSE] * bary$where[ok, 3] +
        mesh$loc[simplex[, 4L], , drop = FALSE] * bary$where[ok, 4]
    )
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
  format <- match.arg(format, c("numeric"))
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
      loc[ok] <- (
        mesh$loc[simplex[, 1L]] * bary$where[ok, 1] +
          mesh$loc[simplex[, 2L]] * bary$where[ok, 2]
      )
    }
  }
  loc
}

#' @describeIn fm_bary_loc Extract points on a 2D lattice. Implemented
#' formats are `"matrix"` (default) and `"sf"`.
#' @export
#'
#' @examples
#' m <- fm_lattice_2d(x = 1:3, y = 1:4)
#' head(fm_bary_loc(m))
#' (bary <- fm_bary(m, cbind(1.5, 3.2)))
#' fm_bary_loc(m, bary, format = "matrix")
#' fm_bary_loc(m, bary, format = "sf")
fm_bary_loc.fm_lattice_2d <- function(mesh, bary = NULL, ..., format = NULL) {
  format <- match.arg(format, c("matrix", "sf"))
  if (is.null(bary)) {
    loc <- mesh$loc
  } else if (NROW(bary) == 0L) {
    loc <- matrix(0.0, 0L, ncol(mesh$loc))
  } else {
    loc <- matrix(NA_real_, NROW(bary), ncol(mesh$loc))
    ok <- !is.na(bary$index)
    simplex <- fm_bary_simplex(mesh, bary = bary[ok, ])
    loc[ok, ] <- (
      mesh$loc[simplex[, 1L], , drop = FALSE] * bary$where[ok, 1] +
        mesh$loc[simplex[, 2L], , drop = FALSE] * bary$where[ok, 2] +
        mesh$loc[simplex[, 3L], , drop = FALSE] * bary$where[ok, 3] +
        mesh$loc[simplex[, 4L], , drop = FALSE] * bary$where[ok, 4]
    )
    # if (fm_manifold(mesh, "S2")) {
    #   loc[ok, ] <- loc[ok, ] / rowSums(loc[ok, ]^2)^0.5 *
    #     mean(rowSums(mesh$loc^2)^0.5)
    # }
  }
  if (format == "sf") {
    loc <- sf::st_as_sf(
      as.data.frame(loc),
      coords = seq_len(ncol(loc)),
      crs = fm_crs(loc)
    )
  }
  loc
}


#' @describeIn fm_bary_loc Extract points on a ND lattice.
#' @export
#'
#' @examples
#' m <- fm_lattice_Nd(list(x = 1:3, y = 1:4, z = 1:2))
#' head(fm_bary_loc(m))
#' (bary <- fm_bary(m, cbind(1.5, 3.2, 1.5)))
#' fm_bary_loc(m, bary)
fm_bary_loc.fm_lattice_Nd <- function(mesh, bary = NULL, ..., format = NULL) {
  format <- match.arg(format, c("matrix", "sf"))
  stopifnot(format == "matrix")
  if (is.null(bary)) {
    loc <- mesh$loc
  } else if (NROW(bary) == 0L) {
    loc <- matrix(0.0, 0L, ncol(mesh$loc))
  } else {
    loc <- matrix(NA_real_, NROW(bary), ncol(mesh$loc))
    ok <- !is.na(bary$index)
    simplex <- fm_bary_simplex(mesh, bary = bary[ok, , drop = FALSE])
    d <- length(mesh$dims)
    d_bary <- 2^d
    loc[ok, ] <- mesh$loc[simplex[, 1L], , drop = FALSE] * bary$where[ok, 1]
    for (k in seq_len(d_bary - 1) + 1) {
      loc[ok, ] <- (
        loc[ok, ] +
          mesh$loc[simplex[, k], , drop = FALSE] * bary$where[ok, k]
      )
    }
  }
  loc
}
