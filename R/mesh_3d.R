#' @title Construct a 3D tetrahedralisation
#'
#' @description
#' Constructs a 3D tetrahedralisation object.
#'
#' @param loc Input coordinates that should be part of the mesh. Can be a
#'   matrix, `sf`, `sfc`, `SpatialPoints`, or other object supported by
#'   [fm_unify_coords()].
#' @param tv Tetrahedron indices, as a N-by-4 index vector into `loc`
#' @param ... Currently unused.
#' @returns An `fm_mesh_3d` object
#' @examples
#' (m <- fm_mesh_3d(
#'   matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0), 4, 3, byrow = TRUE),
#'   matrix(c(1, 2, 3, 4), 1, 4, byrow = TRUE)
#' ))
#' @export
fm_mesh_3d <- function(loc = NULL,
                       tv = NULL,
                       ...) {
  loc <- fm_unify_coords(loc)

  loc.n <- max(0L, nrow(loc))

  stopifnot(all(as.vector(tv) >= 1L))
  stopifnot(all(as.vector(tv) <= NROW(loc)))

  result <- fmesher_mesh3d(
    options = list(),
    loc = loc,
    tv = tv - 1L
  )

  idx_C2R <- function(x) {
    x <- x + 1L
    x[x == 0] <- NA
    x
  }

  m <- structure(
    list(
      meta = list(),
      manifold = result[["manifold"]],
      n = nrow(result[["loc"]]),
      loc = result[["loc"]],
      graph = list(
        tv = idx_C2R(result[["tv"]]),
        vt = lapply(
          result[["vt"]],
          function(xx) {
            yy <- idx_C2R(xx)
            colnames(yy) <- c("t", "vi")
            yy
          }
        ),
        tt = idx_C2R(result[["tt"]]),
        tti = idx_C2R(result[["tti"]]),
        vv = fm_as_dgCMatrix(result[["vv"]]),
        mesh_local = fm_rcdt_2d(
          loc = rbind(
            c(1.0, 0, 0),
            c(0, 1, 0),
            c(0, 0, 1),
            c(0, 0, 0)
          ),
          tv = rbind(
            c(4L, 3L, 2L),
            c(3L, 4L, 1L),
            c(2L, 1L, 4L),
            c(1L, 2L, 3L)
          )
        )
      ),
      idx = list(loc = seq_len(NROW(loc)))
    ),
    class = c("fm_mesh_3d")
  )

  remap_unused <- function(mesh) {
    ## Remap indices to remove unused vertices
    if (length(mesh$graph$vt) == 0) {
      # warning("VT information missing from mesh, rebuilding")
      # Old storage mode: mesh$graph$vt <- rep(NA_integer_, nrow(mesh$loc))
      mesh$graph$vt <- list()
      for (vv in seq_len(nrow(mesh$loc))) {
        mesh$graph$vt[[vv]] <- matrix(NA_integer_, 0, 2)
        colnames(mesh$graph$vt[[vv]]) <- c("t", "vi")
      }
      for (tt in seq_len(nrow(mesh$graph$tv))) {
        for (vvi in seq_len(4)) {
          vv <- mesh$graph$tv[tt, vvi]
          mesh$graph$vt[[vv]] <- rbind(mesh$graph$vt[[vv]], c(tt, vvi))
        }
      }
    }
    used <- vapply(mesh$graph$vt, function(x) NROW(x) > 0, logical(1))
    if (!all(used)) {
      used <- which(used)
      idx.map <- rep(NA, nrow(mesh$loc))
      idx.map[used] <- seq_along(used)
      mesh$loc <- mesh$loc[used, , drop = FALSE]
      mesh$n <- nrow(mesh[["loc"]])
      mesh$graph$tv <-
        matrix(idx.map[as.vector(mesh$graph$tv)], nrow(mesh$graph$tv), 4)
      mesh$graph$vt <- mesh$graph$vt[used]
      ## graph$tt  ## No change needed
      ## graph$tti ## No change needed
      mesh$graph$vv <- mesh$graph$vv[used, used, drop = FALSE]
      if (!is.null(mesh$idx$loc)) {
        mesh$idx$loc <- idx.map[mesh$idx$loc]
      }
    }
    mesh
  }

  m <- remap_unused(m)

  m
}

#' @describeIn fm_mesh_3d Construct a plain Delaunay triangulation in 3D.
#' Requires the `geometry` package.
#' @export
#' @examplesIf requireNamespace("geometry", quietly = TRUE)
#' (m <- fm_delaunay_3d(matrix(rnorm(30), 10, 3)))
fm_delaunay_3d <- function(loc, ...) {
  stopifnot(requireNamespace("geometry"))
  tv <- geometry::delaunayn(loc)
  fm_mesh_3d(
    loc = loc,
    tv = tv,
    ...
  )
}


#' @title Convert objects to `fm_mesh_3d`
#' @describeIn fm_as_mesh_3d Convert an object to `fm_mesh_3d`.
#' @param x Object to be converted.
#' @param ... Arguments passed on to submethods
#' @returns An `fm_mesh_3d` or `fm_mesh_3d_list` object
#' @export
#' @family object creation and conversion
#' @export
#' @examples
#' (m <- fm_mesh_3d(
#'   matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0), 4, 3, byrow = TRUE),
#'   matrix(c(1, 2, 3, 4), 1, 4, byrow = TRUE)
#' ))
#' fm_as_mesh_3d_list(list(m))
fm_as_mesh_3d <- function(x, ...) {
  if (is.null(x)) {
    return(NULL)
  }
  UseMethod("fm_as_mesh_3d")
}
#' @describeIn fm_as_mesh_3d Convert each element of a list
#' @export
fm_as_mesh_3d_list <- function(x, ...) {
  fm_as_list(x, ..., .class_stub = "mesh_3d")
}
#' @rdname fm_as_mesh_3d
#' @param x Object to be converted
#' @export
fm_as_mesh_3d.fm_mesh_3d <- function(x, ...) {
  #  class(x) <- c("fm_mesh_3d", setdiff(class(x), "fm_mesh_3d"))
  x
}


#' @title Convert a 3D mesh to a 3D rgl triangulation
#' @rawNamespace S3method(rgl::as.triangles3d, fm_mesh_3d)
#' @description Extracts a matrix of coordinates of triangles, suitable for
#'   passing to `rgl::triangles3d()`.
#' @param obj An `fm_mesh_3d` object
#' @param subset Character string specifying which triangles to extract. Either
#'   "all" (default) or "boundary".
#' @param \dots Currently unused
#' @returns A 3-column matrix of coordinates of triangles, suitable for
#'   passing to `rgl::triangles3d()`.
#' @examples
#' # Protect against unavailable rgl device by only running interactively
#' if (interactive() &&
#'   requireNamespace("geometry", quietly = TRUE) &&
#'   requireNamespace("rgl", quietly = TRUE)) {
#'   (m <- fm_delaunay_3d(matrix(rnorm(30), 10, 3)))
#'   rgl::open3d()
#'   rgl::triangles3d(rgl::as.triangles3d(m, "boundary"), col = "blue")
#'   rgl::axes3d()
#' }
#'
as.triangles3d.fm_mesh_3d <- function(obj, subset = NULL, ...) {
  subset <- match.arg(subset, c("all", "boundary"))
  tv <- rbind(
    obj$graph$tv[, obj$graph$mesh_local$graph$tv[1, ], drop = FALSE],
    obj$graph$tv[, obj$graph$mesh_local$graph$tv[2, ], drop = FALSE],
    obj$graph$tv[, obj$graph$mesh_local$graph$tv[3, ], drop = FALSE],
    obj$graph$tv[, obj$graph$mesh_local$graph$tv[4, ], drop = FALSE]
  )
  if (identical(subset, "boundary")) {
    keep <- is.na(obj$graph$tt)
    tv <- tv[as.vector(keep), , drop = FALSE]
  }
  loc <- obj$loc[t(tv), ]
  loc
}

#' @describeIn fm_as_mesh_2d Construct a 2D mesh of the boundary of a 3D mesh
#' @export
fm_as_mesh_2d.fm_mesh_3d <- function(x, ...) {
  tv <- rbind(
    x$graph$tv[, x$graph$mesh_local$graph$tv[1, ], drop = FALSE],
    x$graph$tv[, x$graph$mesh_local$graph$tv[2, ], drop = FALSE],
    x$graph$tv[, x$graph$mesh_local$graph$tv[3, ], drop = FALSE],
    x$graph$tv[, x$graph$mesh_local$graph$tv[4, ], drop = FALSE]
  )
  keep <- is.na(x$graph$tt)
  tv <- tv[as.vector(keep), , drop = FALSE]
  fm_rcdt_2d(
    loc = x$loc,
    tv = tv
  )
}
