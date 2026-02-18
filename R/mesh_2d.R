#' @include deprecated.R

# {(v, t, vi)} matrix to v-->{(t,ti)} list of matrices
vt_matrix2list <- function(vt) {
  if (is.null(vt)) {
    return(NULL)
  }
  vt_new <- list()
  for (v in seq_len(NROW(vt))) {
    vt_new[[v]] <- vt[vt[, 1] == v, -1, drop = FALSE]
  }
  vt_new
}

vt_list2matrix <- function(vt) {
  if (is.null(vt)) {
    return(NULL)
  }
  v <- rep(seq_along(vt), vapply(vt, NROW, 1L))
  vt_new <- cbind(v, do.call(rbind, vt))
  vt_new
}

vt_counts <- function(vt) {
  if (is.list(vt)) {
    # list version
    vapply(vt, NROW, 1L)
  } else {
    vt_counts(vt_matrix2list(vt))
    ## Slower than converting to list first:
    # as.vector(Matrix::sparseMatrix(
    #   i = vt[, 1],
    #   j = rep(1L, NROW(vt)),
    #   x = 1L,
    #   dims = c(NROW(vt), 1)
    # ))
  }
}

# fm_mesh_2d ####

#' @title Unify coordinates to 3-column matrix
#'
#' @description Convert coordinate information to a 3-column matrix.
#' This is mainly an internal function, and the interface may change.
#'
#' @param x A object with coordinate information
#' @param crs A optional crs object to convert the coordinates to
#' @returns A coordinate matrix
#' @keywords internal
#' @export
#' @examples
#' fm_unify_coords(fmexample$loc_sf)
#'
fm_unify_coords <- function(x, crs = NULL) {
  UseMethod("fm_unify_coords")
}

#' @rdname fm_unify_coords
#' @usage
#' ## S3 method for class 'NULL'
#' fm_unify_coords(x, crs = NULL)
#' @export
fm_unify_coords.NULL <- function(x, crs = NULL) {
  matrix(0.0, 0, 3)
}

#' @rdname fm_unify_coords
#' @export
fm_unify_coords.default <- function(x, crs = NULL) {
  if (!is.matrix(x)) {
    if (is.vector(x)) {
      x <- matrix(x, 1, length(x))
    } else {
      x <- as.matrix(x)
    }
  }
  if (ncol(x) < 3) {
    if (nrow(x) > 0) {
      while (ncol(x) < 3) {
        x <- cbind(x, 0.0)
      }
    } else {
      x <- matrix(0.0, 0, 3)
    }
  } else if (ncol(x) > 3) {
    stop("Coordinates can have at most 3 columns.")
  }
  colnames(x) <- NULL
  x
}

#' @rdname fm_unify_coords
#' @export
fm_unify_coords.Spatial <- function(x, crs = NULL) {
  fm_safe_sp(force = TRUE)
  x <- fm_transform(
    sp::coordinates(x),
    crs0 = fm_crs(x),
    crs = crs,
    passthrough = TRUE
  )
  if (ncol(x) < 3) {
    if (nrow(x) > 0) {
      while (ncol(x) < 3) {
        x <- cbind(x, 0.0)
      }
    } else {
      x <- matrix(0.0, 0, 3)
    }
  } else if (ncol(x) > 3) {
    stop("Coordinates can have at most 3 columns.")
  }
  colnames(x) <- NULL
  x
}

#' @rdname fm_unify_coords
#' @export
fm_unify_coords.sf <- function(x, crs = NULL) {
  fm_unify_coords.sfc(sf::st_geometry(x), crs = crs)
}

#' @rdname fm_unify_coords
#' @export
fm_unify_coords.sfc <- function(x, crs = NULL) {
  loc <- fm_zm(x)
  loc <- sf::st_coordinates(loc)
  loc <- loc[, intersect(colnames(loc), c("X", "Y", "Z")), drop = FALSE]
  x <- fm_transform(
    loc,
    crs0 = fm_crs(x),
    crs = crs,
    passthrough = TRUE
  )
  if (ncol(x) < 3) {
    if (nrow(x) > 0) {
      while (ncol(x) < 3) {
        x <- cbind(x, 0.0)
      }
    } else {
      x <- matrix(0.0, 0, 3)
    }
  } else if (ncol(x) > 3) {
    stop("Coordinates can have at most 3 columns.")
  }
  colnames(x) <- NULL
  x
}


unify_segm_coords <- function(segm, crs = NULL) {
  if (is.null(segm)) {
    return(NULL)
  }
  unify.one.segm <- function(segm, crs = NULL) {
    segm <- fm_transform(segm, crs, passthrough = TRUE)
    if (ncol(segm$loc) == 2) {
      segm$loc <- cbind(segm$loc, 0.0)
    }
    segm
  }
  if (inherits(segm, "fm_segm")) {
    segm <- unify.one.segm(segm, crs = crs)
  } else {
    for (j in seq_along(segm)) {
      if (!is.null(segm[[j]])) {
        segm[[j]] <- unify.one.segm(segm[[j]], crs = crs)
      }
    }
  }
  segm
}


handle_rcdt_options_inla <- function(
  ...,
  quality.spec = NULL,
  cutoff = 1e-12,
  extend = NULL,
  refine = NULL,
  delaunay = TRUE,
  .n,
  .loc
) {
  options <- list(cutoff = cutoff, delaunay = delaunay)
  if (is.null(quality.spec)) {
    quality <- NULL
  } else {
    quality <- rep(NA, .n$segm + .n$lattice + .n$loc)
    ## Order must be same as in .loc
    if (!is.null(quality.spec$segm)) {
      quality[seq_len(.n$segm)] <- quality.spec$segm
    }
    if (!is.null(quality.spec$lattice)) {
      quality[.n$segm + seq_len(.n$lattice)] <- quality.spec$lattice
    }
    if (!is.null(quality.spec$loc)) {
      quality[.n$segm + .n$lattice + seq_len(.n$loc)] <- quality.spec$loc
    }
    ## NA:s will be replaced with max.edge settings below.

    options <- c(options, list(quality = quality))
  }

  cet_sides <- NULL
  cet_margin <- NULL
  if (isTRUE(extend)) {
    extend <- list()
  }
  if (inherits(extend, "list")) {
    cet_sides <- ifelse(is.null(extend$n), 16L, as.integer(extend$n))
    cet_margin <- ifelse(is.null(extend$offset), -0.1, extend$offset)
  }
  options <- c(options, list(cet_sides = cet_sides, cet_margin = cet_margin))

  if (isTRUE(refine)) {
    refine <- list()
  }
  if (inherits(refine, "list")) {
    # Override possible delaunay=FALSE option
    options[["delaunay"]] <- TRUE

    rcdt_min_angle <- 0
    rcdt_max_edge <- 0
    # Multiply by 2 to cover S2; could remove if supplied manifold info
    max.edge.default <- fm_diameter(.loc) * 2
    if ((inherits(extend, "list")) && (!is.null(extend$offset))) {
      max.edge.default <- (max.edge.default +
        max(0, 2 * extend$offset))
      max.edge.default <- (max.edge.default *
        (1 + max(0, -2 * extend$offset)))
    }
    rcdt_min_angle <- ifelse(is.null(refine$min.angle), 21, refine$min.angle)
    rcdt_max_edge <- ifelse(
      is.null(refine$max.edge) ||
        is.na(refine$max.edge),
      max.edge.default,
      refine$max.edge
    )
    max_edge_extra <- ifelse(
      is.null(refine$max.edge.extra) ||
        is.na(refine$max.edge.extra),
      rcdt_max_edge,
      refine$max.edge.extra
    )

    if (!is.null(refine[["max.n.strict"]]) &&
      !is.na(refine$max.n.strict)) {
      rcdt_max_n0 <- as.integer(refine$max.n.strict)
    } else {
      rcdt_max_n0 <- -1L
    }
    if (!is.null(refine[["max.n"]]) &&
      !is.na(refine$max.n)) {
      rcdt_max_n1 <- as.integer(refine$max.n)
    } else {
      rcdt_max_n1 <- -1L
    }

    if (!is.null(options[["quality"]])) {
      options[["quality"]][is.na(options[["quality"]])] <- max_edge_extra
    }

    options <-
      c(
        options,
        list(
          rcdt_min_angle = rcdt_min_angle,
          rcdt_max_edge = rcdt_max_edge,
          rcdt_max_n0 = rcdt_max_n0,
          rcdt_max_n1 = rcdt_max_n1
        )
      )
  }

  options
}


#' @title Refined Constrained Delaunay Triangulation
#'
#' @description
#' Computes a refined constrained Delaunay triangulation on R2 or S2.
#'
#' @param loc Input coordinates that should be part of the mesh. Can be a
#'   matrix, `sf`, `sfc`, `SpatialPoints`, or other object supported by
#'   [fm_unify_coords()].
#' @param tv Initial triangulation, as a N-by-3 index vector into `loc`
#' @param boundary,interior Objects supported by [fm_as_segm()].
#' If `boundary` is `numeric`, `fm_nonconvex_hull(loc, convex = boundary)` is
#' used.
#' @param extend `logical` or `list` specifying whether to extend the
#' data region, with parameters \describe{ \item{list("n")}{the number of edges
#' in the extended boundary (default=16)} \item{list("offset")}{the extension
#' distance.  If negative, interpreted as a factor relative to the approximate
#' data diameter (default=-0.10)} } Setting to `FALSE` is only useful in
#' combination `lattice` or `boundary`.
#' @param refine `logical` or `list` specifying whether to refine the
#' triangulation, with parameters \describe{ \item{list("min.angle")}{the
#' minimum allowed interior angle in any triangle.  The algorithm is guaranteed
#' to converge for `min.angle` at most 21 (default=`21`)}
#' \item{list("max.edge")}{the maximum allowed edge length in any triangle.  If
#' negative, interpreted as a relative factor in an ad hoc formula depending on
#' the data density (default=`Inf`)} \item{list("max.n.strict")}{the
#' maximum number of vertices allowed, overriding `min.angle` and
#' `max.edge` (default=-1, meaning no limit)} \item{list("max.n")}{the
#' maximum number of vertices allowed, overriding `max.edge` only
#' (default=-1, meaning no limit)} }
#' @param lattice An `fm_lattice_2d` object, generated by
#' [fm_lattice_2d()], specifying points on a regular lattice.
#' @param cutoff The minimum allowed distance between points.  Point at most as
#' far apart as this are replaced by a single vertex prior to the mesh
#' refinement step.
#' @param globe If non-NULL, an integer specifying the level of subdivision
#' for global mesh points, used with [fmesher_globe_points()]
#' @param quality.spec List of vectors of per vertex `max.edge` target
#' specification for each location in `loc`, `boundary/interior`
#' (`segm`), and `lattice`.  Only used if refining the mesh.
#' @param crs Optional crs object
#' @param delaunay logical; If `FALSE`, `refine` is `FALSE`, and a ready-made
#'   mesh is provided, only creates the mesh data structure. Default `TRUE`, for
#'   ensuring a Delaunay triangulation.
#' @param ... Currently passed on to `fm_mesh_2d_inla` or converted to
#' [fmesher_rcdt()] options.
#' @returns An `fm_mesh_2d` object
#' @examples
#' (m <- fm_rcdt_2d_inla(
#'   boundary = fm_nonconvex_hull(cbind(0, 0), convex = 5)
#' ))
#'
#' @export
fm_rcdt_2d <-
  function(...) {
    fm_rcdt_2d_inla(...)
  }

#' @describeIn fm_rcdt_2d Legacy method for the `INLA::inla.mesh.create()`
#' interface
#' @inheritSection fm_mesh_2d INLA compatibility
#' @export
fm_rcdt_2d_inla <- function(loc = NULL,
                            tv = NULL,
                            boundary = NULL,
                            interior = NULL,
                            extend = (missing(tv) || is.null(tv)),
                            refine = FALSE,
                            lattice = NULL,
                            globe = NULL,
                            cutoff = 1e-12,
                            quality.spec = NULL,
                            crs = NULL,
                            delaunay = TRUE,
                            ...) {
  crs.target <- crs
  if (!fm_crs_is_null(crs) &&
    fm_crs_is_geocent(crs)) {
    ## Build all geocentric meshes on a sphere, and transform afterwards,
    ## to allow general geoids.
    crs <- fm_crs("sphere")
  }

  if (!is.null(loc) && !is.matrix(loc)) {
    crs.loc <- fm_crs(loc)
  } else {
    crs.loc <- NULL
  }
  loc <- fm_unify_coords(loc)
  if (!fm_crs_is_null(crs.loc) && !fm_crs_is_null(crs)) {
    loc <- fm_transform(loc, crs = crs, passthrough = TRUE, crs0 = crs.loc)
    loc <- fm_unify_coords(loc)
  }

  if (!is.null(globe)) {
    loc.globe <- fmesher_globe_points(globe = globe)
    crs.globe <- fm_crs("sphere")
    if (!fm_crs_is_null(crs.globe) && !fm_crs_is_null(crs)) {
      loc.globe <- fm_transform(loc.globe,
        crs = crs,
        passthrough = TRUE,
        crs0 = crs.globe
      )
      loc.globe <- fm_unify_coords(loc.globe)
    }
    loc <- rbind(loc, loc.globe)
  }
  loc.n <- max(0L, nrow(loc))

  lattice.boundary <- NULL
  if (is.null(lattice) || !is.null(tv)) {
    if (!is.null(lattice)) {
      warning("Both 'lattice' and 'tv' specified.  Ignoring 'lattice'.")
    }
    lattice <- list(loc = NULL, segm = NULL)
    lattice.n <- 0L
  } else {
    lattice <- fm_as_lattice_2d(lattice)

    if (!fm_crs_is_null(fm_crs(lattice))) {
      lattice <- fm_transform(
        lattice,
        crs = crs,
        passthrough = TRUE
      )
    }
    if (NCOL(lattice$loc) == 2) {
      lattice$loc <- cbind(lattice$loc, 0.0)
    }
    if (NCOL(lattice$segm$loc) == 2) {
      lattice$segm$loc <- cbind(lattice$segm$loc, 0.0)
    }

    if (is.logical(extend) && !extend) {
      lattice.boundary <- lattice$segm
    }
  }
  lattice.n <- max(0L, nrow(lattice$loc))

  segm.n <- 0L
  if (!is.null(lattice.boundary) && is.null(boundary)) {
    boundary <- lattice.boundary
    lattice.boundary <- NULL
  }
  if (is.null(boundary)) {
    bnd <- NULL
    bnd_grp <- NULL
    loc.bnd <- matrix(0.0, 0, 3)
  } else {
    if (is.numeric(boundary)) {
      boundary <- fm_nonconvex_hull(loc, convex = boundary)
    } else {
      boundary <- fm_as_segm(boundary)
    }
    if (is.null(boundary$loc)) {
      boundary$loc <- loc
      boundary$crs <- fm_crs(crs)
    } else if (!fm_crs_is_null(crs)) {
      boundary <- fm_transform(boundary, crs = crs, passthrough = TRUE)
    }
    if (!is.null(lattice.boundary)) {
      boundary <-
        fm_segm_join(fm_as_segm_list(list(boundary, lattice.boundary)))
    }

    bnd <- segm.n + boundary$idx
    bnd_grp <- boundary$grp
    if (ncol(boundary$loc) == 2) {
      boundary$loc <- cbind(boundary$loc, 0.0)
    }
    segm.n <- segm.n + max(0L, nrow(boundary$loc))
    loc.bnd <- boundary$loc
  }

  if (is.null(interior)) {
    int <- NULL
    int_grp <- NULL
    loc.int <- matrix(0.0, 0, 3)
  } else {
    interior <- fm_as_segm(interior)
    if (is.null(interior$loc)) {
      interior$loc <- loc
      interior$crs <- fm_crs(crs)
    } else if (!fm_crs_is_null(crs)) {
      interior <- fm_transform(interior, crs = crs, passthrough = TRUE)
    }
    int <- segm.n + interior$idx
    int_grp <- interior$grp
    if (ncol(interior$loc) == 2) {
      interior$loc <- cbind(interior$loc, 0.0)
    }
    segm.n <- segm.n + max(0L, nrow(interior$loc))
    loc.int <- interior$loc
  }

  if (!is.null(tv)) {
    stopifnot(all(as.vector(tv) >= 1L))
    stopifnot(all(as.vector(tv) <= NROW(loc)))
  }
  loc <- rbind(loc.bnd, loc.int, lattice$loc, loc)

  options <- handle_rcdt_options_inla(
    extend = extend,
    refine = refine,
    cutoff = cutoff,
    quality.spec = quality.spec,
    delaunay = delaunay,
    ...,
    .n = list(
      segm = segm.n,
      lattice = lattice.n,
      loc = loc.n
    ),
    .loc = loc
  )

  if (!is.null(tv)) {
    tv <- tv + segm.n + lattice.n - 1L
  }
  if (!is.null(bnd)) {
    bnd <- bnd - 1L
  }
  if (!is.null(int)) {
    int <- int - 1L
  }
  result <- fmesher_rcdt(
    options = options,
    loc = loc, tv = tv,
    boundary = bnd, interior = int,
    boundary_grp = bnd_grp, interior_grp = int_grp
  )

  idx_C2R <- function(x) {
    x <- x + 1L
    x[x == 0] <- NA
    x
  }

  if (!fm_crs_is_null(crs) &&
    !fm_crs_is_identical(crs, crs.target)) {
    ## Target is a non-spherical geoid
    result[["s"]] <- fm_transform(result[["s"]], crs0 = crs, crs = crs.target)
    crs <- crs.target
  }

  split_idx <- function(idx, splits) {
    cumulative_splits <- c(0L, cumsum(splits))
    idx <- lapply(
      seq_along(splits),
      function(k) {
        if (splits[k] > 0) {
          idx[cumulative_splits[k] + seq_len(splits[k])]
        } else {
          NULL
        }
      }
    )
    names(idx) <- names(splits)
    idx
  }
  idx.all <- idx_C2R(result[["idx"]])
  idx <- split_idx(idx.all, c(segm = segm.n, lattice = lattice.n, loc = loc.n))

  m <- structure(
    list(
      meta = list(
        is.refined = !is.null(options[["rcdt_max_edge"]])
      ),
      manifold = result[["manifold"]],
      n = nrow(result[["s"]]),
      loc = result[["s"]],
      graph = list(
        tv = idx_C2R(result[["tv"]]),
        # Note: triangle vt indexing will be sorted out in remap_unused.
        vt = lapply(result[["vt"]], idx_C2R),
        tt = idx_C2R(result[["tt"]]),
        tti = idx_C2R(result[["tti"]]),
        vv = fm_as_dgCMatrix(result[["vv"]])
      ),
      segm = list(
        int = fm_segm(
          idx = idx_C2R(result[["segm.int.idx"]]),
          grp = result[["segm.int.grp"]],
          is.bnd = FALSE
        ),
        bnd = fm_segm(
          idx = idx_C2R(result[["segm.bnd.idx"]]),
          grp = result[["segm.bnd.grp"]],
          is.bnd = TRUE
        )
      ),
      idx = idx,
      crs = fm_crs(crs)
    ),
    class = c("fm_mesh_2d", "inla.mesh")
  )

  remap_unused <- function(mesh) {
    ## Remap indices to remove unused vertices
    if (length(mesh$graph$vt) > 0) {
      for (vv in seq_len(nrow(mesh$loc))) {
        vt <- mesh$graph$vt[[vv]]
        # Need to do the C->R index conversion for the triangle indices here!
        mesh$graph$vt[[vv]] <-
          matrix(c(as.integer(names(vt)) + 1L, vt), length(vt), 2,
                 dimnames = list(NULL, c("t", "vi")))
      }
    } else {
      # warning("VT information missing from mesh, rebuilding")
      # Old storage mode: mesh$graph$vt <- rep(NA_integer_, nrow(mesh$loc))
      mesh$graph$vt <- list()
      for (vv in seq_len(nrow(mesh$loc))) {
        mesh$graph$vt[[vv]] <- matrix(NA_integer_, 0, 2,
                                      dimnames = list(NULL, c("t", "vi")))
      }
      for (tt in seq_len(nrow(mesh$graph$tv))) {
        for (vvi in seq_len(3)) {
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
        matrix(idx.map[as.vector(mesh$graph$tv)], nrow(mesh$graph$tv), 3)
      mesh$graph$vt <- mesh$graph$vt[used]
      ## graph$tt  ## No change needed
      ## graph$tti ## No change needed
      mesh$graph$vv <- mesh$graph$vv[used, used, drop = FALSE]
      if (!is.null(mesh$idx$loc)) {
        mesh$idx$loc <- idx.map[mesh$idx$loc]
      }
      if (!is.null(mesh$idx$lattice)) {
        mesh$idx$lattice <- idx.map[mesh$idx$lattice]
      }
      if (!is.null(mesh$idx$segm)) {
        mesh$idx$segm <- idx.map[mesh$idx$segm]
      }
      mesh$segm$bnd$idx <-
        matrix(idx.map[mesh$segm$bnd$idx], nrow(mesh$segm$bnd$idx), 2)
      mesh$segm$int$idx <-
        matrix(idx.map[mesh$segm$int$idx], nrow(mesh$segm$int$idx), 2)
      mesh$segm$bnd$idx[mesh$segm$bnd$idx == 0L] <- NA
      mesh$segm$int$idx[mesh$segm$int$idx == 0L] <- NA
    }
    mesh
  }

  # Note: this also handles the C->R conversion for triangle indexing in vt.
  m <- remap_unused(m)

  m
}

#' @describeIn fm_rcdt_2d Construct a plain Delaunay triangulation.
#' @export
#' @examples
#' fm_delaunay_2d(matrix(rnorm(30), 15, 2))
#'
fm_delaunay_2d <- function(loc, crs = NULL, ...) {
  if (is.null(crs) && !is.matrix(loc)) {
    crs <- fm_crs(loc)
  }
  loc <- fm_unify_coords(loc, crs = crs)

  hull <- grDevices::chull(loc[, 1], loc[, 2])
  bnd <- fm_segm(
    loc = loc[hull[rev(seq_along(hull))], , drop = FALSE],
    is.bnd = TRUE
  )
  mesh <- fm_rcdt_2d_inla(
    loc = loc,
    boundary = bnd,
    extend = list(n = 3),
    refine = FALSE,
    crs = crs,
    ...
  )
  mesh
}


#' @title Make a 2D mesh object
#' @export
#' @param ... Currently passed on to `fm_mesh_2d_inla`
#' @family object creation and conversion
#' @section INLA compatibility:
#' For mesh and curve creation, the [fm_rcdt_2d_inla()], [fm_mesh_2d_inla()],
#' and [fm_nonconvex_hull_inla()] methods will keep the interface syntax used by
#' `INLA::inla.mesh.create()`, `INLA::inla.mesh.2d()`, and
#' `INLA::inla.nonconvex.hull()` functions, respectively, whereas the
#' [fm_rcdt_2d()], [fm_mesh_2d()], and [fm_nonconvex_hull()] interfaces may be
#' different, and potentially change in the future.
#'
#' @examples
#' fm_mesh_2d_inla(boundary = fm_extensions(cbind(2, 1), convex = 1, 2))
#'
fm_mesh_2d <- function(...) {
  fm_mesh_2d_inla(...)
}

#' @describeIn fm_mesh_2d Legacy method for `INLA::inla.mesh.2d()`
#' Create a triangle mesh based on initial point locations, specified or
#' automatic boundaries, and mesh quality parameters.
#' @export
#'
#' @param loc Matrix of point locations to be used as initial triangulation
#' nodes.  Can alternatively be a `sf`, `sfc`, `SpatialPoints` or
#' `SpatialPointsDataFrame` object.
#' @param loc.domain Matrix of point locations used to determine the domain
#' extent.  Can alternatively be a `SpatialPoints` or
#' `SpatialPointsDataFrame` object.
#' @param offset The automatic extension distance.  One or two values, for an
#' inner and an optional outer extension.  If negative, interpreted as a factor
#' relative to the approximate data diameter (default=-0.10???)
#' @param n The number of initial nodes in the automatic extensions
#' (default=16)
#' @param boundary one or more (as list) of [fm_segm()] objects, or objects
#' supported by [fm_as_segm()]
#' @param interior one object supported by [fm_as_segm()], or (from version
#' `0.2.0.9016`) a list of such objects.  If a list, the objects are joined
#' into a single object.
#' @param max.edge The largest allowed triangle edge length.  One or two
#' values.
#' @param min.angle The smallest allowed triangle angle.  One or two values.
#' (Default=21)
#' @param cutoff The minimum allowed distance between points.  Point at most as
#' far apart as this are replaced by a single vertex prior to the mesh
#' refinement step.
#' @param max.n.strict The maximum number of vertices allowed, overriding
#' `min.angle` and `max.edge` (default=-1, meaning no limit).  One or
#' two values, where the second value gives the number of additional vertices
#' allowed for the extension.
#' @param max.n The maximum number of vertices allowed, overriding
#' `max.edge` only (default=-1, meaning no limit).  One or two values,
#' where the second value gives the number of additional vertices allowed for
#' the extension.
# @param plot.delay On Linux (and Mac if appropriate X11 libraries are
# installed), specifying a nonnegative numeric value activates a rudimentary
# plotting system in the underlying `fmesher` program, showing the
# triangulation algorithm at work, with waiting time factor `plot.delay`
# between each step.
#' @param plot.delay If logical `TRUE` or a negative numeric value,
#' activates displaying the
#' result after each step of the multi-step domain extension algorithm.
#' @param crs An optional [fm_crs()], `sf::crs` or `sp::CRS` object
#' @returns An `fm_mesh_2d` object.
#' @author Finn Lindgren <Finn.Lindgren@@gmail.com>
#' @seealso [fm_rcdt_2d()], [fm_mesh_2d()], [fm_delaunay_2d()],
#' [fm_nonconvex_hull()], [fm_extensions()], [fm_refine()]
fm_mesh_2d_inla <- function(loc = NULL,
                            loc.domain = NULL,
                            offset = NULL,
                            n = NULL,
                            boundary = NULL,
                            interior = NULL,
                            max.edge = NULL,
                            min.angle = NULL,
                            cutoff = 1e-12,
                            max.n.strict = NULL,
                            max.n = NULL,
                            plot.delay = NULL,
                            crs = NULL,
                            ...) {
  ## plot.delay: Do plotting.
  ## NULL --> No plotting
  ## <0  --> Intermediate meshes displayed at the end
  ## TRUE  --> Intermediate meshes displayed at the end
  ## >0   --> Dynamical fmesher X11 plotting is not available in the R interface
  if (is.null(plot.delay)) {
    plot.intermediate <- FALSE
  } else if (is.logical(plot.delay)) {
    plot.intermediate <- plot.delay
  } else {
    plot.intermediate <- plot.delay < 0
  }

  if ((missing(max.edge) || is.null(max.edge)) &&
    (missing(max.n.strict) || is.null(max.n.strict)) &&
    (missing(max.n) || is.null(max.n))) {
    max.edge <- NA
  }

  if (!is.null(crs)) {
    issphere <- fm_crs_is_identical(crs, fm_crs("sphere"))
    isgeocentric <- fm_crs_is_geocent(crs)
    if (isgeocentric) {
      crs.target <- crs
      crs <- fm_CRS("sphere")
    }
  }

  loc <- fm_unify_coords(loc, crs = crs)
  loc.domain <- fm_unify_coords(loc.domain, crs = crs)

  boundary <- fm_as_segm_list(boundary)
  if (is.list(interior)) {
    interior <- fm_segm_join(fm_as_segm_list(interior))
  } else {
    interior <- fm_as_segm(interior)
  }

  if (length(boundary) == 0) {
    list(NULL)
  }

  if (missing(offset) || is.null(offset)) {
    if (length(boundary) < 2) {
      offset <- -0.05
    } else {
      offset <- c(-0.05, -0.15)
    }
  }
  if (missing(n) || is.null(n)) {
    n <- c(8)
  }
  if (missing(max.edge) || is.null(max.edge)) {
    max.edge <- c(NA)
  }
  if (missing(min.angle) || is.null(min.angle)) {
    min.angle <- c(21)
  }
  if (missing(max.n.strict) || is.null(max.n.strict)) {
    max.n.strict <- c(NA)
  }
  if (missing(max.n) || is.null(max.n)) {
    max.n <- c(NA)
  }
  if (missing(cutoff) || is.null(cutoff)) {
    cutoff <- 1e-12
  }

  num.layers <-
    max(c(
      length(boundary), length(offset), length(n),
      length(min.angle), length(max.edge),
      length(max.n.strict), length(max.n)
    ))
  if (num.layers > 2) {
    warning(paste("num.layers=", num.layers, " > 2 detected.  ",
      "Excess information ignored.",
      sep = ""
    ))
    num.layers <- 2
  }

  if (length(boundary) < num.layers) {
    boundary <- c(boundary, rep(list(NULL), num.layers - length(boundary)))
  }
  if (length(min.angle) < num.layers) {
    min.angle <- c(min.angle, min.angle)
  }
  if (length(max.n.strict) < num.layers) {
    max.n0 <- c(max.n.strict, max.n.strict)
  }
  if (length(max.n) < num.layers) {
    max.n <- c(max.n, max.n)
  }
  if (length(max.edge) < num.layers) {
    max.edge <- c(max.edge, max.edge)
  }
  if (length(offset) < num.layers) {
    offset <- c(offset, -0.15)
  }
  if (length(n) < num.layers) {
    n <- c(n, 16)
  }
  if (length(n) < num.layers) {
    n <- c(n, 16)
  }

  if (fm_diameter(loc) +
    fm_diameter(loc.domain) +
    fm_diameter(interior) == 0.0) {
    bnd_diam <- 0.0
    for (k in seq_len(num.layers)) {
      if ((length(boundary) >= k) && !is.null(boundary[[k]])) {
        bnd_diam <- bnd_diam + fm_diameter(boundary[[k]])
      }
      if (offset[k] < 0) {
        if (bnd_diam == 0.0) {
          offset[k] <- 1
        }
      }
    }
  }

  ## Unify the dimensionality of the boundary&interior segments input
  ## and optionally transform coordinates.
  boundary <- unify_segm_coords(boundary, crs = crs)
  interior <- unify_segm_coords(interior, crs = crs)

  if (is.null(boundary[[1]])) {
    ## Triangulate to get inner domain boundary
    ## Constraints included only to get proper domain extent
    ## First, attach the loc points to the domain definition set
    if (!is.null(loc) && !is.null(loc.domain)) {
      loc.domain <- rbind(loc.domain, loc)
    }
    mesh1 <-
      fm_rcdt_2d(
        loc = loc.domain,
        boundary = boundary[[1]],
        interior = interior,
        cutoff = cutoff,
        extend = list(n = n[1], offset = offset[1]),
        refine = FALSE,
        crs = crs
      )

    ## Save the resulting boundary
    boundary1 <- fm_segm(mesh1, boundary = TRUE)
    interior1 <- fm_segm(mesh1, boundary = FALSE)

    if (plot.intermediate) {
      plot(mesh1)
    }

    ## Triangulate inner domain
    mesh2 <-
      fm_rcdt_2d(
        loc = loc,
        boundary = boundary1,
        interior = interior1,
        cutoff = cutoff,
        extend = if (fm_manifold(mesh1, "S2")) {
          list(n = n[1], offset = offset[1])
        } else {
          FALSE ## Should have no effect
        },
        refine =
          list(
            min.angle = min.angle[1],
            max.edge = max.edge[1],
            max.edge.extra = max.edge[1],
            max.n.strict = max.n.strict[1],
            max.n = max.n[1]
          ),
        crs = crs
      )
  } else {
    mesh2 <-
      fm_rcdt_2d(
        loc = loc,
        boundary = boundary[[1]],
        interior = interior,
        cutoff = cutoff,
        extend = FALSE, # Should have no effect
        refine =
          list(
            min.angle = min.angle[1],
            max.edge = max.edge[1],
            max.edge.extra = max.edge[1],
            max.n.strict = max.n.strict[1],
            max.n = max.n[1]
          ),
        crs = crs
      )
  }

  boundary2 <- fm_segm(mesh2, boundary = TRUE)
  interior2 <- fm_segm(mesh2, boundary = FALSE)

  if (plot.intermediate) {
    plot(mesh2)
  }

  if (num.layers == 1) {
    if (!is.null(crs) && isgeocentric && !issphere) {
      mesh2$loc <- fm_transform(mesh2$loc, crs0 = mesh2$crs, crs = crs.target)
      mesh2$crs <- crs.target
    }
    return(mesh2)
  }

  ## Triangulate inner+outer domain
  mesh3 <-
    fm_rcdt_2d(
      loc = rbind(loc, mesh2$loc),
      boundary = boundary[[2]],
      interior = fm_segm(boundary2, interior2, is.bnd = FALSE),
      cutoff = cutoff,
      extend = list(n = n[2], offset = offset[2]),
      refine =
        list(
          min.angle = min.angle[2],
          max.edge = max.edge[2],
          max.edge.extra = max.edge[2],
          max.n.strict = mesh2$n + max.n.strict[2],
          max.n = mesh2$n + max.n[2]
        ),
      crs = crs
    )

  ## Hide generated points, to match regular fm_rcdt_2d_inla output
  mesh3$idx$loc <- mesh3$idx$loc[seq_len(nrow(loc))]

  ## Obtain the corresponding segm indices.
  segm.loc <- matrix(0.0, 0, 3)
  for (k in seq_along(boundary)) {
    if (!is.null(boundary[[k]])) {
      segm.loc <- rbind(segm.loc, boundary[[k]]$loc)
    }
  }
  if (!is.null(interior)) {
    segm.loc <- rbind(segm.loc, interior$loc)
  }
  if (nrow(segm.loc) > 0) {
    proj <- fm_basis(mesh3, loc = segm.loc, full = TRUE)
    mesh3$idx$segm <- rep(NA, nrow(segm.loc))
    if (any(proj$ok)) {
      t.idx <- proj$bary$index[proj$ok]
      tv.idx <- max.col(proj$bary$where[proj$ok, , drop = FALSE],
        ties.method = "first"
      )
      mesh3$idx$segm[proj$ok] <-
        mesh3$graph$tv[t.idx + nrow(mesh3$graph$tv) * (tv.idx - 1)]
    }
  } else {
    mesh3$idx$segm <- NULL
  }

  if (!is.null(crs) && isgeocentric && !issphere) {
    mesh3$loc <-
      fm_transform(mesh3$loc, crs0 = mesh3$crs, crs = crs.target)
    mesh3$crs <- crs.target
  }

  if (plot.intermediate) {
    plot(mesh3)
  }

  mesh3
}

#' @title Convert objects to `fm_mesh_2d`
#' @describeIn fm_as_mesh_2d Convert an object to `fm_mesh_2d`.
#' @param x Object to be converted.
#' @param ... Arguments passed on to submethods
#' @returns An `fm_mesh_2d` or `fm_mesh_2d_list` object
#' @export
#' @family object creation and conversion
#' @export
#' @examples
#' fm_as_mesh_2d_list(list(fm_mesh_2d(cbind(2, 1))))
fm_as_mesh_2d <- function(x, ...) {
  if (is.null(x)) {
    return(NULL)
  }
  UseMethod("fm_as_mesh_2d")
}
#' @describeIn fm_as_mesh_2d Convert each element of a list
#' @export
fm_as_mesh_2d_list <- function(x, ...) {
  fm_as_list(x, ..., .class_stub = "mesh_2d")
}
#' @rdname fm_as_mesh_2d
#' @param x Object to be converted
#' @export
fm_as_mesh_2d.fm_mesh_2d <- function(x, ...) {
  #  class(x) <- c("fm_mesh_2d", setdiff(class(x), "fm_mesh_2d"))
  x
}
#' @rdname fm_as_mesh_2d
#' @export
#' @method fm_as_mesh_2d inla.mesh
fm_as_mesh_2d.inla.mesh <- function(x, ...) {
  x[["crs"]] <- fm_crs(x[["crs"]])
  if (!is.null(x$segm$bnd)) {
    x$segm$bnd <- fm_as_fm(x$segm$bnd)
  }
  if (!is.null(x$segm$int)) {
    x$segm$int <- fm_as_fm(x$segm$int)
  }
  class(x) <- c("fm_mesh_2d", class(x))
  x
}


# Hex points ####

#' @title Create hexagon lattice points
#' @description Create hexagon lattice points within a boundary
#' @param bnd Boundary object
#' @param x_bin Number of bins in x axis
#' @param edge_len_n Number of edge length of mesh from the boundary to create
#'   hexagon mesh using x_bin
#' @return A list with lattice points, edge length, and inner boundary
#' @author Man Ho Suen <M.H.Suen@@sms.ed.ac.uk>
#' @keywords internal
fm_hexagon_lattice_orig <- function(bnd,
                                    x_bin = 250, # 300 then running forever
                                    edge_len_n = 1) {
  stopifnot(x_bin / 2 > edge_len_n)
  crs <- fm_crs(bnd)
  fm_crs(bnd) <- NA
  # two separate grid and combine
  edge_len <- as.numeric((fm_bbox(bnd)[[1]][2] - fm_bbox(bnd)[[1]][1])) /
    (x_bin + 2 * edge_len_n) # two ends
  # sf_buffer to work on negative buffer to stay a distance from the boundary
  # Turn off S2 to avoid zig zag
  # suppressMessages(sf::sf_use_s2(FALSE))
  #  # st_buffer for edge_len x1
  bnd_inner <- sf::st_buffer(bnd, dist = -edge_len_n * edge_len)
  y_diff <- fm_bbox(bnd_inner)[[2]][2] - fm_bbox(bnd_inner)[[2]][1]
  x_diff <- fm_bbox(bnd_inner)[[1]][2] - fm_bbox(bnd_inner)[[1]][1]
  y_bin <- as.integer(y_diff / (sqrt(3) / 2 * edge_len))
  # TODO rep n, n-1, length
  h <- (sqrt(3) / 2 * edge_len) # height
  x_adj <- .5 * (x_diff - x_bin * edge_len)
  y_adj <- .5 * (y_diff - y_bin * h)
  # x
  x_1_ <- seq(
    fm_bbox(bnd_inner)[[1]][1] + x_adj,
    fm_bbox(bnd_inner)[[1]][2] - x_adj, edge_len
  )
  x_2_ <- seq(
    (fm_bbox(bnd_inner)[[1]][1] + x_adj + .5 * edge_len),
    (fm_bbox(bnd_inner)[[1]][2] - x_adj - .5 * edge_len),
    edge_len
  )
  y_1_ <- seq(
    fm_bbox(bnd_inner)[[2]][1] + y_adj,
    fm_bbox(bnd_inner)[[2]][2] - y_adj,
    by = 2 * h
  )
  y_2_ <- seq(
    fm_bbox(bnd_inner)[[2]][1] + y_adj + h,
    fm_bbox(bnd_inner)[[2]][2] - y_adj + h, 2 * h
  )

  x_1 <- rep(x_1_, times = length(y_1_))
  x_2 <- rep(x_2_, times = length(y_2_))
  y_1 <- rep(y_1_, each = length(x_1_))
  y_2 <- rep(y_2_, each = length(x_2_))

  mesh_df <- data.frame(x = c(x_1, x_2), y = c(y_1, y_2))
  # turn the mesh nodes into lattice sf
  lattice_sf <- sf::st_as_sf(mesh_df,
    coords = c("x", "y"),
    crs = sf::st_crs(bnd)
  )
  lattice_sfc <- sf::st_as_sfc(lattice_sf)
  pts_inside <- lengths(sf::st_intersects(lattice_sfc, bnd_inner)) != 0
  pts_lattice_sfc <- lattice_sfc[pts_inside]
  fm_crs(pts_lattice_sfc) <- fm_crs(bnd_inner) <- crs
  list(
    lattice = pts_lattice_sfc,
    edge_len = edge_len,
    bnd_inner = bnd_inner
  )
}


#' @title Create hexagon lattice points
#' @description `r lifecycle::badge("experimental")` from `0.3.0.9001`. Create
#'   hexagon lattice points within a boundary. By default, the hexagonal lattice
#'   is anchored at the coordinate system origin, so that grids with different
#'   but overlapping boundaries will have matching points.
#' @param bnd Boundary object (`sf` polygon or boundary `fm_segm` object)
#' @param edge_len Triangle edge length. Default `diff(fm_bbox(bnd)[[1]]) /
#'   250`.
#' @param buffer_n Number of triangle height multiples for buffer inside the
#'   boundary object to the start of the lattice. Default 0.49.
#' @param align Alignment of the hexagon lattice, either a length-2 numeric, or
#'   character, a `sf`/`sfc`/`sfg` object containing a single point), or
#'   `character`, default `"origin"`:
#' \describe{
#' \item{"origin"}{align the lattice with the coordinate system origin}
#' \item{"bbox"}{align the lattice with the midpoint of the bounding box of
#' `bnd`}
#' \item{"centroid"}{align the lattice with the centroid of the boundary,
#' `sf::st_centroid(bnd)`}
#' }
#' @param meta logical; if `TRUE`, return a list with diagnostic information
#' from the lattice construction (including the points themselves in `lattice`)
#' @return An `sfc` object with points, if `meta` is `FALSE` (default), or if
#' `meta=TRUE`, a list:
#' \describe{
#' \item{lattice}{`sfc` with lattice points}
#' \item{edge_len}{`numeric` with edge length}
#' \item{bnd_inner}{`sf` object with the inner boundary used to filter points
#'   outside of a `edge_len * buffer_n` distance from the boundary}
#' \item{grid_n}{`integer` with the number of points in each direction prior to
#' filtering}
#' \item{align}{`numeric` with the alignment coordinates of the hexagon lattice}
#' }
#' @author Man Ho Suen <M.H.Suen@@sms.ed.ac.uk>,
#'  Finn Lindgren <Finn.Lindgren@@gmail.com>
#' @seealso [fm_mesh_2d()]
#' @export
#' @examples
#' (m <- fm_mesh_2d(
#'   fm_hexagon_lattice(
#'     fmexample$boundary_sf[[1]],
#'     edge_len = 0.1 * 5
#'   ),
#'   max.edge = c(0.2, 1) * 5,
#'   boundary = fmexample$boundary_sf
#' ))
#'
#' (m2 <- fm_mesh_2d(
#'   fm_hexagon_lattice(
#'     fmexample$boundary_sf[[1]],
#'     edge_len = 0.1 * 5,
#'     align = "centroid"
#'   ),
#'   max.edge = c(0.2, 1) * 5,
#'   boundary = fmexample$boundary_sf
#' ))
#'
#' if (require("ggplot2", quietly = TRUE) &&
#'   require("patchwork", quietly = TRUE)) {
#'   ((ggplot() +
#'     geom_fm(data = m) +
#'     geom_point(aes(0, 0), col = "red")) |
#'     (ggplot() +
#'       geom_fm(data = m2) +
#'       geom_point(aes(0, 0), col = "red") +
#'       geom_sf(data = sf::st_centroid(fmexample$boundary_sf[[1]]))
#'     )
#'   )
#' }
fm_hexagon_lattice <- function(bnd,
                               edge_len = NULL,
                               buffer_n = 0.49,
                               align = "origin",
                               meta = FALSE) {
  #  stopifnot(x_bin / 2 > edge_len_n)
  if (inherits(bnd, "fm_segm")) {
    bnd <- fm_as_sfc(bnd)
  }
  crs <- fm_crs(bnd)
  # Avoid longlat S2 issues by removing the CRS information
  fm_crs(bnd) <- NA

  bbox <- fm_bbox(bnd)
  if (is.null(edge_len)) {
    edge_len <- diff(bbox[[1]]) / 250
  }

  if (is.character(align)) {
    align <- match.arg(align, c("origin", "bbox", "centroid"))
    if (align == "bbox") {
      # Align the hexagon lattice with the bounding box
      align <- c(
        (bbox[[1]][1] + bbox[[1]][2]) / 2,
        (bbox[[2]][1] + bbox[[2]][2]) / 2
      )
    } else if (align == "centroid") {
      # Align the hexagon lattice with the bounding box
      align <- sf::st_centroid(bnd)
    } else {
      # Align the hexagon lattice with the coordinate system origin
      align <- c(0, 0)
    }
  }

  if (inherits(align, c("sf", "sfc", "sfg"))) {
    align <- sf::st_coordinates(sf::st_centroid(align))
    align <- align[, intersect(colnames(align), c("X", "Y", "Z")), drop = FALSE]
  }
  stopifnot(is.numeric(align))
  origin <- align

  # Find covering rectangular grid extent
  h <- edge_len * sqrt(3) / 2
  grid_start <-
    c(
      floor((bbox[[1]][1] - origin[1]) / edge_len),
      floor((bbox[[2]][1] - origin[2]) / (2 * h)) * 2L
    )
  grid_end <-
    c(
      ceiling((bbox[[1]][2] - origin[1]) / edge_len),
      ceiling((bbox[[2]][2] - origin[2]) / (2 * h)) * 2L
    )
  grid_n <- grid_end - grid_start + 1L

  if (any(grid_n == 0L)) {
    # Empty lattice
    return()
  }

  # x
  x_1_ <- origin[1] + seq(
    grid_start[1] * edge_len,
    grid_end[1] * edge_len,
    length.out = grid_n[1]
  )
  x_2_ <- origin[1] + seq(
    (grid_start[1] + 0.5) * edge_len,
    (grid_end[1] - 0.5) * edge_len,
    length.out = grid_n[1] - 1L
  )
  y_1_ <- origin[2] + seq(
    grid_start[2] * h,
    grid_end[2] * h,
    length.out = (grid_n[2] + 1L) / 2L
  )
  y_2_ <- origin[2] + seq(
    (grid_start[2] + 1) * h,
    (grid_end[2] - 1) * h,
    length.out = (grid_n[2] + 1L) / 2L - 1L
  )

  x_1 <- rep(x_1_, times = length(y_1_))
  x_2 <- rep(x_2_, times = length(y_2_))
  y_1 <- rep(y_1_, each = length(x_1_))
  y_2 <- rep(y_2_, each = length(x_2_))

  mesh_df <- data.frame(x = c(x_1, x_2), y = c(y_1, y_2))
  # turn the mesh nodes into lattice sf
  lattice_sf <- sf::st_as_sf(mesh_df, coords = c("x", "y"), crs = crs)
  lattice_sfc <- sf::st_as_sfc(lattice_sf)

  # sf_buffer to work on negative buffer to stay a distance from the boundary
  # Turn off S2 to avoid zig zag
  # suppressMessages(sf::sf_use_s2(FALSE))
  #  # st_buffer for edge_len x1
  bnd_inner <- sf::st_buffer(bnd, dist = -buffer_n * h)
  fm_crs(bnd_inner) <- crs

  pts_inside <- lengths(sf::st_intersects(lattice_sfc, bnd_inner)) != 0
  pts_lattice_sfc <- lattice_sfc[pts_inside]

  if (meta) {
    return(list(
      lattice = pts_lattice_sfc,
      edge_len = edge_len,
      bnd_inner = bnd_inner,
      grid_n = grid_n,
      align = origin
    ))
  }

  pts_lattice_sfc
}


circle_mesh <- function(
  centre = c(0, 0),
  radius = 1,
  max.edge = NULL,
  layers = ceiling(sqrt(2) * radius / max.edge),
  crs = NULL,
  cumulative_shifts = FALSE,
  ...
) {
  centre <- fm_unify_coords(centre, crs = fm_crs(crs))

  layers <- max(1L, layers)
  if (is.null(max.edge)) {
    max.edge <- sqrt(2) * radius / layers
  }

  kk <- seq_len(layers)
  m <- ceiling(pi * kk * radius / (3 * layers * max.edge))
  n <- c(1, 6 * m)
  radii <- radius * c(0, kk) / layers
  shift <- (c(0, 0, 1 * (diff(m) == 0)))
  if (cumulative_shifts) {
    shift <- cumsum(shift)
  }
  angles <- lapply(seq_along(n), function(k) {
    (seq_len(n[k]) - 1 + 0.5 * shift[k]) * 2 * pi / n[k]
  })
  loc <- lapply(seq_along(n), function(k) {
    theta <- angles[[k]]
    x <- centre[1] + radii[k] * cos(theta)
    y <- centre[2] + radii[k] * sin(theta)
    cbind(x, y)
  })
  loc <- do.call(rbind, loc)
  bnd <- fm_segm(
    loc[sum(n[-(layers + 1L)]) + seq_len(n[layers + 1L]), , drop = FALSE],
    is.bnd = TRUE
  )

  mesh <- fm_mesh_2d(
    loc = loc,
    boundary = bnd,
    max.edge = max.edge * 2,
    crs = fm_crs(crs),
    ...
  )
  mesh$radius <- sqrt((mesh$loc[, 1] - centre[1])^2 +
    (mesh$loc[, 2] - centre[2])^2)
  mesh
}
