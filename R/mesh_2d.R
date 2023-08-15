#' @include deprecated.R

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
  return(matrix(0.0, 0, 3))
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
    while (ncol(x) < 3) {
      x <- cbind(x, 0.0)
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
  x <- fm_transform(
    sp::coordinates(x),
    crs0 = fm_crs(x),
    crs = crs,
    passthrough = TRUE
  )
  if (ncol(x) < 3) {
    while (ncol(x) < 3) {
      x <- cbind(x, 0.0)
    }
  } else if (ncol(x) > 3) {
    stop("Coordinates can have at mots 3 columns.")
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
fm_unify_coords.sf <- function(x, crs = NULL) {
  fm_unify_coords.sfc(sf::st_sfc(x), crs = crs)
}
#' @rdname fm_unify_coords
#' @export
fm_unify_coords.sfc <- function(x, crs = NULL) {
  loc <- sf::st_coordinates(x)
  loc <- loc[, intersect(colnames(loc), c("X", "Y", "Z")), drop = FALSE]
  x <- fm_transform(
    loc,
    crs0 = fm_crs(x),
    crs = crs,
    passthrough = TRUE
  )
  if (ncol(x) < 3) {
    while (ncol(x) < 3) {
      x <- cbind(x, 0.0)
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
    .n,
    .loc) {
  options <- list(cutoff = cutoff)
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
  if (inherits(extend, "list")) {
    cet_sides <- ifelse(is.null(extend$n), 16L, as.integer(extend$n))
    cet_margin <- ifelse(is.null(extend$offset), -0.1, extend$offset)
  }
  options <- c(options, list(cet_sides = cet_sides, cet_margin = cet_margin))

  if (inherits(refine, "list")) {
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
#' @param loc Input coordinates that should be part of the mesh
#' @param tv Initial triangulation, as a N-by-3 indec vector into `loc`
#' @param boundary,interior Objects supported by [fm_as_segm()].
#' If `boundary` is `numeric`, `fm_nonconvex_hull(loc, convex = boundary)` is
#' used.
#' @param extend `logical` or `list` specifying whether to extend the
#' data region, with parameters \describe{ \item{list("n")}{the number of edges
#' in the extended boundary (default=8)} \item{list("offset")}{the extension
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
      loc.globe <- fm_transform(loc.globe, crs = crs, passthrough = TRUE, crs0 = crs.globe)
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

  loc <- rbind(loc.bnd, loc.int, lattice$loc, loc)

  options <- handle_rcdt_options_inla(
    extend = extend,
    refine = refine,
    cutoff = cutoff,
    qulity.spec = quality.spec,
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
        vt = idx_C2R(result[["vt"]]),
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
      crs = fm_crs(crs),
      n = nrow(result[["s"]])
    ),
    class = c("fm_mesh_2d", "inla.mesh")
  )

  remap_unused <- function(mesh) {
    ## Remap indices to remove unused vertices
    used <- !is.na(mesh$graph$vt)
    if (!all(used)) {
      used <- which(used)
      idx.map <- rep(NA, nrow(mesh$loc))
      idx.map[used] <- seq_len(length(used))
      mesh$loc <- mesh$loc[used, , drop = FALSE]
      mesh$graph$tv <-
        matrix(idx.map[as.vector(mesh$graph$tv)], nrow(mesh$graph$tv), 3)
      mesh$graph$vt <- mesh$graph$vt[used, , drop = FALSE]
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
  return(mesh)
}


#' @title Make a 2D mesh object
#' @export
#' @param ... Currently passed on to `fm_mesh_2d_inla`
#' @family object creation and conversion
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
#' nodes.  Can alternatively be a `SpatialPoints` or
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
#' @param interior one object supported by [fm_as_segm()]
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
#' @param plot.delay On Linux (and Mac if appropriate X11 libraries are
#' installed), specifying a nonnegative numeric value activates a rudimentary
#' plotting system in the underlying `fmesher` program, showing the
#' triangulation algorithm at work, with waiting time factor `plot.delay`
#' between each step.
#'
#' On all systems, specifying any negative value activates displaying the
#' result after each step of the multi-step domain extension algorithm.
#' @param crs An optional [fm_crs()], `sf::crs` or `sp::CRS` object
#' @return An `inla.mesh` object.
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @seealso [fm_rcdt_2d()], [fm_mesh_2d()], [fm_delaunay_2d()],
#' [fm_nonconvex_hull()], [fm_extensions()], [fm_refine()]
fm_mesh_2d_inla <- function(loc = NULL, ## Points to include in final triangulation
                            loc.domain = NULL, ## Points that determine the automatic domain
                            offset = NULL, ## Size of automatic extensions
                            n = NULL, ## Sides of automatic extension polygons
                            boundary = NULL, ## User-specified domains (list of length 2)
                            interior = NULL, ## User-specified constraints for the inner domain
                            max.edge = NULL,
                            min.angle = NULL, ## Angle constraint for the entire domain
                            cutoff = 1e-12, ## Only add input points further apart than this
                            max.n.strict = NULL,
                            max.n = NULL,
                            plot.delay = NULL,
                            crs = NULL,
                            ...) {
  ## plot.delay: Do plotting.
  ## NULL --> No plotting
  ## <0  --> Intermediate meshes displayed at the end
  ## >0   --> Dynamical fmesher plotting


  if ((missing(max.edge) || is.null(max.edge)) &&
    (missing(max.n.strict) || is.null(max.n.strict)) &&
    (missing(max.n) || is.null(max.n))) {
    max.edge <- NA
    #    stop("At least one of max.edge, max.n.strict, and max.n must be specified")
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
  interior <- fm_as_segm(interior)

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
  if (any(offset < 0) &&
    (fm_diameter(loc) +
      fm_diameter(loc.domain) +
      fm_diameter(interior) == 0.0)) {
    offset[offset < 0] <- 1
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
  if (missing(plot.delay) || is.null(plot.delay)) {
    plot.delay <- NULL
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

  ## Unify the dimensionality of the boundary&interior segments input
  ## and optionally transform coordinates.
  boundary <- unify_segm_coords(boundary, crs = crs)
  interior <- unify_segm_coords(interior, crs = crs)

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
      plot.delay = plot.delay,
      crs = crs
    )

  ## Save the resulting boundary
  boundary1 <- fm_segm(mesh1, boundary = TRUE)
  interior1 <- fm_segm(mesh1, boundary = FALSE)

  if (!is.null(plot.delay) && (plot.delay < 0)) {
    plot(mesh1)
  }

  ## Triangulate inner domain
  mesh2 <-
    fm_rcdt_2d(
      loc = loc,
      boundary = boundary1,
      interior = interior1,
      cutoff = cutoff,
      extend = FALSE, ## Should have no effect
      refine =
        list(
          min.angle = min.angle[1],
          max.edge = max.edge[1],
          max.edge.extra = max.edge[1],
          max.n.strict = max.n.strict[1],
          max.n = max.n[1]
        ),
      plot.delay = plot.delay,
      crs = crs
    )

  boundary2 <- fm_segm(mesh2, boundary = TRUE)
  interior2 <- fm_segm(mesh2, boundary = FALSE)

  if (!is.null(plot.delay) && (plot.delay < 0)) {
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
      interior = fm_segm(boundary2, interior2),
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
      plot.delay = plot.delay,
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
  for (k in seq_along(interior)) {
    if (!is.null(interior[[k]])) {
      segm.loc <- rbind(segm.loc, interior[[k]]$loc)
    }
  }
  if (nrow(segm.loc) > 0) {
    proj <- fm_evaluator(mesh3, loc = segm.loc)$proj
    mesh3$idx$segm <- rep(NA, nrow(segm.loc))
    if (any(proj$ok)) {
      t.idx <- proj$t[proj$ok]
      tv.idx <- max.col(proj$bary[proj$ok, , drop = FALSE],
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

  if (!is.null(plot.delay) && (plot.delay < 0)) {
    plot(mesh3)
  }

  return(mesh3)
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
