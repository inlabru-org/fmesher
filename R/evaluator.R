# Point/mesh connection methods ####

#' @title Methods for projecting to/from an inla.mesh
#'
#' @description Calculate evaluation information and/or evaluate a function
#' defined on a mesh or function space.
#'
#' @param mesh An `inla.mesh` or `inla.mesh.1d` object.
#' @param loc Projection locations.  Can be a matrix, `SpatialPoints`,
#' `SpatialPointsDataFrame`, `sf`, `sfc`, or `sfg` object.
#' @param lattice An [fm_lattice_2d()] object.
#' @param xlim X-axis limits for a lattice. For R2 meshes, defaults to covering
#' the domain.
#' @param ylim Y-axis limits for a lattice. For R2 meshes, defaults to covering
#' the domain.
#' @param dims Lattice dimensions.
#' @param projector An `fm_evaluator` object.
#' @param field Basis function weights, one per mesh basis function, describing
#' the function to be evaluated at the projection locations
#' @param projection One of `c("default", "longlat", "longsinlat",
#' "mollweide")`.
#' @param crs An optional CRS or inla.CRS object associated with `loc`
#' and/or `lattice`.
#' @param \dots Additional arguments passed on to methods.
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @seealso [fm_mesh_2d()], [fm_mesh_1d()],
#' [fm_lattice_2d()]
#' @examples
#' if (fm_safe_inla()) {
#'   n <- 20
#'   loc <- matrix(runif(n * 2), n, 2)
#'   mesh <- INLA::inla.mesh.create(loc, refine = list(max.edge = 0.05))
#'   proj <- fm_evaluator(mesh)
#'   field <- cos(mesh$loc[, 1] * 2 * pi * 3) * sin(mesh$loc[, 2] * 2 * pi * 7)
#'   image(proj$x, proj$y, fm_evaluate(proj, field))
#' }
#' \donttest{
#' if (fm_safe_inla() &&
#'   require("ggplot2") &&
#'   fm_safe_sp() &&
#'   require("ggpolypath") &&
#'   require("inlabru")) {
#'   ggplot() +
#'     gg(mesh, col = field)
#' }
#' }
#'
#' @name fm_evaluate
#' @rdname fm_evaluate
NULL

#' @describeIn fm_evaluate
#' Returns the field function evaluated at the locations determined by an
#' `fm_evaluator` object. `fm_evaluate(mesh, field = field, ...)` is a
#' shortcut to `fm_evaluate(fm_evaluator(mesh, ...), field = field)`.
#' @export fm_evaluate
fm_evaluate <- function(...) {
  UseMethod("fm_evaluate")
}

#' @export
#' @rdname fm_evaluate
fm_evaluate.inla.mesh <- function(mesh, field, ...) {
  fm_evaluate.fm_mesh_2d(fm_as_mesh_2d(mesh), field = field, ...)
}
#' @export
#' @rdname fm_evaluate
fm_evaluate.fm_mesh_2d <- function(mesh, field, ...) {
  if (missing(field) || is.null(field)) {
    lifecycle::deprecate_stop(
      "0.0.1",
      "fm_evaluate(field = ' must not be missing or NULL.')",
      "fm_evaluator()"
    )
  }

  proj <- fm_evaluator(mesh, ...)
  fm_evaluate(proj, field = field)
}
#' @export
#' @rdname fm_evaluate
fm_evaluate.fm_tensor <- function(mesh, field, ...) {
  if (missing(field) || is.null(field)) {
    lifecycle::deprecate_stop(
      "0.0.1",
      "fm_evaluate(field = ' must not be missing or NULL.')",
      "fm_evaluator()"
    )
  }

  proj <- fm_evaluator(mesh, ...)
  fm_evaluate(proj, field = field)
}


#' @export
#' @rdname fm_evaluate
fm_evaluate.inla.mesh.1d <- function(mesh, field, ...) {
  fm_evaluate.fm_mesh_1d(fm_as_mesh_1d(mesh), field = field, ...)
}
#' @export
#' @rdname fm_evaluate
fm_evaluate.fm_mesh_1d <- function(mesh, field, ...) {
  if (missing(field) || is.null(field)) {
    lifecycle::deprecate_stop(
      "0.0.1",
      "fm_evaluate(field = ' must not be missing or NULL.')",
      "fm_evaluator()"
    )
  }

  proj <- fm_evaluator(mesh, ...)
  fm_evaluate(proj, field)
}


#' @export
#' @rdname fm_evaluate
fm_evaluate.fm_evaluator <-
  function(projector, field, ...) {
    if (is.data.frame(field)) {
      field <- as.matrix(field)
    }

    if (is.null(dim(field))) {
      if (is.null(projector$lattice)) {
        data <- as.vector(projector$proj$A %*% as.vector(field))
        data[!projector$proj$ok] <- NA
        return(data)
      } else {
        data <- as.vector(projector$proj$A %*% as.vector(field))
        data[!projector$proj$ok] <- NA
        return(matrix(
          data,
          projector$lattice$dims[1],
          projector$lattice$dims[2]
        ))
      }
    } else if (inherits(field, "sparseMatrix")) {
      data <- projector$proj$A %*% field
      data[!projector$proj$ok, ] <- NA
      return(data)
    } else {
      data <- as.matrix(projector$proj$A %*% field)
      data[!projector$proj$ok, ] <- NA
      return(data)
    }
  }


#' @describeIn fm_evaluate
#' Returns an `fm_evaluator` list object with evaluation information.
#' The `proj` element contains a mapping matrix `A` and a logical vector `ok`,
#' that indicates which locations were mappable to the input mesh.
#' For `fm_mesh_2d` and `inla.mesh`
#' input, `proj` also contains a matrix `bary` and vector `t`, with the
#' barycentric coordinates within the triangle each input location falls in.
#' @export
fm_evaluator <- function(...) {
  UseMethod("fm_evaluator")
}


#' @export
#' @rdname fm_evaluate
fm_evaluator_mesh_2d <- function(mesh, loc = NULL, crs = NULL, ...) {
  smorg <- fm_bary(mesh, loc = loc, crs = crs)
  ti <- matrix(0L, NROW(loc), 1)
  ti[, 1L] <- smorg$t
  b <- smorg$bary

  ok <- !is.na(ti[, 1L])

  ii <- which(ok)
  A <- (Matrix::sparseMatrix(
    dims = c(NROW(loc), mesh$n),
    i = rep(ii, 3),
    j = as.vector(mesh$graph$tv[ti[ii, 1L], ]),
    x = as.vector(b[ii, ])
  ))

  list(t = ti, bary = b, A = A, ok = ok)
}


#' @export
#' @rdname fm_evaluate
fm_evaluator_mesh_1d <- function(mesh, loc, ...) {
  A <- INLA::inla.mesh.1d.A(mesh, loc = loc)

  return(list(
    A = A,
    ok = (loc >= mesh$interval[1]) & (loc <= mesh$interval[2])
  ))
}




#' @describeIn fm_evaluate
#' Creates an [fm_lattice_2d()] object, by default covering the input mesh.
#' @export
fm_evaluator_lattice <- function(mesh,
                                 xlim = NULL,
                                 ylim = NULL,
                                 dims = c(100, 100),
                                 projection = NULL,
                                 crs = NULL,
                                 ...) {
  stopifnot(inherits(mesh, c("fm_mesh_2d", "inla.mesh")))
  if (identical(mesh$manifold, "R2") &&
    (is.null(mesh$crs) || is.null(crs))) {
    units <- "default"
    lim <- list(
      xlim = if (is.null(xlim)) range(mesh$loc[, 1]) else xlim,
      ylim = if (is.null(ylim)) range(mesh$loc[, 2]) else ylim
    )
  } else if (identical(mesh$manifold, "S2") &&
    (is.null(mesh$crs) || is.null(crs))) {
    projection <-
      match.arg(projection, c(
        "longlat", "longsinlat",
        "mollweide"
      ))
    units <- projection
    lim <- INLA::inla.mesh.map.lim(loc = mesh$loc, projection = projection)
  } else {
    lim <- fm_crs_bounds(crs)
    if (identical(mesh$manifold, "R2")) {
      lim0 <- list(
        xlim = if (is.null(xlim)) range(mesh$loc[, 1]) else xlim,
        ylim = if (is.null(ylim)) range(mesh$loc[, 2]) else ylim
      )
      lim$xlim[1] <- max(lim$xlim[1], lim0$xlim[1])
      lim$xlim[2] <- min(lim$xlim[2], lim0$xlim[2])
      lim$ylim[1] <- max(lim$ylim[1], lim0$ylim[1])
      lim$ylim[2] <- min(lim$ylim[2], lim0$ylim[2])
    }
  }
  if (missing(xlim) && is.null(xlim)) {
    xlim <- lim$xlim
  }
  if (missing(ylim) && is.null(ylim)) {
    ylim <- lim$ylim
  }
  x <- seq(xlim[1], xlim[2], length.out = dims[1])
  y <- seq(ylim[1], ylim[2], length.out = dims[2])
  if (is.null(mesh$crs) || is.null(crs)) {
    lattice <- fm_lattice_2d(x = x, y = y, units = units)
  } else {
    lattice <- fm_lattice_2d(x = x, y = y, crs = crs)
  }
  lattice
}

#' @export
#' @describeIn fm_evaluate The `...` arguments are passed on to `fm_evaluator_lattice()`
#' if no `loc` or `lattice` is provided.
fm_evaluator.inla.mesh <- function(mesh,
                                   loc = NULL,
                                   lattice = NULL,
                                   crs = NULL,
                                   ...) {
  fm_evaluator.fm_mesh_2d(
    fm_as_mesh_2d(mesh),
    loc = loc,
    lattice = lattice,
    crs = crs,
    ...
  )
}
#' @export
#' @describeIn fm_evaluate The `...` arguments are passed on to `fm_evaluator_lattice()`
#' if no `loc` or `lattice` is provided.
fm_evaluator.fm_mesh_2d <- function(mesh,
                                    loc = NULL,
                                    lattice = NULL,
                                    crs = NULL,
                                    ...) {
  if (missing(loc) || is.null(loc)) {
    if (missing(lattice) || is.null(lattice)) {
      lattice <- fm_evaluator_lattice(mesh,
        crs = crs,
        ...
      )
    }
    dims <- lattice$dims
    x <- lattice$x
    y <- lattice$y
    crs <- lattice$crs

    if (is.null(mesh$crs) || is.null(lattice$crs)) {
      proj <- fm_evaluator_mesh_2d(mesh, lattice$loc)
    } else {
      proj <- fm_evaluator_mesh_2d(mesh,
        loc = lattice$loc,
        crs = lattice$crs
      )
    }
    projector <- list(x = x, y = y, lattice = lattice, loc = NULL, proj = proj, crs = crs)
    class(projector) <- "fm_evaluator"
  } else {
    proj <- fm_evaluator_mesh_2d(mesh, loc = loc, crs = crs)
    projector <- list(x = NULL, y = NULL, lattice = NULL, loc = loc, proj = proj, crs = crs)
    class(projector) <- "fm_evaluator"
  }

  return(projector)
}


#' @export
#' @rdname fm_evaluate
fm_evaluator.inla.mesh.1d <- function(mesh,
                                      loc = NULL,
                                      xlim = mesh$interval,
                                      dims = 100,
                                      ...) {
  fm_evaluator.fm_mesh_1d(fm_as_mesh_1d(mesh), loc = loc, xlim = xlim, dims = dims, ...)
}
#' @export
#' @rdname fm_evaluate
fm_evaluator.fm_mesh_1d <- function(mesh,
                                    loc = NULL,
                                    xlim = mesh$interval,
                                    dims = 100,
                                    ...) {
  if (missing(loc) || is.null(loc)) {
    loc <- seq(xlim[1], xlim[2], length.out = dims[1])
  }

  proj <- fm_evaluator_mesh_1d(mesh, loc)
  projector <- list(x = loc, lattice = NULL, loc = loc, proj = proj)
  class(projector) <- "fm_evaluator"

  return(projector)
}


#' @param x [fm_tensor()] object
#' @export
#' @rdname fm_evaluate
fm_evaluator.fm_tensor <- function(x,
                                   loc,
                                   ...) {
  if (length(loc) != length(x[["fun_spaces"]])) {
    stop("")
  }
  if (is.null(names(loc))) {
    names(loc) <- names(x[["fun_spaces"]])
  } else if (!setequal(names(x[["fun_spaces"]]), names(loc))) {
    stop("")
  }
  proj <- lapply(
    names(x[["fun_spaces"]]),
    function(k) {
      fm_evaluator(x[["fun_spaces"]][[k]], loc = loc[[k]])$proj
    }
  )

  # Combine the matrices
  # (A1, A2, A3) -> rowkron(A3, rowkron(A2, A1))
  A <- proj[[1]][["A"]]
  ok <- proj[[1]][["ok"]]
  for (k in seq_len(length(x[["fun_spaces"]]) - 1)) {
    A <- fm_row_kron(proj[[k + 1]][["A"]], A)
    ok <- proj[[k + 1]][["ok"]] & ok
  }

  structure(
    list(proj = list(A = A, ok = ok)),
    class = "fm_evaluator"
  )
}





#' Check which mesh triangles are inside a polygon
#'
#' Wrapper for the [sf::st_contains()] (previously `sp::over()`) method to find triangle centroids
#' or vertices inside `sf` or `sp` polygon objects
#'
#' @param x geometry (typically an `sf` or `sp::SpatialPolygons` object) for the queries
#' @param y an [fm_mesh_2d()] or `inla.mesh` object
#' @param \dots Passed on to other methods
#' @param type the query type; either `'centroid'` (default, for triangle centroids),
#' or `'vertex'` (for mesh vertices)
#'
#' @return List of vectors of triangle indices (when `type` is `'centroid'`) or
#' vertex indices (when `type` is `'vertex'`). The list has one entry per row of the `sf` object.
#' Use `unlist(fm_contains(...))` if the combined union is needed.
#'
#' @author Haakon Bakka, \email{bakka@@r-inla.org}, and Finn Lindgren \email{finn.lindgren@@gmail.com}
#'
#' @examples
#' if (fm_safe_inla() &&
#'   fm_safe_sp()) {
#'   # Create a polygon and a mesh
#'   obj <- sp::SpatialPolygons(
#'     list(sp::Polygons(
#'       list(sp::Polygon(rbind(
#'         c(0, 0),
#'         c(50, 0),
#'         c(50, 50),
#'         c(0, 50)
#'       ))),
#'       ID = 1
#'     )),
#'     proj4string = fm_CRS("longlat_globe")
#'   )
#'   mesh <- INLA::inla.mesh.create(globe = 2, crs = fm_crs("sphere"))
#'
#'   ## 3 vertices found in the polygon
#'   fm_contains(obj, mesh, type = "vertex")
#'
#'   ## 3 triangles found in the polygon
#'   fm_contains(obj, mesh)
#'
#'   ## Multiple transformations can lead to slightly different results due to edge cases
#'   ## 4 triangles found in the polygon
#'   fm_contains(
#'     obj,
#'     fm_transform(mesh, crs = fm_crs("mollweide_norm"))
#'   )
#' }
#'
#' @export
fm_contains <- function(x, y, ...) {
  UseMethod("fm_contains")
}

#' @rdname fm_contains
#' @export
fm_contains.Spatial <- function(x, y, ...) {
  fm_contains(sf::st_as_sf(x), y = y, ...)
}

#' @rdname fm_contains
#' @export
fm_contains.sf <- function(x, y, ...) {
  fm_contains(sf::st_geometry(x), y = y, ...)
}

#' @rdname fm_contains
#' @export
fm_contains.sfc <- function(x, y, ..., type = c("centroid", "vertex")) {
  if (!inherits(y, c("fm_mesh_2d", "inla.mesh"))) {
    stop(paste0(
      "'y' must be an 'fm_mesh_2d' or 'inla.mesh' object, not '",
      paste0(class(y), collapse = ", "),
      "'."
    ))
  }

  type <- match.arg(type)
  if (identical(type, "centroid")) {
    ## Extract triangle centroids
    points <- (y$loc[y$graph$tv[, 1], , drop = FALSE] +
      y$loc[y$graph$tv[, 2], , drop = FALSE] +
      y$loc[y$graph$tv[, 3], , drop = FALSE]) / 3
  } else if (identical(type, "vertex")) {
    ## Extract vertices
    points <- y$loc
  }
  if (identical(y$manifold, "S2")) {
    points <- points / rowSums(points^2)^0.5
  }
  ## Convert to sf points
  ## Extract coordinate system information
  if (identical(y$manifold, "S2")) {
    crs <- fm_crs("sphere")
  } else {
    crs <- fm_crs(y)
  }
  crs_x <- fm_crs(x)
  ## Create sfc_POINT object and transform the coordinates.
  points <- sf::st_as_sf(as.data.frame(points),
    coords = seq_len(ncol(points)),
    crs = crs
  )
  if (!fm_crs_is_null(crs) &&
    !fm_crs_is_null(crs_x)) {
    ## Convert to the target object CRS
    points <- fm_transform(points, crs = crs_x)
  }

  ## Find indices:
  ids <- sf::st_contains(x, points, sparse = TRUE)

  ids
}


#' @title Query if points are inside a mesh
#'
#' @description
#'  Queries whether each input point is within a mesh or not.
#'
#' @param x A set of points of a class supported by `fm_evaluator(y, loc = x)`
#' @param y An `inla.mesh`
#' @param \dots Currently unused
#' @returns A logical vector
#' @examples
#' \dontrun{
#' if (fm_safe_inla(quietly = TRUE)) {
#'   # Load Gorilla data
#'
#'   data("gorillas", package = "inlabru")
#'
#'   # Check if all Gorilla nests are inside the mesh
#'
#'   all(fm_is_within(gorillas$nests, gorillas$mesh))
#'
#'   # Also works for locations not stored as SpatialPoints object
#'
#'   loc <- coordinates(gorillas$nests)
#'   all(fm_is_within(loc, gorillas$mesh))
#' }
#' }
#' @export
fm_is_within <- function(x, y, ...) {
  UseMethod("fm_is_within")
}

#' @rdname fm_is_within
#' @export
fm_is_within.default <- function(x, y, ...) {
  fm_evaluator(y, loc = x)$proj$ok
}


#' @title Compute mapping matrix between mesh function space and points
#'
#' @description
#'  Computes the mapping matrix between a function space on a mesh, and locations.
#'
#' @param x An mesh object supported by an [fm_evaluator()] class
#' @param y A set of points of a class supported by `fm_evaluator(x, loc = y)`
#' @param \dots Currently unused
#' @returns A `sparseMatrix`
#' @examples
#' \dontrun{
#' if (fm_safe_inla(quietly = TRUE)) {
#'   # Load Gorilla data
#'
#'   data("gorillas", package = "inlabru")
#'
#'   # Compute mapping matrix
#'   str(fm_A(gorillas$mesh, gorillas$nests))
#' }
#' }
#' @export
fm_A <- function(x, ...) {
  UseMethod("fm_A")
}

#' @rdname fm_A
#' @export
fm_A.default <- function(x, y, ...) {
  fm_evaluator(x, loc = y)$proj$A
}
