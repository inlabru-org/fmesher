# Point/mesh connection methods ####

# fm_evaluate ####

#' @title Methods for projecting to/from mesh objects
#'
#' @description Calculate evaluation information and/or evaluate a function
#' defined on a mesh or function space.
#'
#' @param mesh An [fm_mesh_1d], [fm_mesh_2d], or other object supported by a
#' sub-method.
#' @param loc Projection locations.  Can be a matrix, `SpatialPoints`,
#' `SpatialPointsDataFrame`, `sf`, `sfc`, or `sfg` object.
#' @param lattice An [fm_lattice_2d()] object.
#' @param xlim X-axis limits for a lattice. For R2 meshes, defaults to covering
#' the domain.
#' @param ylim Y-axis limits for a lattice. For R2 meshes, defaults to covering
#' the domain.
#' @param dims Lattice dimensions.
#' @param projector An `fm_evaluator` object.
#' @param basis An [fm_basis] object.
#' @param field Basis function weights, one per mesh basis function, describing
#' the function to be evaluated at the projection locations
#' @param projection One of `c("default", "longlat", "longsinlat",
#' "mollweide")`.
#' @param crs An optional CRS or inla.CRS object associated with `loc`
#' and/or `lattice`.
#' @param \dots Additional arguments passed on to methods.
#' @author Finn Lindgren <Finn.Lindgren@@gmail.com>
#' @seealso [fm_mesh_2d()], [fm_mesh_1d()],
#' [fm_lattice_2d()]
#' @examples
#' if (TRUE) {
#'   n <- 20
#'   loc <- matrix(runif(n * 2), n, 2)
#'   mesh <- fm_rcdt_2d_inla(loc, refine = list(max.edge = 0.05))
#'   proj <- fm_evaluator(mesh)
#'   field <- cos(mesh$loc[, 1] * 2 * pi * 3) * sin(mesh$loc[, 2] * 2 * pi * 7)
#'   image(proj$x, proj$y, fm_evaluate(proj, field))
#' }
#' \donttest{
#' # if (require("ggplot2") &&
#' #  require("ggpolypath")) {
#' #  ggplot() +
#' #    gg(data = fm_as_sfc(mesh), col = field)
#' # }
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
#' @returns A vector or matrix of the evaluated function
fm_evaluate <- function(...) {
  UseMethod("fm_evaluate")
}

#' @export
#' @describeIn fm_evaluate The default method calls
#' `proj = fm_evaluator(mesh, ...)`, followed by `fm_evaluate(proj, field)`.
fm_evaluate.default <- function(mesh, field, ...) {
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
fm_evaluate.fm_evaluator <-
  function(projector, field, ...) {
    if (is.data.frame(field)) {
      field <- as.matrix(field)
    }

    data <- fm_evaluate(fm_basis(projector, full = TRUE), field = field)

    if (is.null(dim(field)) &&
      !is.null(projector$lattice)) {
      return(array(
        data,
        dim = projector$lattice$dims
      ))
    }

    data
  }

#' @export
#' @rdname fm_evaluate
fm_evaluate.fm_basis <-
  function(basis, field, ...) {
    if (is.data.frame(field)) {
      field <- as.matrix(field)
    }

    if (is.null(dim(field))) {
      data <- as.vector(basis$A %*% as.vector(field))
      data[!basis$ok] <- NA
    } else if (inherits(field, "sparseMatrix")) {
      data <- basis$A %*% field
      data[!basis$ok, ] <- NA
    } else {
      data <- as.matrix(basis$A %*% field)
      data[!basis$ok, ] <- NA
    }

    data
  }


# fm_evaluator ####

#' @describeIn fm_evaluate
#' Returns an `fm_evaluator` list object with evaluation information.
#' The `proj` element is a `fm_basis` object, containing (at least)
#' a mapping matrix `A` and a logical vector `ok`, that indicates which
#' locations were mappable to the input mesh.
#' For `fm_mesh_2d`
#' input, `proj` also contains a `bary` [fm_bary] object, with the
#' barycentric coordinates within the triangle each input location falls in.
#' @export
#' @returns An `fm_evaluator` object
fm_evaluator <- function(...) {
  UseMethod("fm_evaluator")
}

#' @export
#' @describeIn fm_evaluate The default method calls `fm_basis` and creates
#' a basic `fm_evaluator` object
fm_evaluator.default <- function(...) {
  structure(
    list(proj = fm_basis(..., full = TRUE)),
    class = "fm_evaluator"
  )
}


#' @export
#' @describeIn fm_evaluate The `...` arguments are passed on to
#'   `fm_evaluator_lattice()` if no `loc` or `lattice` is provided.
fm_evaluator.fm_mesh_3d <- function(mesh,
                                    loc = NULL,
                                    lattice = NULL,
                                    dims = NULL,
                                    ...) {
  if (missing(loc) || is.null(loc)) {
    if (missing(lattice) || is.null(lattice)) {
      lattice <- fm_evaluator_lattice(mesh,
        dims = dims,
        ...
      )
    }
    proj <- fm_basis(mesh, lattice$loc, full = TRUE)
    projector <-
      structure(
        list(
          loc = NULL,
          lattice = lattice,
          proj = proj
        ),
        class = "fm_evaluator"
      )
  } else {
    proj <- fm_basis(mesh, loc, full = TRUE)
    projector <-
      structure(
        list(
          loc = loc,
          lattice = NULL,
          proj = proj
        ),
        class = "fm_evaluator"
      )
  }

  projector
}


#' @export
#' @describeIn fm_evaluate The `...` arguments are passed on to
#'   `fm_evaluator_lattice()` if no `loc` or `lattice` is provided.
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
      proj <- fm_basis_mesh_2d(mesh, lattice$loc)
    } else {
      proj <- fm_basis_mesh_2d(mesh,
        loc = lattice$loc,
        crs = lattice$crs
      )
    }
    projector <-
      structure(
        list(
          x = x,
          y = y,
          lattice = lattice,
          loc = NULL,
          proj = proj,
          crs = crs
        ),
        class = "fm_evaluator"
      )
  } else {
    proj <- fm_basis_mesh_2d(mesh, loc = loc, crs = crs)
    projector <-
      structure(
        list(
          x = NULL,
          y = NULL,
          lattice = NULL,
          loc = loc,
          proj = proj,
          crs = crs
        ),
        class = "fm_evaluator"
      )
  }

  projector
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

  proj <- fm_basis_mesh_1d(mesh, loc)
  projector <-
    structure(
      list(
        x = loc,
        lattice = NULL,
        loc = loc,
        proj = proj
      ),
      class = "fm_evaluator"
    )

  projector
}


#' @describeIn fm_evaluate
#' Create a lattice object by default covering the input mesh.
#' @export
fm_evaluator_lattice <- function(mesh,
                                 ...) {
  UseMethod("fm_evaluator_lattice")
}

#' @describeIn fm_evaluate
#' Creates an [fm_lattice_2d()] object, by default covering the input mesh.
#' @export
fm_evaluator_lattice.default <- function(mesh,
                                         dims = 100,
                                         ...) {
  bbox <- fm_bbox(mesh)
  if (length(dims) == 1L) {
    dims <- rep(dims, length(bbox))
  }
  if (length(bbox) != length(dims)) {
    stop("The length of 'dims' must match the length of 'fm_bbox(mesh)'.")
  }
  fm_lattice_Nd(bbox, dims = dims)
}

#' @describeIn fm_evaluate
#' Creates an [fm_lattice_Nd()] object, by default covering the input mesh.
#' @export
fm_evaluator_lattice.fm_bbox <- function(mesh,
                                         dims = 100,
                                         ...) {
  bbox <- mesh
  if (length(dims) == 1L) {
    dims <- rep(dims, length(bbox))
  }
  if (length(bbox) != length(dims)) {
    stop("The length of 'dims' must match the length of 'fm_bbox(mesh)'.")
  }
  fm_lattice_Nd(bbox, dims = dims)
}

#' @describeIn fm_evaluate
#' Creates an [fm_lattice_2d()] object, by default covering the input mesh.
#' @export
fm_evaluator_lattice.fm_mesh_2d <- function(mesh,
                                            xlim = NULL,
                                            ylim = NULL,
                                            dims = c(100, 100),
                                            projection = NULL,
                                            crs = NULL,
                                            ...) {
  if (fm_manifold(mesh, "R2") &&
    (is.null(mesh$crs) || is.null(crs))) {
    units <- "default"
    lim <- list(
      xlim = if (is.null(xlim)) range(mesh$loc[, 1]) else xlim,
      ylim = if (is.null(ylim)) range(mesh$loc[, 2]) else ylim
    )
  } else if (fm_manifold(mesh, "S2") &&
    (is.null(mesh$crs) || is.null(crs))) {
    projection <-
      match.arg(projection, c(
        "longlat", "longsinlat",
        "mollweide"
      ))
    units <- projection
    lim <- fm_mesh_2d_map_lim(loc = mesh$loc, projection = projection)
  } else {
    lim <- fm_crs_bounds(crs)
    if (fm_manifold(mesh, "R2")) {
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


# fm_contains ####

#' Check which mesh triangles are inside a polygon
#'
#' Wrapper for the [sf::st_contains()] (previously `sp::over()`) method to find
#' triangle centroids or vertices inside `sf` or `sp` polygon objects
#'
#' @param x geometry (typically an `sf` or `sp::SpatialPolygons` object) for the
#'   queries
#' @param y an [fm_mesh_2d()] object
#' @param \dots Passed on to other methods
#' @param type the query type; either `'centroid'` (default, for triangle
#'   centroids), or `'vertex'` (for mesh vertices)
#'
#' @returns List of vectors of triangle indices (when `type` is `'centroid'`) or
#'   vertex indices (when `type` is `'vertex'`). The list has one entry per row
#'   of the `sf` object. Use `unlist(fm_contains(...))` if the combined union is
#'   needed.
#'
#' @author Haakon Bakka, <bakka@@r-inla.org>, and Finn Lindgren
#'   <Finn.Lindgren@@gmail.com>
#'
#' @examples
#' # Create a polygon and a mesh
#' obj <- sf::st_sfc(
#'   sf::st_polygon(
#'     list(rbind(
#'       c(0, 0),
#'       c(50, 0),
#'       c(50, 50),
#'       c(0, 50),
#'       c(0, 0)
#'     ))
#'   ),
#'   crs = fm_crs("longlat_globe")
#' )
#' mesh <- fm_rcdt_2d_inla(globe = 2, crs = fm_crs("sphere"))
#'
#' ## 2 vertices found in the polygon
#' fm_contains(obj, mesh, type = "vertex")
#'
#' ## 3 triangles found in the polygon
#' fm_contains(obj, mesh)
#'
#' ## Multiple transformations can lead to slightly different results
#' ## due to edge cases:
#' ## 4 triangles found in the polygon
#' fm_contains(
#'   obj,
#'   fm_transform(mesh, crs = fm_crs("mollweide_norm"))
#' )
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
  if (!inherits(y, "fm_mesh_2d")) {
    stop(paste0(
      "'y' must be an 'fm_mesh_2d' object, not '",
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
  if (fm_manifold(y, "S2")) {
    points <- points / rowSums(points^2)^0.5
  }
  ## Convert to sf points
  ## Extract coordinate system information
  if (fm_manifold(y, "S2")) {
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
