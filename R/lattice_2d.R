#' @include deprecated.R

# fm_lattice_2d ####

#' @title Special coordinate mappings for `fm_mesh_2d` projections.
#'
#' @description
#' Calculates coordinate mappings for spherical `fm_mesh_2d` projections.
#' This is an internal function not intended for general use.
#'
#' @keywords internal
#' @param loc Coordinates to be mapped.
#' @param projection The projection type.
#' @param inverse If `TRUE`, `loc` are map coordinates and
#' coordinates in the spherical domain are calculated.  If `FALSE`, `loc`
#' are coordinates in the spherical domain and the forward map projection is
#' calculated. Default: `TRUE`
#' @return For `fm_mesh_2d_map_lim`, a list:
#' \item{xlim }{X axis limits in the map domain}
#' \item{ylim }{Y axis limits in the map domain}
#' No attempt is
#' made to find minimal limits for partial spherical domains.
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @seealso [fm_evaluator()]
#' @export
#' @examples
#' (loc <- fm_mesh_2d_map(cbind(20, 10), "longlat"))
#' fm_mesh_2d_map(loc, "longlat", inverse = FALSE)
#'
fm_mesh_2d_map <- function(loc,
                           projection =
                             c("default", "longlat", "longsinlat", "mollweide"),
                           inverse = TRUE) {
  projection <- match.arg(projection)
  if (identical(projection, "default")) {
    return(loc)
  } else if (identical(projection, "longlat")) {
    if (inverse) {
      proj <-
        cbind(
          cos(loc[, 1] * pi / 180) * cos(loc[, 2] * pi / 180),
          sin(loc[, 1] * pi / 180) * cos(loc[, 2] * pi / 180),
          sin(loc[, 2] * pi / 180)
        )
    } else {
      proj <-
        cbind(
          atan2(loc[, 2], loc[, 1]) * 180 / pi,
          asin(pmax(-1, pmin(+1, loc[, 3]))) * 180 / pi
        )
    }
  } else if (identical(projection, "longsinlat")) {
    if (inverse) {
      coslat <- sqrt(pmax(0, 1 - loc[, 2]^2))
      proj <-
        cbind(
          cos(loc[, 1] * pi / 180) * coslat,
          sin(loc[, 1] * pi / 180) * coslat,
          loc[, 2]
        )
    } else {
      proj <-
        cbind(
          atan2(loc[, 2], loc[, 1]) * 180 / pi,
          loc[, 3]
        )
    }
  } else if (identical(projection, "mollweide")) {
    if (inverse) {
      ok <- ((loc[, 1]^2 + 4 * loc[, 2]^2) <= 4)
      cos.theta <- sqrt(pmax(0, 1 - loc[ok, 2]^2))
      theta <- atan2(loc[ok, 2], cos.theta)
      sin.lat <- (2 * theta + sin(2 * theta)) / pi
      cos.lat <- sqrt(pmax(0, 1 - sin.lat^2))
      lon <- loc[ok, 1] * pi / 2 / (cos.theta + (cos.theta == 0))
      lon[cos.theta == 0] <- pi / 2 * sign(theta[cos.theta == 0])
      proj <- matrix(NA, nrow(loc), 3)
      proj[ok, ] <- cbind(cos(lon) * cos.lat, sin(lon) * cos.lat, sin.lat)
    } else {
      lon <- atan2(loc[, 2], loc[, 1])
      z <- pmin(1, pmax(-1, loc[, 3]))
      sin.theta <- z
      cos.theta <- sqrt(pmax(0, 1 - sin.theta^2))
      ## NR-solver for sin.theta.
      ## Typically finishes after at most 7 iterations.
      ## When cos.theta=0, sin.theta is already correct, +/- 1.
      nook <- (cos.theta > 0)
      for (k in 1:20) {
        if (any(nook)) {
          delta <-
            (atan2(sin.theta[nook], cos.theta[nook]) +
              sin.theta[nook] * cos.theta[nook] - pi / 2 * z[nook]) /
              (2 * cos.theta[nook])
          sin.theta[nook] <- sin.theta[nook] - delta
          cos.theta[nook] <- sqrt(1 - sin.theta[nook]^2)
          nook[nook] <- (abs(delta) > 1e-14)
        }
      }
      proj <- cbind(2 * lon / pi * cos.theta, sin.theta)
    }
  } else {
    stop(paste("Unknown projection '", projection, "'.", sep = ""))
  }
  return(proj)
}


#' @export
#' @describeIn fm_mesh_2d_map Projection extent limit calculations
fm_mesh_2d_map_lim <- function(loc = NULL,
                               projection =
                                 c("default", "longlat", "longsinlat", "mollweide")) {
  projection <- match.arg(projection)
  if (identical(projection, "default")) {
    if (is.null(loc)) {
      lim <- list(xlim = c(0, 1), ylim = c(0, 1))
    } else {
      lim <-
        list(
          xlim = range(loc[, 1], na.rm = TRUE),
          ylim = range(loc[, 2], na.rm = TRUE)
        )
    }
  } else if (identical(projection, "longlat")) {
    lim <- list(xlim = c(-180, 180), ylim = c(-90, 90))
  } else if (identical(projection, "longsinlat")) {
    lim <- list(xlim = c(-180, 180), ylim = c(-1, 1))
  } else if (identical(projection, "mollweide")) {
    lim <- list(xlim = c(-2, 2), ylim = c(-1, 1))
  } else {
    stop(paste("Unknown projection '", projection, "'.", sep = ""))
  }
  return(lim)
}




#' @title Make a lattice object
#' @export
#' @param ... Passed on to submethods
#' @family object creation and conversion
fm_lattice_2d <- function(...) {
  UseMethod("fm_lattice_2d")
}

#' Lattice grids for inla.mesh
#'
#' Construct a lattice grid for [fm_mesh_2d()]
#'
#' @param x vector or grid matrix of x-values. Vector values are sorted before use.
#' Matrix input is assumed to be a grid of x-values with the same ordering convention
#' of `as.vector(x)` as `rep(x, times = dims[2])` for vector input.
#' @param y vector of grid matrix of y-values. Vector values are sorted before use.
#' Matrix input is assumed to be a grid of y-values with the same ordering convention
#' of `as.vector(y)` as `rep(y, each = dims[1])` for vector input.
#' @param z if x is a matrix, a grid matrix of z-values, with the same ordering as `x`
#' and `y`. If `x` is a vector, `z` is ignored.
#' @param dims the size of the grid, length 2 vector
#' @param units One of `c("default", "longlat", "longsinlat", "mollweide")`
#' or NULL (equivalent to `"default"`).
#' @param crs An optional `fm_crs`, `sf::st_crs`, or `sp::CRS` object
#' @return An `fm_lattice_2d` object with elements
#' \describe{
#' \item{dims}{integer vector}
#' \item{x}{x-values for original vector input}
#' \item{y}{y-values for original vector input}
#' \item{loc}{matrix of `(x, y)` values or `(x, y, z)` values. May be altered by [fm_transform()]}
#' \item{segm}{`fm_segm` object}
#' \item{crs}{`fm_crs` object or `NULL`}
#' }
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @seealso [fm_mesh_2d()]
#' @examples
#' lattice <- fm_lattice_2d(
#'   seq(0, 1, length.out = 17),
#'   seq(0, 1, length.out = 10)
#' )
#'
#' ## Use the lattice "as-is", without refinement:
#' mesh <- fm_rcdt_2d_inla(lattice = lattice, boundary = lattice$segm)
#' mesh <- fm_rcdt_2d_inla(lattice = lattice, extend = FALSE)
# plot(mesh)
#'
#' ## Refine the triangulation, with limits on triangle angles and edges:
#' mesh <- fm_rcdt_2d(
#'   lattice = lattice,
#'   refine = list(max.edge = 0.08),
#'   extend = FALSE
#' )
# plot(mesh)
#'
#' ## Add an extension around the lattice, but maintain the lattice edges:
#' mesh <- fm_rcdt_2d(
#'   lattice = lattice,
#'   refine = list(max.edge = 0.08),
#'   interior = lattice$segm
#' )
# plot(mesh)
#'
#' ## Only add extension:
#' mesh <- fm_rcdt_2d(lattice = lattice, refine = list(max.edge = 0.08))
# plot(mesh)
#' @rdname fm_lattice_2d
#' @export
fm_lattice_2d.default <- function(
    x = seq(0, 1, length.out = 2),
    y = seq(0, 1, length.out = 2),
    z = NULL,
    dims =
      if (is.matrix(x)) {
        dim(x)
      } else {
        c(length(x), length(y))
      },
    units = NULL,
    crs = NULL,
    ...) {
  if (is.null(crs)) {
    units <- match.arg(units, c("default", "longlat", "longsinlat", "mollweide"))

    lim <- fm_mesh_2d_map_lim(projection = units)
    xlim <- lim$xlim
    ylim <- lim$ylim
  } else { ## !is.null(crs)
    if (!is.null(units)) {
      stop("Only one of 'units' and 'crs' can be non-null.")
    }
    bounds <- fm_crs_bounds(crs, warn.unknown = TRUE)
    xlim <- bounds$xlim
    ylim <- bounds$ylim
  }

  if (missing(x) && !missing(dims)) {
    x <- seq(xlim[1], xlim[2], length.out = dims[1])
  }
  if (missing(y) && !missing(dims)) {
    y <- seq(ylim[1], ylim[2], length.out = dims[2])
  }
  dims <- as.integer(dims)

  if (is.matrix(x)) {
    if (!identical(dims, dim(x)) ||
      !identical(dims, dim(y)) ||
      (is.matrix(z) && !identical(dims, dim(z)))) {
      stop("The size of matrices 'x', 'y', and 'z' must match 'dims'.")
    }
    loc <- cbind(as.vector(x), as.vector(y), as.vector(z))
    x <- NULL
    y <- NULL
  } else {
    if (!identical(dims[1], length(x)) ||
      !identical(dims[2], length(y))) {
      stop(paste("The lengths of vectors 'x' and 'y' (",
        length(x), ",", length(y),
        ") must match 'dims' (", dims[1], ",", dims[2], ").",
        sep = ""
      ))
    }

    # Ensure correct point ordering
    x <- sort(unique(x))
    y <- sort(unique(y))

    # Expand coordinates
    loc <- (cbind(
      rep(x, times = dims[2]),
      rep(y, each = dims[1])
    ))
  }
  if (!is.double(loc)) {
    storage.mode(loc) <- "double"
  }

  if (is.null(crs)) {
    loc <- fm_mesh_2d_map(loc = loc, projection = units, inverse = TRUE)
  }

  ## Construct lattice boundary
  segm.idx <- (c(
    1:(dims[1] - 1),
    dims[1] * (1:(dims[2] - 1)),
    dims[1] * dims[2] - (0:(dims[1] - 2)),
    dims[1] * ((dims[2] - 1):1) + 1
  ))
  segm.grp <- (c(
    rep(1L, dims[1] - 1),
    rep(2L, dims[2] - 1),
    rep(3L, dims[1] - 1),
    rep(4L, dims[2] - 1)
  ))

  segm <- fm_segm(
    loc = loc[segm.idx, , drop = FALSE],
    grp = segm.grp,
    is.bnd = TRUE,
    crs = crs
  )

  lattice <- list(dims = dims, x = x, y = y, loc = loc, segm = segm, crs = crs)
  class(lattice) <- c("fm_lattice_2d", "inla.mesh.lattice")
  return(lattice)
}

#' @title Convert objects to `fm_lattice_2d`
#' @describeIn fm_as_lattice_2d Convert an object to `fm_lattice_2d`.
#' @param x Object to be converted.
#' @param ... Arguments passed on to submethods
#' @returns An `fm_lattice_2d` or `fm_lattice_2d_list` object
#' @export
#' @family object creation and conversion
#' @export
#' @examples
#' str(fm_as_lattice_2d_list(list(fm_lattice_2d(), fm_lattice_2d())))
#'
fm_as_lattice_2d <- function(...) {
  UseMethod("fm_as_lattice_2d")
}
#' @describeIn fm_as_lattice_2d Convert each element of a list
#' @export
fm_as_lattice_2d_list <- function(x, ...) {
  fm_as_list(x, ..., .class_stub = "lattice_2d")
}
#' @rdname fm_as_lattice_2d
#' @param x Object to be converted
#' @export
fm_as_lattice_2d.fm_lattice_2d <- function(x, ...) {
  #  class(x) <- c("fm_lattice_2d", setdiff(class(x), "fm_lattice_2d"))
  x
}
#' @rdname fm_as_lattice_2d
#' @param x Object to be converted
#' @export
#' @method fm_as_lattice_2d inla.mesh.lattice
fm_as_lattice_2d.inla.mesh.lattice <- function(x, ...) {
  x[["crs"]] <- fm_crs(x[["crs"]])
  x[["segm"]] <- fm_as_fm(x[["segm"]])
  class(x) <- c("fm_lattice_2d", class(x))
  x
}
