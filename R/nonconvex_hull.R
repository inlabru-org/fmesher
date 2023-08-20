#' @include deprecated.R

# fm_nonconvex_hull_inla ####

#' @title Contour segment
#'
#' @description
#' Helper from legacy `INLA::inla.contour.segment()`
#' @returns An `fm_segm` object
#' @export
#' @keywords internal
#' @family nonconvex inla legacy support
#' @examples
#' fm_segm_contour_helper(z = matrix(1:16, 4, 4))
#'
fm_segm_contour_helper <- function(x = seq(0, 1, length.out = nrow(z)),
                                   y = seq(0, 1, length.out = ncol(z)),
                                   z, nlevels = 10,
                                   levels = pretty(range(z, na.rm = TRUE), nlevels),
                                   groups = seq_len(length(levels)),
                                   positive = TRUE,
                                   eps = NULL,
                                   eps_rel = NULL,
                                   crs = NULL) {
  ## Input checking from contourLines:
  if (missing(z)) {
    if (!missing(x)) {
      if (is.list(x)) {
        z <- x$z
        y <- x$y
        x <- x$x
      } else {
        z <- x
        x <- seq.int(0, 1, length.out = nrow(z))
      }
    } else {
      stop("no 'z' matrix specified")
    }
  } else if (is.list(x)) {
    y <- x$y
    x <- x$x
  }

  if (is.null(eps)) {
    eps <- min(c(min(diff(x)), min(diff(y)))) / 8
  }
  ## End of input checking.

  ## Get contour pieces
  curves <- grDevices::contourLines(x, y, z, levels = levels)

  ## Make a mesh for easy gradient interpolation:
  latt <- fm_lattice_2d(x, y)
  mesh <-
    fm_rcdt_2d(
      lattice = latt,
      boundary = latt$segm,
      extend = (list(
        n = 3,
        offset = (max(
          diff(range(x)),
          diff(range(y))
        ) * 0.1)
      ))
    )
  ## Map function values to mesh indexing:
  zz <- rep(0, mesh$n)
  ## 2020-07-15:
  ## Bug in fmesher may lead to some lattice points missing. Ignore those.
  lattice_points_ok <- !is.na(mesh$idx$lattice)
  zz[mesh$idx$lattice[lattice_points_ok]] <- as.vector(z)[lattice_points_ok]

  ## Mapping from level to group value:
  level2grp <- function(level) {
    if (length(groups) == 1) {
      return(groups)
    }
    for (k in seq_along(groups)) {
      if (levels[k] == level) {
        return(groups[k])
      }
    }
    return(0)
  }

  ## Join all contour pieces into a single mesh.segment
  ## Different levels can later be identified via the grp indices.
  loc <- matrix(0, 0, 2)
  idx <- matrix(0, 0, 2)
  grp <- c()
  for (k in seq_len(length(curves))) {
    curve.loc <- cbind(curves[[k]]$x, curves[[k]]$y)
    curve.n <- nrow(curve.loc)

    ## Extract the rotated gradients along the curve
    curve.mid <- (curve.loc[1:(curve.n - 1), , drop = FALSE] +
      curve.loc[2:curve.n, , drop = FALSE]) / 2
    A <- fm_basis(mesh, loc = curve.mid, derivatives = TRUE)
    ## Gradients rotated 90 degrees CW, i.e. to the direction
    ## of CCW curves around positive excursions:
    grid.diff <- cbind(A$dy %*% zz, -A$dx %*% zz)

    ## Determine the CCW/CW orientation
    curve.diff <- diff(curve.loc)
    ccw <- (sum(curve.diff * grid.diff) >= 0) ## True if in CCW direction
    if ((ccw && positive) || (!ccw && !positive)) {
      curve.idx <- seq_len(curve.n)
    } else {
      curve.idx <- rev(seq_len(curve.n))
    }

    ## Filter short line segments:
    curve.idx <-
      fm_simplify_helper(curve.loc,
        curve.idx,
        eps = eps,
        eps_rel = eps_rel
      )

    ## Reorder, making sure any unused points are removed:
    curve.loc <- curve.loc[curve.idx, , drop = FALSE]
    curve.n <- nrow(curve.loc)
    curve.idx <- cbind(seq_len(curve.n - 1L), seq_len(curve.n - 1L) + 1L)

    ## Check if the curve is closed, and adjust if it is:
    if (max(abs(curve.loc[1, , drop = FALSE] - curve.loc[curve.n, , drop = FALSE])) < 1e-12) {
      curve.loc <- curve.loc[-curve.n, , drop = FALSE]
      curve.n <- nrow(curve.loc)
      curve.idx <-
        cbind(seq_len(curve.n), c(seq_len(curve.n - 1L) + 1L, 1L))
    }

    ## Add the curve:
    offset <- nrow(loc)
    loc <- rbind(loc, curve.loc)
    idx <- rbind(idx, curve.idx + offset)
    grp <- c(grp, rep(level2grp(curves[[k]]$level), curve.n - 1L))
  }

  fm_segm(loc = loc, idx = idx, grp = grp, is.bnd = FALSE, crs = crs)
}

#' @title Non-convex hull computation
#' @description Legacy method for `INLA::inla.nonconvex.hull()`
#' @seealso [fm_nonconvex_hull()]
#' @param resolution The internal computation resolution.  A warning will be
#' issued when this needs to be increased for higher accuracy, with the
#' required resolution stated.
#' @param eps,eps_rel The polygonal curve simplification tolerances used for
#' simplifying the resulting boundary curve.  See [fm_simplify_helper()] for
#' details.
#' @param \dots Unused.
#' @inheritParams fm_nonconvex_hull
#' @returns `fm_nonconvex_hull_inla()` returns an `fm_segm`/`inla.mesh.segment`
#' object, for compatibility with `inla.nonconvex.hull()`.
#' @export
#' @family nonconvex inla legacy support
#' @examples
#' fm_nonconvex_hull_inla(cbind(0, 0), convex = 1)
#'
fm_nonconvex_hull_inla <- function(x,
                                   convex = -0.15,
                                   concave = convex,
                                   resolution = 40,
                                   eps = NULL,
                                   eps_rel = NULL,
                                   crs = NULL,
                                   ...) {
  stopifnot(!is.null(x))
  if (inherits(x, c("SpatialPoints", "SpatialPointsDataFrame"))) {
    x <- fm_transform(
      sp::coordinates(x),
      crs0 = fm_crs(x),
      crs = fm_crs(crs),
      passthrough = TRUE
    )
  } else if (inherits(x, c("sf", "sfc"))) {
    x <- fm_transform(
      sf::st_coordinates(x),
      crs0 = fm_crs(x),
      crs = fm_crs(x),
      passthrough = TRUE
    )
  }

  if (length(resolution) == 1) {
    resolution <- rep(resolution, 2)
  }
  lim <- rbind(range(x[, 1]), range(x[, 2]))

  approx.diam <- max(diff(lim[1, ]), diff(lim[2, ]))
  if (convex < 0) {
    convex <- -convex * approx.diam
  }
  if (concave < 0) {
    concave <- -concave * approx.diam
  }
  if (concave == 0) {
    return(fm_nonconvex_hull_inla_basic(x, convex, resolution, eps,
      crs = crs
    ))
  }

  ex <- convex + concave
  domain <- c(diff(lim[1, ]), diff(lim[2, ])) + 2 * ex
  dif <- domain / (resolution - 1)
  if (max(dif) > min(convex, concave)) {
    req.res <- ceiling(domain / min(convex, concave) + 1)
    warning(paste("Resolution (",
      paste(resolution, collapse = ","),
      ") too small for convex/concave radius (",
      convex, ",", concave,
      ").\n",
      "Resolution >=(",
      paste(req.res, collapse = ","),
      ") required for more accurate results.",
      sep = ""
    ))
  }
  ax <-
    list(
      seq(lim[1, 1] - ex, lim[1, 2] + ex, length.out = resolution[1]),
      seq(lim[2, 1] - ex, lim[2, 2] + ex, length.out = resolution[2])
    )
  xy <- as.matrix(expand.grid(ax[[1]], ax[[2]]))

  fm_require_stop("splancs")
  z <- (matrix(
    splancs::nndistF(x, xy),
    resolution[1], resolution[2]
  ))
  segm.dilation <-
    fm_segm_contour_helper(
      ax[[1]], ax[[2]], z,
      levels = c(convex + concave),
      positive = TRUE,
      eps = 0
    ) ## Don't simplify curve at this stage
  mesh.dilation <-
    fm_rcdt_2d(
      loc = xy,
      boundary = segm.dilation,
      extend = (list(
        n = 3,
        offset = (max(
          diff(ax[[1]]),
          diff(ax[[2]])
        ) * 0.1)
      ))
    )

  z <- (matrix(
    splancs::nndistF(mesh.dilation$loc, xy),
    resolution[1], resolution[2]
  ))
  segm.closing <-
    fm_segm_contour_helper(
      ax[[1]], ax[[2]], z,
      levels = c(concave),
      positive = TRUE,
      eps = eps,
      eps_rel = eps_rel
    )

  segm.closing$crs <- crs

  fm_as_segm(segm.closing)
}



#' @details Requires `splancs::nndistF()`
#' @export
#' @describeIn fm_nonconvex_hull_inla Special method for `convex = 0`.
## Based on an idea from Elias Teixeira Krainski
fm_nonconvex_hull_inla_basic <- function(x, convex = -0.15, resolution = 40,
                                         eps = NULL, crs = NULL) {
  stopifnot(!is.null(x))
  if (inherits(x, c("SpatialPoints", "SpatialPointsDataFrame"))) {
    x <- fm_transform(
      sp::coordinates(x),
      crs0 = fm_crs(x),
      crs = fm_crs(crs),
      passthrough = TRUE
    )
  } else if (inherits(x, c("sf", "sfc"))) {
    x <- fm_transform(
      sf::st_coordinates(x),
      crs0 = fm_crs(x),
      crs = fm_crs(x),
      passthrough = TRUE
    )
  }

  if (length(convex) == 1) {
    convex <- rep(convex, 2)
  }
  if (length(resolution) == 1) {
    resolution <- rep(resolution, 2)
  }

  lim <- rbind(range(x[, 1]), range(x[, 2]))
  ex <- convex
  if (convex[1] < 0) {
    ex[1] <- -convex[1] * diff(lim[1, ])
  }
  if (convex[2] < 0) {
    ex[2] <- -convex[2] * diff(lim[2, ])
  }

  domain <- c(diff(lim[1, ]), diff(lim[2, ])) + 2 * ex
  dif <- domain / (resolution - 1)
  if (any(dif > min(convex))) {
    req.res <- ceiling(domain / convex + 1)
    warning(paste("Resolution (",
      paste(resolution, collapse = ","),
      ") too small for convex (",
      paste(convex, collapse = ","),
      ").\n",
      "Resolution >=(",
      paste(req.res, collapse = ","),
      ") required for more accurate results.",
      sep = ""
    ))
  }

  ax <- list(
    seq(lim[1, 1] - ex[1], lim[1, 2] + ex[1], length.out = resolution[1]),
    seq(lim[2, 1] - ex[2], lim[2, 2] + ex[2], length.out = resolution[2])
  )
  xy <- as.matrix(expand.grid(ax[[1]], ax[[2]]))
  tr <- diag(c(1 / ex[1], 1 / ex[2]), nrow = 2, ncol = 2)

  fm_require_stop("splancs")
  z <- matrix(
    splancs::nndistF(x %*% tr, xy %*% tr),
    resolution[1], resolution[2]
  )
  segm <- fm_segm_contour_helper(
    ax[[1]],
    ax[[2]],
    z,
    levels = c(1),
    positive = FALSE,
    eps = eps,
    crs = crs
  )
  return(segm)
}




# fm_nonconvex_hull ####

#' @title Compute an extension of a spatial object
#'
#' @description
#' Constructs a potentially nonconvex extension of a spatial object by
#' performing dilation by `convex + concave` followed by
#' erosion by `concave`. This is equivalent to dilation by `convex` followed
#' by closing (dilation + erosion) by `concave`.
#'
#' @describeIn fm_nonconvex_hull Basic nonconvex hull method.
#'
#' @details
#' Morphological dilation by `convex`, followed by closing by
#' `concave`, with minimum concave curvature radius `concave`.  If
#' the dilated set has no gaps of width between \deqn{2 \text{convex} (\sqrt{1+2
#' \text{concave}/\text{convex}} - 1)}{2*convex*(sqrt(1+2*concave/convex) - 1)}
#' and \eqn{2\text{concave}}{2*concave}, then the minimum convex curvature radius is
#' `convex`.
#'
#' The implementation is based on the identity \deqn{\text{dilation}(a) \&
#' \text{closing}(b) = \text{dilation}(a+b) \& \text{erosion}(b)}{
#' dilation(a) & closing(b) = dilation(a+b) & erosion(b)} where all operations
#' are with respect to disks with the specified radii.
#'
#' @param x A spatial object
#' @param x A spatial object
#' @param convex numeric vector; How much to extend
#' @param concave numeric vector; The minimum allowed reentrant curvature. Default equal to `convex`
#' @param preserveTopology logical; argument to `sf::st_simplify()`
#' @param dTolerance If not zero, controls the `dTolerance` argument to `sf::st_simplify()`.
#' The default is `pmin(convex, concave) / 40`, chosen to
#' give approximately 4 or more subsegments per circular quadrant.
#' @param crs Options crs object for the resulting polygon
#' @param ... Arguments passed on to the [fm_nonconvex_hull()] sub-methods
#' @details When `convex`, `concave`, or `dTolerance` are negative,
#' `fm_diameter * abs(...)` is used instead.
#' @returns `fm_nonconvex_hull()` returns an extended object as an `sfc`
#' polygon object (regardless of the `x` class).
#' @references Gonzalez and Woods (1992), Digital Image Processing
#' @seealso [fm_nonconvex_hull_inla()]
#' @export
#' @examples
#' inp <- matrix(rnorm(20), 10, 2)
#' out <- fm_nonconvex_hull(inp, convex = 1)
#' plot(out)
#' points(inp, pch = 20)
fm_nonconvex_hull <- function(x, ...) {
  UseMethod("fm_nonconvex_hull")
}

#' @rdname fm_nonconvex_hull
#' @details Differs from `sf::st_buffer(x, convex)` followed by
#' `sf::st_concave_hull()` (available from GEOS 3.11)
#' in how the amount of allowed concavity is controlled.
#' @export
fm_nonconvex_hull.sfc <- function(x,
                                  convex = -0.15,
                                  concave = convex,
                                  preserveTopology = TRUE,
                                  dTolerance = NULL,
                                  crs = fm_crs(x),
                                  ...) {
  diameter_bound <- fm_diameter(x)
  scale_fun <- function(val) {
    if (val < 0) {
      val <- diameter_bound * abs(val)
    }
    val
  }
  convex <- scale_fun(convex)
  concave <- scale_fun(concave)
  if (is.null(dTolerance)) {
    dTolerance <- min(convex, concave) / 40
  } else {
    dTolerance <- scale_fun(dTolerance)
  }

  nQuadSegs <- 64
  y <- sf::st_buffer(x, dist = convex + concave, nQuadSegs = nQuadSegs)
  y <- sf::st_union(y)
  if (concave > 0) {
    if (sf::sf_use_s2() && isTRUE(sf::st_is_longlat(x))) {
      # s2 gives empty result for negative buffers
      # Use bounding set trick to get around this
      y_box <- sf::st_buffer(y, dist = 100000, nQuadSegs = nQuadSegs)
      y_box <- sf::st_union(y_box)
      y_inverse <- sf::st_difference(y_box, y)
      y_inverse_expanded <- sf::st_buffer(y_inverse, concave, nQuadSegs = nQuadSegs)
      y_inverse_expanded <- sf::st_union(y_inverse_expanded)
      y <- sf::st_difference(y, y_inverse_expanded)
    } else {
      y <- sf::st_buffer(y, dist = -concave, nQuadSegs = nQuadSegs)
    }
  }

  if (dTolerance > 0) {
    y <- sf::st_simplify(y,
      preserveTopology = preserveTopology,
      dTolerance = dTolerance
    )
  }
  y <- sf::st_union(y)
  y <- sf::st_sfc(y, crs = fm_crs(x))
  if (!fm_crs_is_identical(fm_crs(y), crs)) {
    y <- fm_transform(y, crs = crs)
  }
  y
}


#' @describeIn fm_nonconvex_hull
#' Constructs a potentially nonconvex extension of a spatial object by
#' performing dilation by `convex + concave` followed by
#' erosion by `concave`. This is equivalent to dilation by `convex` followed
#' by closing (dilation + erosion) by `concave`.
#'
#' @returns `fm_extensions()` returns a list of `sfc` objects.
#' @export
#' @examples
#' if (TRUE) {
#'   inp <- sf::st_as_sf(as.data.frame(matrix(1:6, 3, 2)), coords = 1:2)
#'   bnd <- fm_extensions(inp, convex = c(0.75, 2))
#'   plot(fm_mesh_2d(boundary = bnd, max.edge = c(0.25, 1)), asp = 1)
#' }
fm_extensions <- function(x,
                          convex = -0.15,
                          concave = convex,
                          dTolerance = NULL,
                          ...) {
  if (any(convex < 0) || any(concave < 0) || any(dTolerance < 0)) {
    diameter_bound <- fm_diameter(x)
  }
  len <- max(length(convex), length(concave), length(dTolerance))
  scale_fun <- function(val) {
    if (any(val < 0)) {
      val[val < 0] <- diameter_bound * abs(val[val < 0])
    }
    if (length(val) < len) {
      val <- c(val, rep(val[length(val)], len - length(val)))
    }
    val
  }
  convex <- scale_fun(convex)
  concave <- scale_fun(concave)
  if (is.null(dTolerance)) {
    dTolerance <- pmin(convex, concave) / 40
  } else {
    dTolerance <- scale_fun(dTolerance)
  }

  y <- lapply(
    seq_along(convex),
    function(k) {
      fm_nonconvex_hull(
        x,
        convex = convex[k],
        concave = concave[k],
        dTolerance = dTolerance[k],
        ...
      )
    }
  )
  y
}



#' @rdname fm_nonconvex_hull
#' @export
fm_nonconvex_hull.matrix <- function(x, ...) {
  fm_nonconvex_hull.sfc(sf::st_multipoint(x), ...)
}

#' @rdname fm_nonconvex_hull
#' @export
fm_nonconvex_hull.sf <- function(x, ...) {
  fm_nonconvex_hull.sfc(sf::st_geometry(x), ...)
}

#' @rdname fm_nonconvex_hull
#' @export
fm_nonconvex_hull.Spatial <- function(x, ...) {
  fm_nonconvex_hull.sfc(sf::st_as_sfc(x), ...)
}

#' @rdname fm_nonconvex_hull
#' @export
fm_nonconvex_hull.sfg <- function(x, ...) {
  fm_nonconvex_hull.sfc(sf::st_sfc(x), ...)
}
