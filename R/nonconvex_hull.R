#' @include deprecated.R

#' @title Contour segment
#'
#' @description
#' Helper from legacy `INLA::inla.contour.segment()`
#' @param x,y,z The `x` and `y` coordinates of the grid, and the `z` values
#' @param nlevels Number of contour levels
#' @param levels The contour levels. If `NULL`,
#'   `pretty(range(z, na.rm = TRUE), nlevels)` is used.
#' @param groups The group values for each contour level.
#'   If `NULL`, `seq_len(length(levels))` is used.
#' @param positive Logical; if `TRUE`, the contour lines are made to be
#' CCW around positive excursions
#' @param eps,eps_rel Polygonal curve simplification tolerances
#' @param crs A coordinate reference system
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
                                   levels = NULL,
                                   groups = NULL,
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

  if (is.null(levels)) {
    levels <- pretty(range(z, na.rm = TRUE), nlevels)
  }
  if (is.null(groups)) {
    groups <- seq_len(length(levels))
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
    0
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
    A <- fm_basis(mesh, loc = curve.mid, derivatives = TRUE, full = TRUE)
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
    if (max(abs(curve.loc[1, , drop = FALSE] -
      curve.loc[curve.n, , drop = FALSE])) < 1e-12) {
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



# fm_nonconvex_hull ####

#' @title Compute an extension of a spatial object
#'
#' @description
#' Constructs a potentially nonconvex extension of a spatial object by
#' performing dilation by `convex + concave` followed by
#' erosion by `concave`. This is equivalent to dilation by `convex` followed
#' by closing (dilation + erosion) by `concave`.
#'
#' @details
#' Morphological dilation by `convex`, followed by closing by
#' `concave`, with minimum concave curvature radius `concave`.  If
#' the dilated set has no gaps of width between \deqn{2 \textrm{convex}
#' (\sqrt{1+2\textrm{concave}/\textrm{convex}} - 1)
#' }{2*convex*(sqrt(1+2*concave/convex) - 1)}
#' and \eqn{2\textrm{concave}}{2*concave}, then the minimum convex curvature
#' radius is `convex`.
#'
#' The implementation is based on the identity \deqn{\textrm{dilation}(a) \&
#' \textrm{closing}(b) = \textrm{dilation}(a+b) \& \textrm{erosion}(b)}{
#' dilation(a) & closing(b) = dilation(a+b) & erosion(b)} where all operations
#' are with respect to disks with the specified radii.
#'
#' @param x A spatial object
#' @param format character specifying the output format; "sf" (default) or "fm"
#' @param method character specifying the construction method; "fm" (default)
#'   or "sf"
#' @param convex numeric vector; How much to extend
#' @param concave numeric vector; The minimum allowed reentrant curvature.
#'   Default equal to `convex`
#' @param crs Optional crs object for the resulting polygon. Default is
#'   `fm_crs(x)`
#' @param ... Arguments passed on to the [fm_nonconvex_hull()] sub-methods
#' @details When `convex`, `concave`, or `dTolerance` are negative,
#' `fm_diameter * abs(...)` is used instead.
#' @returns `fm_nonconvex_hull()` returns an extended object as an `sfc` polygon
#'   object (if `format = "sf"`) or an [fm_segm] object (if `format = "fm")
#' @references Gonzalez and Woods (1992), Digital Image Processing
#' @seealso [fm_nonconvex_hull_inla()]
#' @export
#' @inheritSection fm_mesh_2d INLA compatibility
#' @examples
#' inp <- matrix(rnorm(20), 10, 2)
#' out <- fm_nonconvex_hull(inp, convex = 1, method = "sf")
#' plot(out)
#' points(inp, pch = 20)
#'
#' out <- fm_nonconvex_hull(inp, convex = 1, method = "fm", format = "fm")
#' lines(out, col = 2, add = TRUE)
fm_nonconvex_hull <- function(x, ..., format = "sf", method = "fm") {
  if (match.arg(method, c("fm", "sf")) == "fm") {
    if (!requireNamespace("splancs", quietly = TRUE)) {
      stop("Package 'splancs' is required for method='fm'. Please install it.")
    }
  }
  UseMethod("fm_nonconvex_hull")
}


#' @describeIn fm_nonconvex_hull
#' Constructs a potentially nonconvex extension of a spatial object by
#' performing dilation by `convex + concave` followed by
#' erosion by `concave`. This is equivalent to dilation by `convex` followed
#' by closing (dilation + erosion) by `concave`.
#'
#' The `...` arguments are passed on to `fm_nonconvex_hull_fm()`
#' or `fm_nonconvex_hull_sf()`, depending on the `method` argument.
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
                          ...,
                          format = "sf",
                          method = "fm") {
  if (any(convex < 0) || any(concave < 0)) {
    diameter_bound <- max(fm_diameter(x))
  }
  len <- max(length(convex), length(concave))
  if ("dTolerance" %in% names(list(...))) {
    dTolerance <- list(...)$dTolerance
    len <- max(len, length(dTolerance))
  } else {
    dTolerance <- NULL
  }
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

  args <- list(...)
  args[["dTolerance"]] <- NULL
  y <- lapply(
    seq_along(convex),
    function(k) {
      do.call(
        fm_nonconvex_hull,
        c(
          list(
            x = x,
            convex = convex[k],
            concave = concave[k],
            dTolerance = dTolerance[k],
            format = format,
            method = method
          ),
          args
        )
      )
    }
  )

  if (format == "fm") {
    y <- fm_as_segm_list(y)
  }

  y
}






#' @describeIn fm_nonconvex_hull `fmesher` method for `fm_nonconvex_hull()`,
#'   which uses the `splancs::nndistF()` function to compute nearest-neighbour
#'   distances.
#' @param resolution integer; The internal computation resolution.  A warning
#'   will be issued when this needs to be increased for higher accuracy, with
#'   the required resolution stated. For `method="fm"` only.
#' @param eps,eps_rel The polygonal curve simplification tolerances used for
#'   simplifying the resulting boundary curve.  See [fm_simplify_helper()] for
#'   details. For `method="fm"` only.
#' @export
fm_nonconvex_hull_fm <- function(x,
                                 convex = -0.15,
                                 concave = convex,
                                 resolution = 40,
                                 eps = NULL,
                                 eps_rel = NULL,
                                 crs = fm_crs(x),
                                 ...) {
  stopifnot(!is.null(x))
  diameter_bound <- fm_diameter(x)
  scale_fun <- function(val) {
    if (val < 0) {
      val <- diameter_bound * abs(val)
    }
    val
  }
  convex <- scale_fun(convex)
  concave <- scale_fun(concave)

  if (inherits(x, c("SpatialPoints", "SpatialPointsDataFrame"))) {
    fm_safe_sp(force = TRUE)
    x <- fm_transform(
      sp::coordinates(x),
      crs0 = fm_crs(x),
      crs = fm_crs(crs),
      passthrough = TRUE
    )
    x <- x[, 1:2, drop = FALSE]
  } else if (inherits(x, c("sf", "sfc"))) {
    x <- fm_transform(
      x,
      crs0 = fm_crs(x),
      crs = fm_crs(crs),
      passthrough = TRUE
    )

    z <- sf::st_coordinates(x)
    z <- z[, intersect(colnames(z), c("X", "Y")), drop = FALSE]

    # For both polygons and linestrings, add points along the lines.
    # If polygon, add interior points
    if (inherits(x, c("sfc_POLYGON", "sfc_MULTIPOLYGON"))) {
      x_interior <- fm_hexagon_lattice(x, edge_len = convex)

      z_int <- sf::st_coordinates(x_interior)
      if (NROW(z_int) > 0) {
        z_int <- z_int[, intersect(colnames(z_int), c("X", "Y")), drop = FALSE]

        z <- rbind(z, z_int)
      }
    }

    if (inherits(x, c(
      "sfc_POLYGON", "sfc_MULTIPOLYGON",
      "sfc_LINESTRING", "sfc_MULTILINESTRING"
    ))) {
      # Subdivide the polygon edges
      xx <- fm_as_segm(x)
      z_sub <- do.call(
        rbind,
        lapply(
          seq_len(NROW(xx$idx)),
          function(k) {
            N <- max(1, ceiling(sum((xx$loc[xx$idx[k, 1], ] -
              xx$loc[xx$idx[k, 2], ])^2)^0.5
              / (convex / 2)))
            v <- seq_len(N) / (N + 1)
            cbind(
              xx$loc[xx$idx[k, 1], 1] * (1 - v) +
                xx$loc[xx$idx[k, 2], 1] * v,
              xx$loc[xx$idx[k, 1], 2] * (1 - v) +
                xx$loc[xx$idx[k, 2], 2] * v
            )
          }
        )
      )
      z <- rbind(z, z_sub)
    }

    x <- z
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
    return(fm_nonconvex_hull_fm_basic(x, convex, resolution, eps, crs = crs))
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

  result <- fm_as_segm(segm.closing)
  fm_is_bnd(result) <- TRUE
  result
}


# Special [fm_nonconvex_hull_fm()] method for `concave = 0`.
# Called automatically by fm_nonconvex_hull_fm()
## Based on an idea from Elias Teixeira Krainski
fm_nonconvex_hull_fm_basic <- function(x, convex = -0.15, resolution = 40,
                                       eps = NULL, crs = fm_crs(x)) {
  stopifnot(!is.null(x))
  if (inherits(x, c("SpatialPoints", "SpatialPointsDataFrame"))) {
    fm_safe_sp(force = TRUE)
    x <- fm_transform(
      sp::coordinates(x),
      crs0 = fm_crs(x),
      crs = fm_crs(crs),
      passthrough = TRUE
    )
    x <- x[, 1:2, drop = FALSE]
  } else if (inherits(x, c("sf", "sfc"))) {
    x <- fm_transform(
      x,
      crs0 = fm_crs(x),
      crs = fm_crs(crs),
      passthrough = TRUE
    )
    z <- sf::st_coordinates(x)
    x <- z[, intersect(colnames(z), c("X", "Y")), drop = FALSE]
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



#' @describeIn fm_nonconvex_hull
#' Differs from `sf::st_buffer(x, convex)` followed by
#' `sf::st_concave_hull()` (available from GEOS 3.11)
#' in how the amount of allowed concavity is controlled.
#' @param preserveTopology logical; argument to `sf::st_simplify()`
#'   (for `method="sf"` only)
#' @param dTolerance If not zero, controls the `dTolerance` argument to
#'   `sf::st_simplify()`. The default is `pmin(convex, concave) / 40`, chosen to
#'   give approximately 4 or more subsegments per circular quadrant.
#'   (for `method="sf"` only)
fm_nonconvex_hull_sf <- function(x,
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
  ## This workaround produces spurious extra points.
  # # st_buffer can break for LINESTRING input with large dist,
  # # giving interior holes.
  # # Partial protection obtained by taking the union with an extension of the
  # # points.
  # z <-
  #   sf::st_buffer(
  #     sf::st_cast(x, to = "POINT"),
  #     dist = convex + concave,
  #     nQuadSegs = nQuadSegs
  #   )
  # y <- sf::st_union(y, z)
  if (concave > 0) {
    if (sf::sf_use_s2() && isTRUE(sf::st_is_longlat(x))) {
      # s2 gives empty result for negative buffers
      # Use bounding set trick to get around this
      y_box <- sf::st_buffer(y, dist = 100000, nQuadSegs = nQuadSegs)
      y_box <- sf::st_union(y_box)
      y_inverse <- sf::st_difference(y_box, y)
      y_inverse_expanded <-
        sf::st_buffer(y_inverse, concave, nQuadSegs = nQuadSegs)
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





# Methods ####


#' @rdname fm_nonconvex_hull
#' @export
fm_nonconvex_hull.sfc <- function(x,
                                  ...,
                                  format = "sf",
                                  method = "fm") {
  format <- match.arg(format, c("sf", "fm"))
  method <- match.arg(method, c("sf", "fm"))
  if (method == "sf") {
    result <- fm_nonconvex_hull_sf(x, ...)
    if (format == "fm") {
      result <- fm_as_segm(result)
    }
  } else {
    result <- fm_nonconvex_hull_fm(x, ...)
    if (format == "sf") {
      result <- fm_as_sfc(result)
    }
  }
  result
}


#' @rdname fm_nonconvex_hull
#' @export
fm_nonconvex_hull.matrix <- function(x, ..., format = "sf", method = "fm") {
  fm_nonconvex_hull.sfc(sf::st_multipoint(x), ...,
    format = format, method = method
  )
}

#' @rdname fm_nonconvex_hull
#' @export
fm_nonconvex_hull.sf <- function(x, ..., format = "sf", method = "fm") {
  fm_nonconvex_hull.sfc(
    sf::st_geometry(x), ...,
    format = format, method = method
  )
}

#' @rdname fm_nonconvex_hull
#' @export
fm_nonconvex_hull.Spatial <- function(x, ..., format = "sf", method = "fm") {
  fm_nonconvex_hull.sfc(sf::st_as_sfc(x), ..., format = format, method = method)
}

#' @rdname fm_nonconvex_hull
#' @export
fm_nonconvex_hull.sfg <- function(x, ..., format = "sf", method = "fm") {
  fm_nonconvex_hull.sfc(sf::st_sfc(x), ..., format = format, method = method)
}

#' @rdname fm_nonconvex_hull
#' @export
fm_nonconvex_hull.fm_segm <- function(x, ..., format = "sf", method = "fm") {
  fm_nonconvex_hull.sfc(fm_as_sfc(x), ..., format = format, method = method)
}

#' @rdname fm_nonconvex_hull
#' @export
fm_nonconvex_hull.fm_segm_list <- function(x,
                                           ...,
                                           format = "sf",
                                           method = "fm") {
  fm_nonconvex_hull.sfc(fm_as_sfc(x), ..., format = format, method = method)
}

# Legacy methods ####

#' @title Non-convex hull computation
#' @description `r lifecycle::badge("deprecated")`
#'   Legacy method for `INLA::inla.nonconvex.hull()`.
#'   Use [fm_nonconvex_hull()] with `method = "fm"` instead, with
#'   either `format = "fm"` (for compatibility with code
#'   expecting `fm_segm` output) or `format = "sf"`.
#' @seealso [fm_nonconvex_hull()]
#' @param \dots Unused.
#' @inheritParams fm_nonconvex_hull
#' @returns `fm_nonconvex_hull_inla()` returns an [fm_segm]
#' object, for compatibility with `inla.nonconvex.hull()`.
#' @export
#' @keywords internal
#' @family nonconvex inla legacy support
#' @inheritSection fm_mesh_2d INLA compatibility
#' @examplesIf require("splancs")
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
  lifecycle::deprecate_soft(
    "0.4.0.9002",
    "fm_nonconvex_hull_inla()",
    'fm_nonconvex_hull(format = "fm")',
    paste0(
      " The `fm_nonconvex_hull()` method  with `method = \"fm\"` and",
      "\n",
      " `format = \"fm\"` has replaced `fm_nonconvex_hull_inla()`.\n",
      "Most use cases can use `fm_nonconvex_hull(...)` for `sf` output,",
      " which since version `0.4.0.9002` uses the \"fm\" method by default."
    )
  )

  stopifnot(!is.null(x))
  if (inherits(x, c("SpatialPoints", "SpatialPointsDataFrame"))) {
    fm_safe_sp(force = TRUE)
    x <- fm_transform(
      sp::coordinates(x),
      crs0 = fm_crs(x),
      crs = fm_crs(crs),
      passthrough = TRUE
    )
    x <- x[, 1:2, drop = FALSE]
  } else if (inherits(x, c("sf", "sfc"))) {
    x <- fm_transform(
      x,
      crs0 = fm_crs(x),
      crs = fm_crs(crs),
      passthrough = TRUE
    )

    z <- sf::st_coordinates(x)
    z <- z[, intersect(colnames(z), c("X", "Y")), drop = FALSE]

    x <- z
  }

  lim <- rbind(range(x[, 1]), range(x[, 2]))

  approx.diam <- max(diff(lim[1, ]), diff(lim[2, ]))
  if (convex < 0) {
    convex <- -convex * approx.diam
  }
  if (concave < 0) {
    concave <- -concave * approx.diam
  }

  fm_nonconvex_hull(
    x,
    convex = convex,
    concave = concave,
    resolution = 40,
    eps = eps,
    eps_rel = eps_rel,
    crs = crs,
    ...,
    format = "fm",
    method = "fm"
  )
}

#' @export
#' @describeIn fm_nonconvex_hull_inla Special method [fm_nonconvex_hull_fm()]
#'   method for `concave = 0`. Requires `splancs::nndistF()`.
## Based on an idea from Elias Teixeira Krainski
#' @inheritParams fm_nonconvex_hull
#' @keywords internal
fm_nonconvex_hull_inla_basic <- function(x, convex = -0.15, resolution = 40,
                                         eps = NULL, crs = fm_crs(x)) {
  lifecycle::deprecate_soft(
    "0.4.0.9003",
    "fm_nonconvex_hull_inla_basic()",
    I('fm_nonconvex_hull(..., method = "fm", format = "fm", concave = 0)')
  )
  fm_nonconvex_hull_fm_basic(
    x,
    convex = convex,
    resolution = resolution,
    eps = eps,
    crs = crs
  )
}
