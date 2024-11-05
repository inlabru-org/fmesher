## mapping.R
##
##   Copyright (C) 2015, Finn Lindgren




fm_internal_clip <- function(bounds, coords, eps = 0.05) {
  stopifnot(requireNamespace("sp"))
  ## Clip 2D coordinate matrix of polylines and generate a list of Line objects
  ## bounds is from inla.crs.bounds
  ## This implementation only removes "long" line segments.
  thelines <- list()
  ## Rudimentary cutting:
  toolong <-
    which(c(
      TRUE,
      (diff(coords[, 1]) / diff(bounds$xlim))^2
        + (diff(coords[, 2]) / diff(bounds$ylim))^2
      > eps^2,
      TRUE
    ))
  start <- toolong[-length(toolong)]
  ending <- toolong[-1] - 1
  for (i in seq_along(start)) {
    if (start[i] < ending[i]) {
      thelines <-
        c(
          thelines,
          list(sp::Line(coords[start[i]:ending[i], 1:2]))
        )
    }
  }
  thelines
}



#' @title Plot CRS and fm_crs objects
#'
#' @description `r lifecycle::badge("experimental")` Plot the outline of a `crs`
#'   or [fm_crs()] projection, with optional graticules (transformed parallels
#'   and meridians) and Tissot indicatrices.
#'
#' @param x A `crs` or [fm_crs()] object.
#' @param xlim Optional x-axis limits.
#' @param ylim Optional y-axis limits.
#' @param outline Logical, if `TRUE`, draw the outline of the projection.
#' @param graticule Vector of length at most 3, to plot meridians with spacing
#' `graticule[1]` degrees and parallels with spacing `graticule[2]`
#' degrees. `graticule[3]` optionally specifies the spacing above and
#' below the first and last parallel.  When `graticule[1]==0` no meridians
#' are drawn, and when `graticule[2]==0` no parallels are drawn. Use
#' `graticule=NULL` to skip drawing a graticule.
#' @param tissot Vector of length at most 3, to plot Tissot's indicatrices with
#' spacing `tissot[1]` degrees and parallels with spacing `tissot[2]`
#' degrees. `tissot[3]` specifices a scaling factor.  Use
#' `tissot=NULL` to skip drawing a Tissot's indicatrices.
#' @param asp The aspect ratio for the plot, default 1.
#' @param add If `TRUE`, add the projecton plot to an existing plot.
#' @param eps Clipping tolerance for rudimentary boundary clipping
#' @param \dots Additional arguments passed on to the internal calls to
#' `plot` and `lines`.
#' @returns `NULL`, invisibly
#' @author Finn Lindgren <finn.lindgren@@gmail.com>
#' @seealso [fm_crs()]
#' @examples
#' \donttest{
#' if (require("sf") && require("sp")) {
#'   for (projtype in c(
#'     "longlat_norm",
#'     "lambert_norm",
#'     "mollweide_norm",
#'     "hammer_norm"
#'   )) {
#'     fm_crs_plot(fm_crs(projtype), main = projtype)
#'   }
#' }
#'
#' if (require("sf") && require("sp")) {
#'   oblique <- c(0, 45, 45, 0)
#'   for (projtype in c(
#'     "longlat_norm",
#'     "lambert_norm",
#'     "mollweide_norm",
#'     "hammer_norm"
#'   )) {
#'     fm_crs_plot(
#'       fm_crs(projtype, oblique = oblique),
#'       main = paste("oblique", projtype)
#'     )
#'   }
#' }
#' }
#' @export
fm_crs_plot <- function(x, xlim = NULL, ylim = NULL,
                        outline = TRUE,
                        graticule = c(15, 15, 45),
                        tissot = c(30, 30, 30),
                        asp = 1,
                        add = FALSE,
                        eps = 0.05,
                        ...) {
  stopifnot(requireNamespace("sp"))

  bounds <- fm_crs_bounds(x)
  if (is.null(xlim)) xlim <- bounds$xlim
  if (is.null(ylim)) ylim <- bounds$ylim
  if (!add) {
    args <- list(x = NA, type = "n", xlim = xlim, ylim = ylim, asp = asp, ...)
    args <- args[intersect(names(args), names(formals(graphics::plot.default)))]
    do.call(sp::plot, args)
  }
  ## Outline
  if (outline) {
    args <- list(x = bounds$polygon, ...)
    args <- args[intersect(names(args), union(
      names(formals(graphics::lines.default)),
      names(formals(graphics::plot.xy))
    ))]
    do.call(lines, args)
  }
  ## Graticule
  fm_crs_graticule(
    x,
    by = graticule,
    add = TRUE,
    do.plot = TRUE,
    eps = eps,
    ...
  )
  ## Tissot
  fm_crs_tissot(
    x,
    by = tissot,
    add = TRUE,
    do.plot = TRUE,
    eps = eps,
    ...
  )
  invisible(NULL)
}


#' @describeIn fm_crs_plot
#'   `r lifecycle::badge("experimental")` Constructs graticule
#'   information for a given `CRS` or [fm_crs()] and optionally plots the
#'   graticules.
#' Returns a list with two elements, `meridians` and `parallels`, which are
#'  `SpatialLines` objects.
#' @param by The spacing between `(long, lat, long_at_poles)`
#'   graticules/indicatrices, see the `graticule` and `tissot` arguments.
#' @param do.plot logical; If TRUE, do plotting
#' @export
fm_crs_graticule <- function(x, by = c(15, 15, 45), add = FALSE, do.plot = TRUE,
                             eps = 0.05, ...) {
  stopifnot(requireNamespace("sp"))

  ## Graticule
  if (is.null(by)) {
    return(invisible(list()))
  }
  if (length(by) < 2) {
    by <- by[c(1, 1, 1)]
  } else if (length(by) < 3) {
    by <- by[c(1, 2, 1)]
  }
  n <- c(floor(180 / by[1]), ceiling(90 / by[2]) - 1, floor(180 / by[3]))
  bounds <- fm_crs_bounds(x)
  if (by[1] > 0) {
    special.poles <- (by[1] != by[3]) && (by[2] > 0)
    lon <- ((1 - n[1]):n[1]) * by[1]
    if (special.poles) {
      lat <- seq(-n[2] * by[2], n[2] * by[2], length.out = 91)
    } else {
      lat <- seq(-90 + 1e-6, 90 - 1e-6, length.out = 91)
    }
    meridians <- as.matrix(expand.grid(lat, lon)[, 2:1])
    proj.mer.coords <- fm_transform(meridians,
      crs0 = fmesher::fm_CRS("longlat_norm"),
      crs = x
    )
    proj.mer.coords1 <- matrix(
      proj.mer.coords[, 1], length(lat),
      length(lon)
    )
    proj.mer.coords2 <- matrix(
      proj.mer.coords[, 2], length(lat),
      length(lon)
    )

    mer.coords <-
      unlist(
        lapply(
          seq_along(lon),
          function(k) {
            fm_internal_clip(bounds, cbind(
              proj.mer.coords1[, k, drop = FALSE],
              proj.mer.coords2[, k, drop = FALSE]
            ),
            eps = eps
            )
          }
        ),
        recursive = FALSE
      )

    if (special.poles) {
      if (by[3] > 0) {
        lon <- ((1 - n[3]):n[3]) * by[3]
        lat <- seq(-90 + 1e-6, -n[2] * by[2],
          length.out = ceiling((90 - n[2] * by[2]) / 2) + 1
        )
        meridians <- as.matrix(expand.grid(lat, lon)[, 2:1])
        proj.mer.coords <- fm_transform(meridians,
          crs0 = fmesher::fm_CRS("longlat_norm"),
          crs = x
        )
        proj.mer.coords1 <- matrix(
          proj.mer.coords[, 1], length(lat),
          length(lon)
        )
        proj.mer.coords2 <- matrix(
          proj.mer.coords[, 2], length(lat),
          length(lon)
        )
        mer.coords <-
          c(
            mer.coords,
            unlist(
              lapply(
                seq_along(lon),
                function(k) {
                  fm_internal_clip(bounds, cbind(
                    proj.mer.coords1[, k, drop = FALSE],
                    proj.mer.coords2[, k, drop = FALSE]
                  ),
                  eps = eps
                  )
                }
              ),
              recursive = FALSE
            )
          )

        lat <- seq(n[2] * by[2],
          90 - 1e-6,
          length.out = ceiling((90 - n[2] * by[2]) / 2) + 1
        )
        meridians <- as.matrix(expand.grid(lat, lon)[, 2:1])
        proj.mer.coords <- fm_transform(meridians,
          crs0 = fmesher::fm_CRS("longlat_norm"),
          crs = x
        )
        proj.mer.coords1 <- matrix(
          proj.mer.coords[, 1], length(lat),
          length(lon)
        )
        proj.mer.coords2 <- matrix(
          proj.mer.coords[, 2], length(lat),
          length(lon)
        )
        mer.coords <-
          c(
            mer.coords,
            unlist(
              lapply(
                seq_along(lon),
                function(k) {
                  fm_internal_clip(bounds, cbind(
                    proj.mer.coords1[, k, drop = FALSE],
                    proj.mer.coords2[, k, drop = FALSE]
                  ),
                  eps = eps
                  )
                }
              ),
              recursive = FALSE
            )
          )
      }
    }

    proj.mer <-
      sp::SpatialLines(
        list(sp::Lines(mer.coords, ID = "meridians")),
        proj4string = fm_CRS(x, oblique = NA)
      )
    if (do.plot) {
      args <- list(x = proj.mer, ...)
      if (add) {
        args[["add"]] <- TRUE
      }
      do.call(sp::plot, args)
      add <- TRUE
    }
  } else {
    proj.mer <- NULL
  }
  if (by[2] > 0) {
    lon <- seq(-180 + 1e-6, 180 - 1e-6, length.out = 181)
    lat <- ((-n[2]):n[2]) * by[2]
    parallels <- as.matrix(expand.grid(lon, lat))
    proj.par.coords <- fm_transform(parallels,
      crs0 = fmesher::fm_CRS("longlat_norm"),
      crs = x
    )
    proj.par.coords1 <- matrix(
      proj.par.coords[, 1], length(lon),
      length(lat)
    )
    proj.par.coords2 <- matrix(
      proj.par.coords[, 2], length(lon),
      length(lat)
    )
    proj.par <-
      sp::SpatialLines(
        list(sp::Lines(
          unlist(
            lapply(
              seq_along(lat),
              function(k) {
                fm_internal_clip(bounds, cbind(
                  proj.par.coords1[, k],
                  proj.par.coords2[, k]
                ), eps = eps)
              }
            ),
            recursive = FALSE
          ),
          ID = "parallels"
        )),
        proj4string = fm_CRS(x, oblique = NA)
      )
    if (do.plot) {
      args <- list(x = proj.par, ...)
      if (add) {
        args[["add"]] <- TRUE
      }
      do.call(sp::plot, args)
    }
  } else {
    proj.par <- NULL
  }
  invisible(list(meridians = proj.mer, parallels = proj.par))
}

#' @describeIn fm_crs_plot
#' `r lifecycle::badge("experimental")` Constructs Tissot indicatrix information
#' for a given `CRS` or [fm_crs()] and optionally plots the indicatrices.
#' Returns a list with one element, `tissot`, which is a `SpatialLines` object.
#' @param diff.eps Pre-scaling
#' @export
fm_crs_tissot <- function(x, by = c(30, 30, 30), add = FALSE, do.plot = TRUE,
                          eps = 0.05, diff.eps = 1e-2, ...) {
  stopifnot(requireNamespace("sp"))

  if (is.null(by)) {
    return(invisible(list()))
  }
  if (length(by) < 2) {
    by <- c(by[1], by[1], 30)
  } else if (length(by) < 3) {
    by <- c(by[1:2], 30)
  }
  bounds <- fm_crs_bounds(x)
  n <- c(floor(180 / by[1]), ceiling(90 / by[2]) - 1)
  lon <- ((1 - n[1]):n[1]) * by[1]
  lat <- ((-n[2]):n[2]) * by[2]
  loc0.lon <- loc0.lat <- loc0 <-
    cbind(as.matrix(expand.grid(lat, lon)[, 2:1]), 0)
  loc0.lon[, 1] <- loc0.lon[, 1] + diff.eps / cos(loc0.lat[, 2] * pi / 180)
  loc0.lat[, 2] <- loc0.lat[, 2] + diff.eps
  crs.longlat <- fm_CRS("longlat_norm")

  loc1 <- fm_transform(loc0, crs0 = crs.longlat, crs = x)
  loc1.lon <- fm_transform(loc0.lon, crs0 = crs.longlat, crs = x)
  loc1.lat <- fm_transform(loc0.lat, crs0 = crs.longlat, crs = x)
  ok <- (rowSums(is.na(loc1)) +
    rowSums(is.na(loc1.lon)) +
    rowSums(is.na(loc1.lat)) == 0)
  loc1 <- loc1[ok, , drop = FALSE]
  loc1.lon <- loc1.lon[ok, , drop = FALSE]
  loc1.lat <- loc1.lat[ok, , drop = FALSE]

  diff.lon <- (loc1.lon - loc1) / eps
  diff.lat <- (loc1.lat - loc1) / eps

  scale <- by[3]
  theta <- seq(0, 2 * pi, length.out = 181)
  ct <- cos(theta) * scale
  st <- sin(theta) * scale

  collection <-
    sp::SpatialLines(
      list(sp::Lines(
        unlist(
          lapply(
            seq_len(nrow(loc1)),
            function(k) {
              loc1.ellipse <-
                cbind(
                  loc1[k, 1] + diff.lon[k, 1] * ct + diff.lat[k, 1] * st,
                  loc1[k, 2] + diff.lon[k, 2] * ct + diff.lat[k, 2] * st
                )
              fm_internal_clip(bounds, loc1.ellipse, eps = eps)
            }
          ),
          recursive = FALSE
        ),
        ID = "parallels"
      )),
      proj4string = fm_CRS(x, oblique = NA)
    )
  if (do.plot) {
    args <- list(x = collection, ...)
    if (add) {
      args[["add"]] <- TRUE
    }
    do.call(sp::plot, args)
  }

  invisible(list(tissot = collection))
}
