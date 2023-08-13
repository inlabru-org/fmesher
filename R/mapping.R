## mapping.R
##
##   Copyright (C) 2015, Finn Lindgren
##
##   This program is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   This program is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
##   You should have received a copy of the GNU General Public License
##   along with this program.  If not, see <http://www.gnu.org/licenses/>.

#' @title Old globe projection methods
#' @description Deprecated globe projection methods that may be removed in the
#' future
#' @keywords internal
#' @param type Projection type, see [.globeproj.types()]
#' @param orient long,lat,rotation
#' @param xlim x-axis limits
#' @param ylim y-axis limits
#' @param scale x- and y- scaling factors
#' @author Finn Lindgren
#' @export
globeproj <- function(type = NULL,
                      orient = NULL,
                      xlim = NULL,
                      ylim = NULL,
                      scale = NULL) {
  if (missing(type)) {
    return(.globeproj.types())
  }
  type <- match.arg(type, .globeproj.types())

  name <- type
  if (identical(type, "lambert")) {
    type <- "orthocyl"
    if (is.null(scale)) scale <- c(1, 1)
    if (is.null(orient)) orient <- c(0, 0, 0, 0)
  } else if (identical(type, "gall-peters")) {
    type <- "orthocyl"
    if (is.null(scale)) scale <- c(1, 1)
    scale <- scale * c(1 * cos(pi / 4), 1 / cos(pi / 4))
    if (is.null(orient)) orient <- c(0, 0, 0, 0)
  }

  if (identical(type, "longlat")) {
    if (is.null(orient)) orient <- c(0, 0, 0, 0)
    if (is.null(scale)) scale <- c(1, 1)
    if (is.null(xlim)) xlim <- c(-180, 180) * scale[1]
    if (is.null(ylim)) ylim <- c(-90, 90) * scale[2]
  } else if (identical(type, "orthocyl")) {
    if (is.null(orient)) orient <- c(0, 0, 0, 0)
    if (is.null(scale)) scale <- c(1, 1)
    if (is.null(xlim)) xlim <- c(-pi, pi) * scale[1]
    if (is.null(ylim)) ylim <- c(-1, 1) * scale[2]
  } else if (identical(type, "mollweide")) {
    if (is.null(orient)) orient <- c(0, 0, 0, 0)
    if (is.null(scale)) scale <- c(1, 1) / sqrt(2)
    if (is.null(xlim)) xlim <- c(-2, 2) * sqrt(2) * scale[1]
    if (is.null(ylim)) ylim <- c(-1, 1) * sqrt(2) * scale[2]
  } else if (identical(type, "hammer")) {
    if (is.null(orient)) orient <- c(0, 0, 0, 0)
    if (is.null(scale)) scale <- c(1, 1) / sqrt(2)
    if (is.null(xlim)) xlim <- c(-2, 2) * sqrt(2) * scale[1]
    if (is.null(ylim)) ylim <- c(-1, 1) * sqrt(2) * scale[2]
  }
  while (length(orient) < 4) {
    orient <- c(orient, 0)
  }

  x <- structure(
    list(
      name = name, type = type,
      orient = orient,
      xlim = xlim, ylim = ylim, scale = scale
    ),
    class = "globeproj"
  )
  .validobject(x)
  x
}


##' @describeIn globeproj Types of globe projections
##' @param x A [globeproj] object.
##' @param \dots Not used.
##' @return A vector of names of available projection types; "longlat",
##' "mollweide", "hammer", "orthocyl", "lambert", "gall-peters"
##' @author Finn Lindgren
.globeproj.types <-
  function(x, ...) {
    c(
      "longlat",
      "mollweide", "hammer",
      "orthocyl",
      "lambert", "gall-peters"
    )
  }

setMethodS3(
  "valid", "globeproj",
  function(object, ...) {
    msg <- c()
    if (!(object$type %in% .globeproj.types())) {
      msg <- c(
        msg,
        paste("'type' should be one of (",
          paste("'", .globeproj.types(), "'",
            collapse = ", ", sep = ""
          ),
          "), not '", object$type, "'.",
          sep = ""
        )
      )
    }
    if (length(object$xlim) != 2) {
      msg <- c(
        msg,
        paste("'xlim' should have length 2, not ",
          length(object$xlim), ".",
          sep = ""
        )
      )
    }
    if (length(object$ylim) != 2) {
      msg <- c(
        msg,
        paste("'ylim' should have length 2, not ",
          length(object$ylim), ".",
          sep = ""
        )
      )
    }
    if (length(msg) == 0) {
      TRUE
    } else {
      msg
    }
  }
)

.validobject <- function(x, test = FALSE, ...) {
  ## Mostly from validObject
  errors <- valid(x)
  if (length(errors) > 1 || is.character(errors[1])) {
    if (test) {
      errors
    } else {
      msg <- gettextf("invalid class %s object", dQuote(class(x)))
      if (length(errors) > 1L) {
        stop(paste(paste0(msg, ":"), paste(seq_along(errors),
          errors,
          sep = ": "
        ),
        collapse = "\n"
        ), domain = NA)
      } else {
        stop(msg, ": ", errors, domain = NA)
      }
    }
  } else {
    TRUE
  }
}





.clip <- function(x, coords) {
  ## Clip and generate a list of Line objects
  thelines <- list()
  ## Rudimentary cutting:
  toolong <-
    which(c(
      TRUE,
      diff(coords[, 1])^2 / diff(x$xlim)
        + diff(coords[, 2])^2 / diff(x$ylim)
      > 0.1^2,
      TRUE
    ))
  start <- toolong[-length(toolong)]
  ending <- toolong[-1] - 1
  for (i in seq_along(start)) {
    coords1 <- coords[start[i]:ending[i], 1]
    coords2 <- coords[start[i]:ending[i], 2]
    thelines <-
      c(
        thelines,
        list(sp::Line(cbind(coords1, coords2)))
      )
  }
  thelines
}


#' @title Plot a projected PolySet
#' @param x A PolySet (see `PBSmapping`)
#' @param projection a `globeproj` objcet
#' @param add logical; If TRUE, add to existing plot
#' @param \dots Additional parameters passed on to `sp::plot`
#' @export
#' @examples
#' # data(worldLL, package = "PBSmapping")
#' # plot_PolySet(worldLL, globeproj("longlat"), add = FALSE)
#' @keywords internal
plot_PolySet <- function(x, projection, add = FALSE, ...) {
  coords <- project(projection, cbind(x$X, x$Y))
  proj <-
    sp::SpatialLines(list(sp::Lines(
      unlist(
        lapply(
          unique(x$PID),
          function(k) {
            .clip(
              projection,
              coords[sum(x$PID < k) +
                x$POS[x$PID == k], , drop = FALSE]
            )
          }
        ),
        recursive = FALSE
      ),
      ID = "shape"
    )))
  if (add) {
    sp::plot(proj, add = TRUE, ...)
  } else {
    sp::plot(proj, ...)
  }
  invisible(proj)
}

# outline ####

#' @param x A [globeproj] object
#' @param add logical; If TRUE, add to existing plot
#' @param do.plot logical; If try, do plotting
#' @param \dots Additional parameters passed to other methods
#' @return A
#' @author Finn Lindgren
#' @export
#' @export outline
#' @aliases outline
#' @rdname globeproj

setMethodS3(
  "outline", "globeproj",
  function(x, add = FALSE, do.plot = TRUE,
           ...) {
    thebox <-
      sp::SpatialPolygons(list(sp::Polygons(
        list(sp::Polygon(
          cbind(
            c(
              x$xlim[1], x$xlim[1],
              x$xlim[2], x$xlim[2], x$xlim[1]
            ),
            c(
              x$ylim[1], x$ylim[2],
              x$ylim[2], x$ylim[1], x$ylim[1]
            )
          )
        )),
        ID = "box"
      )))
    thebox <- sf::st_as_sfc(thebox)
    if (x$type %in% c("longlat", "orthocyl")) {
      theoutline <- thebox
    } else if (x$type %in% c("mollweide", "hammer")) {
      angle <- seq(0, 2 * pi, length = 361)
      circle <- cbind(cos(angle), -sin(angle))
      theboundary <-
        sp::SpatialPolygons(list(sp::Polygons(
          list(sp::Polygon(
            cbind(
              circle[, 1] * 2 * sqrt(2) * x$scale[1],
              circle[, 2] * sqrt(2) * x$scale[2]
            )
          )),
          ID = "boundary"
        )))
      theboundary <- sf::st_as_sfc(theboundary)
      theoutline <- sf::st_intersection(theboundary, thebox)
    }
    theoutline <- sf::as_Spatial(theoutline)
    if (do.plot) {
      if (add) {
        sp::plot(theoutline, add = TRUE, ...)
      } else {
        sp::plot(theoutline, ...)
      }
    }
    invisible(theoutline)
  }
)

# graticule ####

#' @rdname globeproj
#' @param x A [globeproj] object
#' @param n The number of graticules (n-long, n-lat) to compute
#' @param add logical; If TRUE, add to existing plot
#' @param do.plot logical; If TRUE, do plotting
#' @param \dots Additional parameters passed on to other methods
#' @return A
#' @author Finn Lindgren
#' @export
#' @export graticule
#' @aliases graticule
setMethodS3(
  "graticule", "globeproj",
  function(x, n = c(24, 12), add = FALSE, do.plot = TRUE,
           ...) {
    ## Graticule
    if (n[1] > 0) {
      graticule1 <- floor(n[1] / 2)
      lon <- ((-graticule1):graticule1) * 2 / n[1] * 180
      lat <- seq(-90, 90, by = 2)
      meridians <- expand.grid(lat, lon)[, 2:1]
      ## TODO: Add oblique support (requires cutting), etc
      proj.mer.coords <- project(x, meridians)
      proj.mer.coords1 <- matrix(
        proj.mer.coords[, 1], length(lat),
        length(lon)
      )
      proj.mer.coords2 <- matrix(
        proj.mer.coords[, 2], length(lat),
        length(lon)
      )

      proj.mer <-
        sp::SpatialLines(list(sp::Lines(
          unlist(
            lapply(
              seq_along(lon),
              function(k) {
                .clip(x, cbind(
                  proj.mer.coords1[, k],
                  proj.mer.coords2[, k]
                ))
              }
            ),
            recursive = FALSE
          ),
          ID = "meridians"
        )))
      if (do.plot) {
        if (add) {
          sp::plot(proj.mer, add = TRUE, ...)
        } else {
          sp::plot(proj.mer, ...)
          add <- TRUE
        }
      }
    } else {
      proj.mer <- NULL
    }
    if (n[2] > 1) {
      lon <- seq(-180, 180, by = 2)
      lat <- seq(-90, 90, length = 1 + n[2])[-c(1, 1 + n[2])]
      parallels <- expand.grid(lon, lat)
      proj.par.coords <- project(x, parallels)
      proj.par.coords1 <- matrix(
        proj.par.coords[, 1], length(lon),
        length(lat)
      )
      proj.par.coords2 <- matrix(
        proj.par.coords[, 2], length(lon),
        length(lat)
      )
      proj.par <-
        sp::SpatialLines(list(sp::Lines(
          unlist(
            lapply(
              seq_along(lat),
              function(k) {
                .clip(x, cbind(
                  proj.par.coords1[, k],
                  proj.par.coords2[, k]
                ))
              }
            ),
            recursive = FALSE
          ),
          ID = "parallels"
        )))
      if (do.plot) {
        if (add) {
          sp::plot(proj.par, add = TRUE, ...)
        } else {
          sp::plot(proj.par, ...)
        }
      }
    } else {
      proj.par <- NULL
    }
    invisible(list(meridians = proj.mer, parallels = proj.par))
  }
)

# tissot ####

##' @rdname globeproj
##' @param x A [globeproj] object
##' @param n The number of Tissot indicatrices (n-long, n-lat) to compute
##' @param add logical; If TRUE, add to existing plot
##' @param do.plot logical; If TRUE, do plotting
##' @param \dots Additional parameters passed on to other methods
##' @return A
##' @author Finn Lindgren
##' @export
##' @export tissot
##' @aliases tissot
setMethodS3(
  "tissot", "globeproj",
  function(x, n = c(12, 6), add = FALSE, do.plot = TRUE,
           ...) {
    ## Tissot
  }
)




#' @title Plot a globeproj object
#' @param x A [globeproj] object
#' @param xlim,ylim The x- and y-axis limits
#' @param outline logical
#' @param graticule The number of graticules (n-long, n-lat) to compute
#' @param tissot The number of Tissot indicatrices (n-long, n-lat) to compute
#' @param asp the aspect ratio. Default = 1
#' @param add logical; If `TRUE`, add to existing plot. Default: `FALSE`
#' @param \dots Additional parameters passed on to other methods
#' @return Nothing
#' @author Finn Lindgren
#' @export
#' @examples
#' proj <- globeproj("moll", orient = c(0, 0, 45))
#' plot_globeproj(proj, graticule = c(24, 12), add = FALSE, asp = 1, lty = 2, lwd = 0.5)
plot_globeproj <-
  function(x, xlim = NULL, ylim = NULL,
           outline = TRUE,
           graticule = c(24, 12),
           tissot = c(12, 6),
           asp = 1,
           add = FALSE,
           ...) {
    if (is.null(xlim)) xlim <- x$xlim
    if (is.null(ylim)) ylim <- x$ylim
    if (!add) {
      plot.new()
      plot(NA, type = "n", xlim = xlim, ylim = ylim, asp = asp, ...)
    }
    ## Outline
    if (outline) {
      outline(x, add = TRUE, ...)
    }
    ## Graticule
    graticule(x, n = graticule, add = TRUE, ...)
    ## Tissot
    tissot(x, n = tissot, add = TRUE, asp = asp, ...)
    invisible(NULL)
  }


# limits ####

##' @describeIn globeproj Calculates projection axis limits
##'
##' @param x A [globeproj] object
##' @param loc Coordinates to be mapped.
##' @param \dots Additional parameters passed on to other methods
##' @return A list:
##' \item{xlim }{X axis limits in the map domain}
##' \item{ylim }{Y axis limits in the map domain}
##' @author Finn Lindgren
##' @export
##' @export limits
##' @aliases limits
setMethodS3(
  "limits", "globeproj",
  function(x,
           loc = NULL,
           ...) {
    if (!is(x, "globeproj")) {
      stop("'x' must be a 'globeproj'")
    }
    lim <- list(xlim = x$xlim, ylim = x$ylim)
    if (!is.null(loc)) {
      locproj <- project(x, loc)
      xlim <- range(locproj[, 1], na.rm = TRUE)
      ylim <- range(locproj[, 2], na.rm = TRUE)
      lim$xlim <- c(max(lim$xlim[1], xlim[1]), min(lim$xlim[2], xlim[2]))
      lim$ylim <- c(max(lim$ylim[1], ylim[1]), min(lim$ylim[2], ylim[2]))
    }
    lim
  }
)


rotmat3213 <- function(rot) {
  cs <- cos(rot[1])
  sn <- sin(rot[1])
  R <- matrix(c(
    cs, -sn, 0,
    sn, cs, 0,
    0, 0, 1
  ), 3, 3)
  cs <- cos(rot[2])
  sn <- sin(rot[2])
  R <- R %*% matrix(c(
    cs, 0, sn,
    0, 1, 0,
    -sn, 0, cs
  ), 3, 3)
  cs <- cos(rot[3])
  sn <- sin(rot[3])
  R <- R %*% matrix(c(
    1, 0, 0,
    0, cs, -sn,
    0, sn, cs
  ), 3, 3)
  cs <- cos(rot[4])
  sn <- sin(rot[4])
  R <- R %*% matrix(c(
    cs, -sn, 0,
    sn, cs, 0,
    0, 0, 1
  ), 3, 3)
  R
}

rotmat3123 <- function(rot) {
  cs <- cos(rot[4])
  sn <- sin(rot[4])
  R <- matrix(c(
    cs, -sn, 0,
    sn, cs, 0,
    0, 0, 1
  ), 3, 3)
  cs <- cos(rot[3])
  sn <- sin(rot[3])
  R <- R %*% matrix(c(
    1, 0, 0,
    0, cs, -sn,
    0, sn, cs
  ), 3, 3)
  cs <- cos(rot[2])
  sn <- sin(rot[2])
  R <- R %*% matrix(c(
    cs, 0, sn,
    0, 1, 0,
    -sn, 0, cs
  ), 3, 3)
  cs <- cos(rot[1])
  sn <- sin(rot[1])
  R <- R %*% matrix(c(
    cs, -sn, 0,
    sn, cs, 0,
    0, 0, 1
  ), 3, 3)
  R
}

# project ####

#' @rdname globeproj
#' @param x A [globeproj] object
#' @param loc Coordinates to be mapped.
#' @param inverse logical
#' @param \dots Additional parameters passed on to other methods
#' @return B
#' @author Finn Lindgren
#'
#' @export
#' @export project
#' @aliases project
setMethodS3(
  "project", "globeproj",
  function(x, loc, inverse = FALSE, ...) {
    if (!is(x, "globeproj")) {
      stop("'x' must be a 'globeproj'")
    }
    if (inverse) {
      ## Convert coordinates to standardised scale
      loc <- loc %*% diag(1 / x$scale, 2)
    } else {
      if (ncol(loc) == 2) {
        ## Standardise (long,lat) to Cartesian coordinates.
        loc <-
          cbind(
            x = cos(loc[, 1] * pi / 180) * cos(loc[, 2] * pi / 180),
            y = sin(loc[, 1] * pi / 180) * cos(loc[, 2] * pi / 180),
            z = sin(loc[, 2] * pi / 180)
          )
      }
      ## Transform to oblique orientation
      ## 1) Rotate -orient[1] around (0,0,1)
      ## 2) Rotate +orient[2] around (0,1,0)
      ## 3) Rotate -orient[3] around (1,0,0)
      ## 3) Rotate -orient[4] around (0,0,1)
      loc <- loc %*% rotmat3213(c(-1, 1, -1, -1) * x$orient * pi / 180)
    }
    if (identical(x$type, "longlat")) {
      if (inverse) {
        proj <-
          cbind(
            x = cos(loc[, 1] * pi / 180) * cos(loc[, 2] * pi / 180),
            y = sin(loc[, 1] * pi / 180) * cos(loc[, 2] * pi / 180),
            z = sin(loc[, 2] * pi / 180)
          )
      } else {
        proj <-
          cbind(
            x = atan2(loc[, 2], loc[, 1]) * 180 / pi,
            y = asin(pmax(-1, pmin(+1, loc[, 3]))) * 180 / pi
          )
      }
    } else if (identical(x$type, "orthocyl")) {
      if (inverse) {
        coslat <- sqrt(pmax(0, 1 - loc[, 2]^2))
        proj <-
          cbind(
            x = cos(loc[, 1] * pi / 180) * coslat,
            y = sin(loc[, 1] * pi / 180) * coslat,
            z = loc[, 2]
          )
      } else {
        proj <-
          cbind(
            x = atan2(loc[, 2], loc[, 1]),
            y = loc[, 3]
          )
      }
    } else if (identical(x$type, "mollweide")) {
      if (inverse) {
        ok <- ((loc[, 1]^2 / 2 + 2 * loc[, 2]^2) <= 4)
        cos.theta <- sqrt(pmax(0, 1 - loc[ok, 2]^2 / 2))
        theta <- atan2(loc[ok, 2] / sqrt(2), cos.theta)
        sin.lat <- (2 * theta + sin(2 * theta)) / pi
        cos.lat <- sqrt(pmax(0, 1 - sin.lat^2))
        lon <- loc[ok, 1] / sqrt(2) * pi / 2 / (cos.theta + (cos.theta == 0))
        lon[cos.theta == 0] <- pi / 2 * sign(theta[cos.theta == 0])
        proj <- matrix(NA, nrow(loc), 3)
        proj[ok, ] <- cbind(
          x = cos(lon) * cos.lat,
          y = sin(lon) * cos.lat,
          z = sin.lat
        )
      } else {
        lon <- atan2(loc[, 2], loc[, 1])
        z <- pmin(1, pmax(-1, loc[, 3]))
        sin.theta <- z
        cos.theta <- sqrt(pmax(0, 1 - sin.theta^2))
        ## NR-solver for sin.theta.
        ## Typically finishes after at most 7 iterations.
        ## When cos.theta=0, sin.theta is already correct, +/- 1.
        notok <- (cos.theta > 0)
        for (k in 1:20) {
          if (!any(notok)) {
            break
          }
          delta <-
            (atan2(sin.theta[notok], cos.theta[notok]) +
              sin.theta[notok] * cos.theta[notok] - pi / 2 * z[notok]) /
              (2 * cos.theta[notok])
          sin.theta[notok] <- sin.theta[notok] - delta
          cos.theta[notok] <- sqrt(1 - sin.theta[notok]^2)
          notok[notok] <- (abs(delta) > 1e-14)
        }
        proj <- cbind(
          x = 2 * lon / pi * cos.theta * sqrt(2),
          y = sin.theta[, drop = TRUE] * sqrt(2)
        )
      }
    } else if (identical(x$type, "hammer")) {
      if (inverse) {
        ok <- ((loc[, 1]^2 / 2 + 2 * loc[, 2]^2) <= 4)
        z2 <- 1 - (loc[ok, 1] / 4)^2 - (loc[ok, 2] / 2)^2
        z <- sqrt(z2)
        lon <- 2 * atan(z * loc[ok, 1] / (2 * (2 * z2 - 1)))
        sin.lat <- z * loc[ok, 2]
        cos.lat <- sqrt(pmax(0, 1 - sin.lat^2))
        proj <- matrix(NA, nrow(loc), 3)
        proj[ok, ] <- cbind(
          x = cos(lon) * cos.lat,
          y = sin(lon) * cos.lat,
          z = sin.lat
        )
      } else {
        lon <- atan2(loc[, 2], loc[, 1])
        sin.lat <- pmin(1, pmax(-1, loc[, 3]))
        cos.lat <- sqrt(pmax(0, 1 - sin.lat^2))
        cos.lon2 <- cos(lon / 2)
        sin.lon2 <- sin(lon / 2)
        scale <- sqrt(2) / sqrt(1 + cos.lat * cos.lon2)
        proj <- cbind(
          x = 2 * cos.lat * sin.lon2 * scale,
          y = sin.lat * scale
        )
      }
    }
    if (inverse) {
      ## Transform back from oblique orientation
      ## 1) Rotate +orient[4] around (0,0,1)
      ## 2) Rotate +orient[3] around (1,0,0)
      ## 3) Rotate -orient[2] around (0,1,0)
      ## 4) Rotate +orient[1] around (0,0,1)
      proj <- proj %*% rotmat3123(c(1, -1, 1, 1) * x$orient * pi / 180)
    } else {
      proj <- cbind(x = proj[, 1] * x$scale[1], y = proj[, 2] * x$scale[2])
    }
    proj
  }
)
