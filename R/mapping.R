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

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title A
##' @param x A \link{globeproj} object.
##' @param ... Not used.
##' @return A vector of names of avaliable projection types.
##' @author Finn Lindgren
##' @name types
.globeproj.types <-
  function(x, ...) {
    c("longlat",
      "mollweide", "hammer",
      "orthocyl",
      "lambert", "gall-peters")
  }

setMethodS3(
    "valid", "globeproj",
    function(object, ...) {
      msg <- c()
      if (!(object$type %in% .globeproj.types())) {
        msg <- c(msg,
                 paste("'type' should be one of (",
                       paste("'", .globeproj.types(), "'",
                             collapse=", ", sep=""),
                       "), not '", object$type, "'.", sep=""))
      }
      if (length(object$xlim) != 2) {
        msg <- c(msg,
                 paste("'xlim' should have length 2, not ",
                       length(object$xlim), ".", sep=""))
      }
      if (length(object$ylim) != 2) {
        msg <- c(msg,
                 paste("'ylim' should have length 2, not ",
                       length(object$ylim), ".", sep=""))
      }
      if (length(msg) == 0) {
        TRUE
      } else {
        msg
      }
    })

.validobject <- function(x, test=FALSE, ...) {
  ## Mostly from validObject
  errors <- valid(x)
  if (length(errors) > 1 || is.character(errors[1])) {
    if (test) {
      errors
    } else {
      msg <- gettextf("invalid class %s object", dQuote(class(x)))
      if (length(errors) > 1L)
        stop(paste(paste0(msg, ":"), paste(seq_along(errors),
                                           errors, sep = ": "),
                   collapse = "\n"), domain = NA)
      else stop(msg, ": ", errors, domain = NA)
    }
  } else {
    TRUE
  }
}


##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title A
##' @param type
##' @param datum long,lat,rotation
##' @param xlim
##' @param ylim
##' @param scale
##' @return B
##' @author Finn Lindgren
##' @export
globeproj <- function(type=NULL,
                      datum=NULL,
                      xlim=NULL,
                      ylim=NULL,
                      scale=NULL)
{
  if (missing(type)) {
    return(.globeproj.types())
  }
  type <- match.arg(type, .globeproj.types())

  name <- type
  if (identical(type, "lambert")) {
    type <- "orthocyl"
    if (is.null(scale)) scale <- c(1, 1)
    if (is.null(datum)) datum <- c(0, 0, 0)
  } else if (identical(type, "gall-peters")) {
    type <- "orthocyl"
    if (is.null(scale)) scale <- c(1, 1)
    scale <- scale * c(1*cos(pi/4), 1/cos(pi/4))
    if (is.null(datum)) datum <- c(0, 0, 0)
  }

  if (identical(type, "longlat")) {
    if (is.null(datum)) datum <- c(0, 0, 0)
    if (is.null(scale)) scale <- c(1,1)
    if (is.null(xlim)) xlim <- c(-180,180)*scale[1]
    if (is.null(ylim)) ylim <- c(-90,90)*scale[2]
  } else if (identical(type, "orthocyl")) {
    if (is.null(datum)) datum <- c(0, 0, 0)
    if (is.null(scale)) scale <- c(1, 1)
    if (is.null(xlim)) xlim <- c(-pi,pi)*scale[1]
    if (is.null(ylim)) ylim <- c(-1,1)*scale[2]
  } else if (identical(type, "mollweide")) {
    if (is.null(datum)) datum <- c(0, 0, 0)
    if (is.null(scale)) scale <- c(1,1)/sqrt(2)
    if (is.null(xlim)) xlim <- c(-2,2)*sqrt(2)*scale[1]
    if (is.null(ylim)) ylim <- c(-1,1)*sqrt(2)*scale[2]
  } else if (identical(type, "hammer")) {
    if (is.null(datum)) datum <- c(0, 0, 0)
    if (is.null(scale)) scale <- c(1,1)/sqrt(2)
    if (is.null(xlim)) xlim <- c(-2,2)*sqrt(2)*scale[1]
    if (is.null(ylim)) ylim <- c(-1,1)*sqrt(2)*scale[2]
  }

  x <- structure(list(name=name, type=type,
                      datum=datum,
                      xlim=xlim, ylim=ylim, scale=scale),
                 class="globeproj")
  .validobject(x)
  x
}



##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title A
##' @param x
##' @param n
##' @param add
##' @param do.plot
##' @param ...
##' @return A
##' @author Finn Lindgren
##' @method graticule globeproj
##' @export graticule
##' @export graticule.globeproj
##' @rdname graticule
##' @name graticule.globeproj
setMethodS3(
    "graticule", "globeproj",
graticule <-     function(x, n=c(24, 12), add=FALSE, do.plot=TRUE,
             ...)
{
      ## Graticule
      if (n[1] > 0) {
        graticule1 <- floor(n[1]/2)
        lon <- ((-graticule1):graticule1) * 2/n[1]*180
        lat <- seq( -90,  90, by=2)
        meridians <- expand.grid(lat, lon)[,2:1]
        ## TODO: Add oblique support, etc
        proj.mer.coords <- project(x, meridians)
        proj.mer.coords1 <- matrix(proj.mer.coords[,1], length(lat),
                                   length(lon))
        proj.mer.coords2 <- matrix(proj.mer.coords[,2], length(lat),
                                   length(lon))
        proj.mer <-
          sp::SpatialLines(
              list(
                  sp::Lines(
                      lapply(seq_along(lon),
                             function(k) sp::Line(cbind(proj.mer.coords1[,k],
                                                        proj.mer.coords2[,k]))),
                      ID="meridians")))
        if (do.plot) {
          if (add) {
            lines(proj.mer, ...)
          } else {
            plot(proj.mer, ...)
            add <- TRUE
          }
        }
      } else {
        proj.mer <- NULL
      }
      if (n[2] > 1) {
        lon <- seq(-180, 180, by=2)
        lat <- seq( -90,  90, length=1+n[2])[-c(1,1+n[2])]
        parallels <- expand.grid(lon, lat)
        proj.par.coords <- project(x, parallels)
        proj.par.coords1 <- matrix(proj.par.coords[,1], length(lon),
                                   length(lat))
        proj.par.coords2 <- matrix(proj.par.coords[,2], length(lon),
                                   length(lat))
        proj.par <-
          sp::SpatialLines(
              list(
                  sp::Lines(
                      lapply(seq_along(lat),
                             function(k) sp::Line(cbind(proj.par.coords1[,k],
                                                        proj.par.coords2[,k]))),
                      ID="meridians")))
        if (do.plot) {
          if (add) {
            lines(proj.par, ...)
          } else {
            plot(proj.par, ...)
          }
        }
      } else {
        proj.par <- NULL
      }
      invisible(list(meridians=proj.mer, parallels=proj.par))
    })


##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title A
##' @param x
##' @param n
##' @param add
##' @param do.plot
##' @param ...
##' @return A
##' @author Finn Lindgren
##' @method tissot globeproj
##' @export tissot
##' @export tissot.globeproj
##' @rdname tissot
##' @name tissot.globeproj
setMethodS3(
    "tissot", "globeproj",
    function(x, n=c(12,6), add=FALSE, do.plot=TRUE,
             ...)
    {
      ## Tissot

    })




##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title The title
##' @param x The object
##' @param outline
##' @param graticule
##' @param tissot
##' @param ...
##' @return Nothing
##' @author Finn Lindgren
##' @usage \method{plot}{globeproj}(x, ...)
##' @rdname plot.globeproj
##' @name plot.globeproj
##' @method plot globeproj
##' @export plot.globeproj
setMethodS3(
    "plot", "globeproj",
    function(x, xlim=NULL, ylim=NULL,
             outline=TRUE,
             graticule=c(24, 12),
             tissot=c(12,6),
             asp=1,
             ...)
    {
      message("Plot something!")
      if (is.null(xlim)) xlim <- x$xlim
      if (is.null(ylim)) ylim <- x$ylim
      plot(NA, type="n", xlim=xlim, ylim=ylim, asp=asp, ...)
      ## Outline
      if (outline) {
        message("'outline' option not implemented")
      }
      ## Graticule
      graticule(x, n=graticule, add=TRUE, ...)
      ## Tissot
      tissot(x, n=tissot, add=TRUE, asp=asp, ...)
    })




##' Calculates projection axis limits
##'
##' @title Projection axis limits
##' @param x
##' @param loc Coordinates to be mapped.
##' @param ...
##' @return A list:
##' \item{xlim }{X axis limits in the map domain}
##' \item{ylim }{Y axis limits in the map domain}
##' @author Finn Lindgren
##' @usage \method{limits}{globeproj}(x, ...)
##' @rdname limits
##' @name limits.globeproj
##' @method limits globeproj
##' @export limits
##' @export limits.globeproj
setMethodS3(
    "limits", "globeproj",
    function(x,
             loc=NULL,
             ...)
    {
      if (!is(x, "globeproj")) {
        stop("'x' must be a 'globeproj'")
      }
      lim <- list(xlim=x$xlim, ylim=x$ylim)
      if (!is.null(loc)) {
        locproj <- project(x, loc)
        xlim <- range(locproj[,1], na.rm=TRUE)
        ylim <- range(locproj[,2], na.rm=TRUE)
        lim$xlim <- c(max(lim$xlim[1], xlim[1]), min(lim$xlim[2], xlim[2]))
        lim$ylim <- c(max(lim$ylim[1], ylim[1]), min(lim$ylim[2], ylim[2]))
      }
      lim
    })


rotmat321 <- function(rot)
{
  cs <- cos(rot[1])
  sn <- sin(rot[1])
  R <- matrix(c(cs, -sn, 0,
                sn, cs, 0,
                0, 0, 1), 3, 3)
  cs <- cos(rot[2])
  sn <- sin(rot[2])
  R <- R %*% matrix(c(cs, 0, -sn,
                      0, 1, 0,
                      sn, 0, cs), 3, 3)
  cs <- cos(rot[3])
  sn <- sin(rot[3])
  R <- R %*% matrix(c(1, 0, 0,
                      0, cs, -sn,
                      0, sn, cs,
                      0, 0, 1), 3, 3)
  R
}

rotmat123 <- function(rot)
{
  cs <- cos(rot[3])
  sn <- sin(rot[3])
  R <- matrix(c(1, 0, 0,
                0, cs, -sn,
                0, sn, cs,
                0, 0, 1), 3, 3)
  cs <- cos(rot[2])
  sn <- sin(rot[2])
  R <- R %*% matrix(c(cs, 0, -sn,
                      0, 1, 0,
                      sn, 0, cs), 3, 3)
  cs <- cos(rot[1])
  sn <- sin(rot[1])
  R <- R %*% matrix(c(cs, -sn, 0,
                      sn, cs, 0,
                      0, 0, 1), 3, 3)
  R
}


##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title A
##' @param x
##' @param loc
##' @param inverse
##' @return B
##' @author Finn Lindgren
##'
##' @usage \method{project}{globeproj}(x, ...)
##' @rdname project
##' @name project.globeproj
##' @method project globeproj
##' @export project
##' @export project.globeproj
setMethodS3(
    "project", "globeproj",
    function(x, loc, inverse=FALSE, ...)
    {
      if (!is(x, "globeproj")) {
        stop("'x' must be a 'globeproj'")
      }
      if (inverse) {
        ## Convert coordinates to standardised scale
        loc <- loc %*% diag(1/x$scale, 2)
      } else {
        if(ncol(loc) == 2) {
          ## Standardise (long,lat) to Cartesian coordinates.
          loc <-
            cbind(x=cos(loc[,1]*pi/180)*cos(loc[,2]*pi/180),
                  y=sin(loc[,1]*pi/180)*cos(loc[,2]*pi/180),
                  z=sin(loc[,2]*pi/180))
        }
        ## Transform to oblique orientation
        ## 1) Rotate -datum[1] around (0,0,1)
        ## 2) Rotate -datum[2] around (0,1,0)
        ## 3) Rotate -datum[3] around (1,0,0)
        loc <- loc %*% rotmat321(-x$datum*pi/180)

      }
      if (identical(x$type, "longlat")) {
        if (inverse) {
          proj <-
            cbind(x=cos(loc[,1]*pi/180)*cos(loc[,2]*pi/180),
                  y=sin(loc[,1]*pi/180)*cos(loc[,2]*pi/180),
                  z=sin(loc[,2]*pi/180))
        } else {
          proj <-
            cbind(x=atan2(loc[,2], loc[,1])*180/pi,
                  y=asin(pmax(-1, pmin(+1, loc[,3])))*180/pi)
        }
      } else if (identical(x$type, "orthocyl")) {
        if (inverse) {
          coslat <- sqrt(pmax(0, 1-loc[,2]^2))
          proj <-
            cbind(x=cos(loc[,1]*pi/180)*coslat,
                  y=sin(loc[,1]*pi/180)*coslat,
                  z=loc[,2]
                  )
        } else {
          proj <-
            cbind(x=atan2(loc[,2], loc[,1]),
                  y=loc[,3])
        }
      } else if (identical(x$type, "mollweide")) {
        if (inverse) {
          ok <- ((loc[,1]^2/2 + 2*loc[,2]^2) <= 4)
          cos.theta <- sqrt(pmax(0, 1-loc[ok,2]^2))
          theta <- atan2(loc[ok,2], cos.theta)
          sin.lat <- (2*theta+sin(2*theta))/pi
          cos.lat = sqrt(pmax(0, 1-sin.lat^2))
          lon <- loc[ok,1]*pi/2/(cos.theta+(cos.theta==0))
          lon[cos.theta==0] <- pi/2*sign(theta[cos.theta==0])
          proj <- matrix(NA, nrow(loc), 3)
          proj[ok,] <- cbind(x=cos(lon)*cos.lat,
                             y=sin(lon)*cos.lat,
                             z=sin.lat
                             )
        } else {
          lon <- atan2(loc[,2], loc[,1])
          z <- pmin(1, pmax(-1, loc[,3]))
          sin.theta <- z
          cos.theta <- sqrt(pmax(0, 1-sin.theta^2))
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
                 sin.theta[notok]*cos.theta[notok] - pi/2*z[notok])/
                (2*cos.theta[notok])
            sin.theta[notok] = sin.theta[notok] - delta
            cos.theta[notok] = sqrt(1-sin.theta[notok]^2)
            notok[notok] = (abs(delta) > 1e-14)
          }
          proj <- cbind(x=2*lon/pi*cos.theta * sqrt(2),
                        y=sin.theta[,drop=TRUE] * sqrt(2))
        }
      } else if (identical(x$type, "hammer")) {
        if (inverse) {
          ok <- ((loc[,1]^2/2 + 2*loc[,2]^2) <= 4)
          z2 <- 1-(loc[ok,1]/4)^2 - (loc[ok,2]/2)^2
          z <- sqrt(z2)
          lon <- 2*atan(z*loc[ok,1]/(2*(2*z2-1)))
          sin.lat <- z*loc[ok,2]
          cos.lat <- sqrt(pmax(0, 1-sin.lat^2))
          proj <- matrix(NA, nrow(loc), 3)
          proj[ok,] <- cbind(x=cos(lon)*cos.lat,
                             y=sin(lon)*cos.lat,
                             z=sin.lat
                             )
        } else {
          lon <- atan2(loc[,2], loc[,1])
          sin.lat <- pmin(1, pmax(-1, loc[,3]))
          cos.lat <- sqrt(pmax(0, 1-sin.lat^2))
          cos.lon2 <- cos(lon/2)
          sin.lon2 <- sin(lon/2)
          scale <- sqrt(2)/sqrt(1 + cos.lat*cos.lon2)
          proj <- cbind(x=2*cos.lat*sin.lon2 * scale,
                        y=sin.lat * scale)
        }
      }
      if (inverse) {
        ## Transform back from oblique orientation
        ## 1) Rotate datum[3] around (1,0,0)
        ## 2) Rotate datum[2] around (0,1,0)
        ## 3) Rotate datum[1] around (0,0,1)
        loc <- loc %*% rotmat123(x$datum*pi/180)
      } else {
        proj <- cbind(x=proj[,1] * x$scale[1], y=proj[,2] * x$scale[2])
      }
      proj
    })
