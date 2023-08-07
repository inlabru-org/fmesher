#' @include deprecated.R

# fm_diameter ####

#' @title Diameter bound for a geometric object
#'
#' @description
#' Find an upper bound to the convex hull of a point set
#'
#' @param x A point set as an \eqn{n\times d}{n x d} matrix, or an
#' `fm_mesh_2d`/`1d`/`sf` related object.
#' @param manifold Character string specifying the manifold type. Default is to
#' treat the point set with Euclidean \eqn{R^d} metrics. Use
#' `manifold="S2"` for great circle distances on the unit sphere (this is
#' set automatically for `fm_fmesh_2d` objects).
#' @param \dots Additional parameters passed on to the submethods.
#' @return A scalar, upper bound for the diameter of the convex hull of the
#' point set.
#' @author Finn Lindgren <finn.lindgren@@gmail.com>
#' @examples
#'
#' fm_diameter(matrix(c(0, 1, 1, 0, 0, 0, 1, 1), 4, 2))
#' @export
fm_diameter <- function(x, ...) {
  if (is.null(x)) {
    return(0.0)
  }
  UseMethod("fm_diameter")
}

#' @rdname fm_diameter
#' @export
fm_diameter.matrix <- function(x, manifold = NULL, ...) {
  if (nrow(x) <= 1) {
    return(0)
  }
  if (ncol(x) == 1) {
    return(diff(range(x)))
  }

  if (identical(manifold, "S2")) {
    radius <- mean(rowSums(x^2)^0.5)
    x <- x / radius
    distance <- function(u, v) {
      2 * asin(pmin(
        1,
        ((u[1] - v[, 1])^2 + (u[2] - v[, 2])^2 + (u[3] - v[, 3])^2)^0.5 / 2
      ))
    }
    center <- colMeans(x)
    tmp <- sqrt(sum(center^2))
    if (tmp < 1e-6 * radius) {
      diam <- pi
    } else {
      center <- center / tmp
      diam <- min(pi, 2 * max(distance(center, x)))
    }
    diam <- diam * radius
    return(diam)
  }

  distance <- function(u, v) {
    d <- 0
    for (k in seq_len(ncol(v))) {
      d <- d + (u[k] - v[, k])^2
    }
    d^0.5
  }
  center <- rep(0, ncol(x))
  for (k in seq_len(ncol(x))) {
    center[k] <- mean(range(x[, k]))
  }
  diam <- 2 * max(distance(center, x))
  return(diam)
}

#' @rdname fm_diameter
#' @export
fm_diameter.sf <- function(x, ...) {
  fm_diameter.sfc(sf::st_geometry(x))
}

#' @rdname fm_diameter
#' @export
fm_diameter.sfg <- function(x, ...) {
  fm_diameter.sfc(sf::st_sfc(x))
}

#' @rdname fm_diameter
#' @export
fm_diameter.sfc <- function(x, ...) {
  fm_diameter.matrix(sf::st_coordinates(x))
}

#' @rdname fm_diameter
#' @export
fm_diameter.fm_lattice_2d <- function(x, ...) {
  fm_diameter.matrix(x$loc, manifold = fm_manifold(x), ...)
}

#' @rdname fm_diameter
#' @export
fm_diameter.fm_segm <- function(x, ...) {
  fm_diameter.matrix(x$loc, manifold = fm_manifold(x), ...)
}

#' @rdname fm_diameter
#' @export
fm_diameter.fm_mesh_2d <- function(x, ...) {
  fm_diameter.matrix(x$loc, manifold = fm_manifold(x), ...)
}

#' @rdname fm_diameter
#' @export
fm_diameter.fm_mesh_1d <- function(x, ...) {
  diff(x[["interval"]])
}

#' @rdname fm_diameter
#' @export
#' @method fm_diameter inla.mesh.1d
fm_diameter.inla.mesh.1d <- function(x, ...) {
  fm_diameter(fm_as_fm(x), ...)
}

#' @rdname fm_diameter
#' @export
#' @method fm_diameter inla.mesh.segment
fm_diameter.inla.mesh.segment <- function(x, ...) {
  fm_diameter(fm_as_fm(x), ...)
}

#' @rdname fm_diameter
#' @export
#' @method fm_diameter inla.mesh.lattice
fm_diameter.inla.mesh.lattice <- function(x, ...) {
  fm_diameter(fm_as_fm(x), ...)
}

#' @rdname fm_diameter
#' @export
#' @method fm_diameter inla.mesh
fm_diameter.inla.mesh <- function(x, ...) {
  fm_diameter(fm_as_fm(x), ...)
}




