#' @include deprecated.R

# fm_diameter ####

#' @title Diameter bound for a geometric object
#'
#' @description
#' Find an upper bound to the convex hull of a point set or function space
#'
#' @param x A point set as an \eqn{n\times d}{n x d} matrix, or an
#' `fm_mesh_2d`/`1d`/`sf` related object.
#' @param manifold Character string specifying the manifold type. Default for
#'   `matrix` input is to treat the point set with Euclidean
#'   \eqn{\mathbb{R}^d}{R^d} metrics.
#'   Use `manifold="S2"` for great circle distances on a sphere centred at the
#'   origin.
#' @param \dots Additional parameters passed on to the submethods.
#' @returns A scalar, upper bound for the diameter of the convex hull of the
#' point set. For multi-domain spaces (e.g. [fm_tensor()] and
#' [fm_collect()]), a vector of upper bounds for each domain is returned.
#' @author Finn Lindgren <Finn.Lindgren@@gmail.com>
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
  diam
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
  z <- sf::st_coordinates(x)
  z <- z[, intersect(colnames(z), c("X", "Y", "Z")), drop = FALSE]
  fm_diameter.matrix(z)
}

#' @rdname fm_diameter
#' @export
fm_diameter.fm_lattice_2d <- function(x, ...) {
  fm_diameter.matrix(x$loc, manifold = fm_manifold(x), ...)
}

#' @rdname fm_diameter
#' @export
fm_diameter.fm_mesh_1d <- function(x, ...) {
  diff(x[["interval"]])
}

#' @rdname fm_diameter
#' @export
fm_diameter.fm_mesh_2d <- function(x, ...) {
  fm_diameter.matrix(x$loc, manifold = fm_manifold(x), ...)
}

#' @rdname fm_diameter
#' @export
fm_diameter.fm_segm <- function(x, ...) {
  fm_diameter.matrix(x$loc, manifold = fm_manifold(x), ...)
}

#' @rdname fm_diameter
#' @export
fm_diameter.fm_mesh_3d <- function(x, ...) {
  fm_diameter.matrix(x$loc, manifold = fm_manifold(x), ...)
}

#' @rdname fm_diameter
#' @export
fm_diameter.fm_tensor <- function(x, ...) {
  vapply(x[["fun_spaces"]], fm_diameter, ..., 1.0)
}

#' @rdname fm_diameter
#' @export
fm_diameter.fm_collect <- function(x, ...) {
  vapply(x[["fun_spaces"]], fm_diameter, ..., 1.0)
}

#' @rdname fm_diameter
#' @export
fm_diameter.fm_list <- function(x, ...) {
  vapply(x, fm_diameter, ..., 1.0)
}
