#' @include deprecated.R

#' @title Recursive curve simplification.
#'
#' @description
#' Helper from legacy `INLA::inla.simplify.curve()`
#'
#' Attempts to simplify a polygonal curve by joining nearly colinear segments.
#'
#' Uses a variation of the binary splitting Ramer-Douglas-Peucker algorithm,
#' with an ellipse of half-width `eps` ellipse instead of a rectangle, motivated
#' by prediction ellipse for Brownian bridge.
#'
#' @param loc Coordinate matrix.
#' @param idx Index vector into `loc` specifying a polygonal curve.
#' @param eps Absolute straightness tolerance. Default `NULL`, no constraint.
#' @param eps_rel Relative straightness tolerance. Default `NULL`, no
#'   constraint.
#' @return An index vector into `loc` specifying the simplified polygonal
#' curve.
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @details
#' Variation of Ramer-Douglas-Peucker.
#' Uses width epsilon ellipse instead of rectangle,
#' motivated by prediction ellipse for Brownian bridge.
#'
#' @examples
#'
#' theta <- seq(0, 2 * pi, length.out = 1000)
#' loc <- cbind(cos(theta), sin(theta))
#' idx <- fm_simplify_helper(loc = loc, idx = 1:nrow(loc), eps = 0.01)
#' print(c(nrow(loc), length(idx)))
#' plot(loc, type = "l")
#' lines(loc[idx, ], col = "red")
#' @export
#' @keywords internal
#' @family nonconvex inla legacy support
fm_simplify_helper <- function(loc, idx, eps = NULL, eps_rel = NULL) {
  n <- length(idx)
  if ((n <= 2) ||
    (is.null(eps) && is.null(eps_rel)) ||
    (min(eps, eps_rel) == 0)) {
    return(idx)
  }
  segm <- loc[idx[n], ] - loc[idx[1], ]
  segm.len <- sum(segm^2)^0.5
  if (segm.len <= 1e-12) {
    ## End point same as start; closed curve.  Split.
    len2 <- ((loc[idx[2:(n - 1)], 1] - loc[idx[1], 1])^2 +
      (loc[idx[2:(n - 1)], 2] - loc[idx[1], 2])^2)
    split <- which.max(len2) + 1L
  } else {
    segm.mid <- (loc[idx[n], ] + loc[idx[1], ]) / 2
    segm <- segm / segm.len
    segm.perp <- c(-segm[2], segm[1])
    vec <- (cbind(
      loc[idx[2:(n - 1)], 1] - segm.mid[1],
      loc[idx[2:(n - 1)], 2] - segm.mid[2]
    ))
    ## Always split if any point is outside the circle
    epsi <- min(c(eps, eps_rel * segm.len / 2, segm.len / 2))
    dist1 <- abs(vec[, 1] * segm[1] + vec[, 2] * segm[2]) /
      (segm.len / 2) * epsi
    dist2 <- abs(vec[, 1] * segm.perp[1] + vec[, 2] * segm.perp[2])
    dist <- (dist1^2 + dist2^2)^0.5

    ## Find the furthest point, in the ellipse metric, and
    ## check if it inside the ellipse (radius=segm.len/2)
    split <- which.max(dist) + 1L
    if (dist[split - 1L] < epsi) {
      ## Flat segment, eliminate.
      return(idx[c(1, n)])
    }
    ## Split at the furthest point.
    split <- which.max(dist) + 1L
  }

  ## Do the split recursively:
  return(c(
    fm_simplify_helper(loc, idx[1L:split], eps = eps, eps_rel = eps_rel),
    fm_simplify_helper(loc, idx[split:n], eps = eps, eps_rel = eps_rel)[-1L]
  ))
}

#' @title Recursive curve simplification.
#'
#' @description `r lifecycle::badge("experimental")`
#' Simplifies polygonal curve segments by joining nearly
#' co-linear segments.
#'
#' Uses a variation of the binary splitting Ramer-Douglas-Peucker algorithm,
#' with an ellipse of half-width `eps` ellipse instead of a rectangle, motivated
#' by prediction ellipse for Brownian bridge.
#'
#' @param x An [fm_segm()] object.
#' @param eps Absolute straightness tolerance. Default `NULL`, no constraint.
#' @param eps_rel Relative straightness tolerance. Default `NULL`, no
#'   constraint.
#' @param ... Currently unused.
#' @return The simplified [fm_segm()] object.
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @details
#' Variation of Ramer-Douglas-Peucker.
#' Uses width epsilon ellipse instead of rectangle,
#' motivated by prediction ellipse for Brownian bridge.
#' @references
#' Ramer, Urs (1972). "An iterative procedure for the polygonal approximation of
#' plane curves". *Computer Graphics and Image Processing*. **1** (3): 244–256.
#' \doi{10.1016/S0146-664X(72)80017-0}
#'
#' Douglas, David; Peucker, Thomas (1973). "Algorithms for the reduction of
#' the number of points required to represent a digitized line or its
#' caricature". *The Canadian Cartographer*. **10** (2): 112–122.
#' \doi{10.3138/FM57-6770-U75U-7727}
#'
#' @examples
#' theta <- seq(0, 2 * pi, length.out = 1000)
#' (segm <- fm_segm(cbind(cos(theta), sin(theta)),
#'   idx = seq_along(theta)
#' ))
#' (segm1 <- fm_simplify(segm, eps_rel = 0.1))
#' (segm2 <- fm_simplify(segm, eps_rel = 0.2))
#' plot(segm)
#' lines(segm1, col = 2)
#' lines(segm2, col = 3)
#'
#' (segm <- fm_segm(cbind(theta, sin(theta * 4)),
#'   idx = seq_along(theta)
#' ))
#' (segm1 <- fm_simplify(segm, eps_rel = 0.1))
#' (segm2 <- fm_simplify(segm, eps_rel = 0.2))
#' plot(segm)
#' lines(segm1, col = 2)
#' lines(segm2, col = 3)
#' @export
#' @family object creation and conversion
fm_simplify <- function(x, eps = NULL, eps_rel = NULL, ...) {
  if ((nrow(x$idx) <= 1) ||
    (is.null(eps) && is.null(eps_rel)) ||
    (min(eps, eps_rel) == 0)) {
    return(x)
  }

  segm_split <- list()
  k <- 0L

  not_handled_seg <- seq_len(nrow(x$idx))
  while (length(not_handled_seg) > 0) {
    next_seg <- not_handled_seg[1]
    seq_seg <- integer(0)
    seq_vtx <- x$idx[next_seg, 1]
    next_vtx <- x$idx[next_seg, 2]
    while (TRUE) {
      final <- next_vtx %in% seq_vtx
      seq_vtx <- c(seq_vtx, next_vtx)
      seq_seg <- c(seq_seg, next_seg)
      not_handled_seg <- setdiff(not_handled_seg, next_seg)
      if (final) {
        break
      }
      if (all(x$is.bnd)) {
        next_seg <-
          not_handled_seg[which(x$idx[not_handled_seg, 1] == next_vtx)]
        if (length(next_seg) == 0) {
          break
        }
        next_seg <- next_seg[1]
        next_vtx <- x$idx[next_seg, 2]
      } else {
        next_seg1 <-
          not_handled_seg[which(x$idx[not_handled_seg, 1] == next_vtx)]
        next_seg2 <-
          not_handled_seg[which(x$idx[not_handled_seg, 2] == next_vtx)]
        if ((length(next_seg1) == 0) && (length(next_seg2) == 0)) {
          break
        }
        if (length(next_seg1) > 0) {
          next_seg <- next_seg1[1]
          next_vtx <- x$idx[next_seg, 2]
        } else {
          next_seg <- next_seg2[1]
          next_vtx <- x$idx[next_seg, 1]
        }
      }
    }

    # seq_vtx and seq_seg
    k <- k + 1
    # TODO: handle geocent data
    idx <- fm_simplify_helper(
      loc = x$loc,
      idx = seq_vtx,
      eps = eps,
      eps_rel = eps_rel
    )
    # TODO: improve granularity of group information.
    segm_split[[k]] <- fm_segm(
      loc = x$loc,
      idx = idx,
      grp = x$grp[seq_seg[1]],
      is.bnd = all(x$is.bnd),
      crs = fm_crs(x)
    )
  }
  segm_split <- fm_segm_join(fm_as_segm_list(segm_split))

  segm_split
}
