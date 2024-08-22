# Point/mesh connection methods ####

# fm_evaluate ####

#' @title Methods for projecting to/from mesh objects
#'
#' @description Calculate evaluation information and/or evaluate a function
#' defined on a mesh or function space.
#'
#' @param mesh An `inla.mesh` or `inla.mesh.1d` object.
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
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
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

    if (is.null(dim(field))) {
      if (is.null(projector$lattice)) {
        data <- as.vector(projector$proj$A %*% as.vector(field))
        data[!projector$proj$ok] <- NA
        return(data)
      } else {
        data <- as.vector(projector$proj$A %*% as.vector(field))
        data[!projector$proj$ok] <- NA
        return(matrix(
          data,
          projector$lattice$dims[1],
          projector$lattice$dims[2]
        ))
      }
    } else if (inherits(field, "sparseMatrix")) {
      data <- projector$proj$A %*% field
      data[!projector$proj$ok, ] <- NA
      return(data)
    } else {
      data <- as.matrix(projector$proj$A %*% field)
      data[!projector$proj$ok, ] <- NA
      return(data)
    }
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
      return(data)
    } else if (inherits(field, "sparseMatrix")) {
      data <- basis$A %*% field
      data[!basis$ok, ] <- NA
      return(data)
    } else {
      data <- as.matrix(basis$A %*% field)
      data[!basis$ok, ] <- NA
      return(data)
    }
  }


# fm_evaluator ####

#' @describeIn fm_evaluate
#' Returns an `fm_evaluator` list object with evaluation information.
#' The `proj` element is a `fm_basis` object, containing (at least)
#' a mapping matrix `A` and a logical vector `ok`, that indicates which
#' locations were mappable to the input mesh.
#' For `fm_mesh_2d` and `inla.mesh`
#' input, `proj` also contains a matrix `bary` and vector `t`, with the
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

#' @title Internal helper functions for mesh field evaluation
#'
#' @description Methods called internally by [fm_evaluator()] methods.
#' @param weights Optional weight vector, one weight for each location
#' @param derivatives logical; If true, also return matrices `dA` and `d2A`
#' for `fm_mesh_1d` objects, and `dx`, `dy`, `dz` for `fm_mesh_2d`.
#' @inheritParams fm_evaluate
#' @export
#' @keywords internal
#' @returns A list of evaluator information objects, at least a matrix `A` and
#' logical vector `ok`.
#' @name fm_evaluator_helpers
#' @examples
#' str(fm_evaluator_mesh_2d(fmexample$mesh, loc = fmexample$loc))
#'
fm_evaluator_mesh_2d <- function(mesh,
                                 loc = NULL,
                                 weights = NULL,
                                 derivatives = NULL,
                                 crs = NULL,
                                 ...) {
  fm_basis_mesh_2d(
    mesh = mesh,
    loc = loc,
    weights = weights,
    derivatives = derivatives,
    crs = crs,
    ...
  )
}

#' @export
#' @rdname fm_evaluator_helpers
fm_evaluator_mesh_1d <- function(mesh,
                                 loc,
                                 weights = NULL,
                                 derivatives = NULL,
                                 ...) {
  fm_basis_mesh_1d(
    mesh = mesh,
    loc = loc,
    weights = weights,
    derivatives = derivatives,
    ...
  )
}

#' @title Internal helper functions for mesh field evaluation
#'
#' @description Methods called internally by [fm_basis()] methods.
#' @param weights Optional weight vector, one weight for each location
#' @param derivatives logical; If true, also return matrices `dA` and `d2A`
#' for `fm_mesh_1d` objects, and `dx`, `dy`, `dz` for `fm_mesh_2d`.
#' @inheritParams fm_basis
#' @export
#' @keywords internal
#' @returns A `fm_basis` object; a list of evaluator information objects,
#' at least a matrix `A` and logical vector `ok`.
#' @name fm_basis_helpers
#' @examples
#' str(fm_basis_mesh_2d(fmexample$mesh, loc = fmexample$loc))
#'
fm_basis_mesh_2d <- function(mesh,
                             loc = NULL,
                             weights = NULL,
                             derivatives = NULL,
                             crs = NULL,
                             ...) {
  smorg <- fm_bary(mesh, loc = loc, crs = crs)
  ti <- matrix(0L, NROW(loc), 1)
  ti[, 1L] <- smorg$t
  b <- smorg$bary

  ok <- !is.na(ti[, 1L])

  if (is.null(weights)) {
    weights <- rep(1.0, NROW(loc))
  } else if (length(weights) == 1) {
    weights <- rep(weights, NROW(loc))
  }

  ii <- which(ok)
  A <- (Matrix::sparseMatrix(
    dims = c(NROW(loc), mesh$n),
    i = rep(ii, 3),
    j = as.vector(mesh$graph$tv[ti[ii, 1L], ]),
    x = as.numeric(as.vector(b[ii, ]) * weights[rep(ii, 3)])
  ))

  mesh_deriv <- function(mesh, info, weights) {
    n.mesh <- mesh$n

    ii <- which(info$ok)
    n.ok <- sum(info$ok)
    tv <- mesh$graph$tv[info$t[ii, 1L], , drop = FALSE]
    e1 <- mesh$loc[tv[, 3], , drop = FALSE] - mesh$loc[tv[, 2], , drop = FALSE]
    e2 <- mesh$loc[tv[, 1], , drop = FALSE] - mesh$loc[tv[, 3], , drop = FALSE]
    e3 <- mesh$loc[tv[, 2], , drop = FALSE] - mesh$loc[tv[, 1], , drop = FALSE]
    n1 <- e2 - e1 * matrix(rowSums(e1 * e2) / rowSums(e1 * e1), n.ok, 3)
    n2 <- e3 - e2 * matrix(rowSums(e2 * e3) / rowSums(e2 * e2), n.ok, 3)
    n3 <- e1 - e3 * matrix(rowSums(e3 * e1) / rowSums(e3 * e3), n.ok, 3)
    g1 <- n1 / matrix(rowSums(n1 * n1), n.ok, 3)
    g2 <- n2 / matrix(rowSums(n2 * n2), n.ok, 3)
    g3 <- n3 / matrix(rowSums(n3 * n3), n.ok, 3)
    x <- cbind(g1[, 1], g2[, 1], g3[, 1])
    y <- cbind(g1[, 2], g2[, 2], g3[, 2])
    z <- cbind(g1[, 3], g2[, 3], g3[, 3])
    dx <- (Matrix::sparseMatrix(
      dims = c(nrow(loc), n.mesh),
      i = rep(ii, 3),
      j = as.vector(tv),
      x = as.vector(x) * weights[rep(ii, 3)]
    ))
    dy <- (Matrix::sparseMatrix(
      dims = c(nrow(loc), n.mesh),
      i = rep(ii, 3),
      j = as.vector(tv),
      x = as.vector(y) * weights[rep(ii, 3)]
    ))
    dz <- (Matrix::sparseMatrix(
      dims = c(nrow(loc), n.mesh),
      i = rep(ii, 3),
      j = as.vector(tv),
      x = as.vector(z) * weights[rep(ii, 3)]
    ))

    return(list(dx = dx, dy = dy, dz = dz))
  }

  info <- list(t = ti, bary = b, A = A, ok = ok)

  if (!is.null(derivatives) && derivatives) {
    info <-
      c(
        info,
        mesh_deriv(
          mesh = mesh,
          info = info,
          weights = weights
        )
      )
  }

  structure(
    info,
    class = "fm_basis"
  )
}


#' @param method character; either "default", "nearest", "linear", or
#' "quadratic". With `NULL` or "default", uses the object definition of the
#' function space. Otherwise overrides the object definition.
#' @export
#' @rdname fm_basis_helpers
fm_basis_mesh_1d <- function(mesh,
                             loc,
                             weights = NULL,
                             derivatives = NULL,
                             method = deprecated(),
                             ...) {
  if (lifecycle::is_present(method)) {
    lifecycle::deprecate_soft(
      "0.0.9.9020",
      "fm_evaluator_mesh_1d(method)",
      details = c("Create a separate fm_mesh_1d() object instead.")
    )
    method <- match.arg(method, c(
      "default",
      "nearest",
      "linear",
      "quadratic"
    ))

    if (!(method %in% "default") &&
        (mesh$degree != c(nearest = 0, linear = 1, quadratic = 2)[method])) {
      deg <- c(nearest = 0, linear = 1, quadratic = 2)[method]
      info <- fm_basis_mesh_1d(
        fm_mesh_1d(mesh$loc,
                   interval = mesh$interval,
                   boundary = mesh$boundary,
                   free.clamped = mesh$free.clamped,
                   degree = deg
        ),
        loc = loc,
        weights = weights,
        derivatives = derivatives
      )
      return(info)
    }
  }

  if (is.null(weights)) {
    weights <- rep(1.0, NROW(loc))
  } else if (length(weights) == 1L) {
    weights <- rep(weights, NROW(loc))
  }

  derivatives <- !is.null(derivatives) && derivatives
  info <- list()

  ## Compute basis based on mesh$degree and mesh$boundary
  if (mesh$degree == 0) {
    info <- fm_bary(mesh, loc = loc, method = "nearest")
    i_ <- seq_along(loc)
    j_ <- info$t[, 1]
    x_ <- info$bary[, 1]
    if (derivatives) {
      if (mesh$cyclic) {
        j_prev <- (j_ - 2L) %% mesh$n + 1L
        j_next <- j_ %% mesh$n + 1L
        ok <- rep(TRUE, length(j_))
        dist <- (mesh$loc[j_next] - mesh$loc[j_prev]) %% diff(mesh$interval)
      } else {
        j_prev <- j_ - 1L
        j_next <- j_ + 1L
        ok <- (j_prev >= 1L) & (j_next <= mesh$n)
        dist <- mesh$loc[j_next] - mesh$loc[j_prev]
      }
      i_d <- c(i_[ok], i_[ok])
      j_d <- c(j_prev[ok], j_next[ok])
      x_d <- c(-x_[ok], x_[ok]) / dist
    }

    if (mesh$boundary[1] == "dirichlet") {
      ok <- j_ > 1L
      i_ <- i_[ok]
      j_ <- j_[ok] - 1L
      x_ <- x_[ok]
      if (derivatives) {
        ok <- j_d > 1L
        i_d <- i_d[ok]
        j_d <- j_d[ok] - 1L
        x_d <- x_d[ok]
      }
    }
    if (mesh$boundary[2] == "dirichlet") {
      ok <- j_ <= mesh$m
      i_ <- i_[ok]
      j_ <- j_[ok]
      x_ <- x_[ok]
      if (derivatives) {
        ok <- j_d <= mesh$m
        i_d <- i_d[ok]
        j_d <- j_d[ok]
        x_d <- x_d[ok]
      }
    }
  } else if (mesh$degree == 1) {
    info <- fm_bary(mesh, loc = loc, method = "linear")
    i_ <- c(seq_along(loc), seq_along(loc))
    j_ <- as.vector(info$t)
    x_ <- as.vector(info$bary)
    if (derivatives) {
      j_curr <- info$t[, 1]
      j_next <- info$t[, 2]
      if (mesh$cyclic) {
        if (mesh$n > 1) {
          dist <- (mesh$loc[j_next] - mesh$loc[j_curr]) %% diff(mesh$interval)
        } else {
          dist <- rep(diff(mesh$interval), length(j_curr))
        }
      } else {
        dist <- mesh$loc[j_next] - mesh$loc[j_curr]
      }
      i_d <- i_
      j_d <- c(j_curr, j_next)
      x_d <- rep(c(-1, 1), each = length(loc)) / rep(dist, times = 2)
    }

    if (mesh$boundary[1] == "dirichlet") {
      ok <- j_ > 1L
      i_ <- i_[ok]
      j_ <- j_[ok] - 1L
      x_ <- x_[ok]
      if (derivatives) {
        ok <- j_d > 1L
        i_d <- i_d[ok]
        j_d <- j_d[ok] - 1L
        x_d <- x_d[ok]
      }
    } else if (mesh$boundary[1] == "neumann") {
      if (derivatives) {
        x_d[(j_ == 1) & (x_ > 1)] <- 0.0
        x_d[(j_ == 2) & (x_ < 0)] <- 0.0
      }
      # Set Anew[, 1] = 1 on the left
      # Set Anew[, 2] = 0 on the left
      x_[(j_ == 1) & (x_ > 1)] <- 1.0
      x_[(j_ == 2) & (x_ < 0)] <- 0.0
    }
    if (mesh$boundary[2] == "dirichlet") {
      ok <- j_ <= mesh$m
      i_ <- i_[ok]
      j_ <- j_[ok]
      x_ <- x_[ok]
      if (derivatives) {
        ok <- j_d <= mesh$m
        i_d <- i_d[ok]
        j_d <- j_d[ok]
        x_d <- x_d[ok]
      }
    } else if (mesh$boundary[2] == "neumann") {
      if (derivatives) {
        x_d[(j_ == mesh$m) & (x_ > 1)] <- 0.0
        x_d[(j_ == mesh$m - 1L) & (x_ < 0)] <- 0.0
      }
      # Set Anew[, m] = 1 on the right
      # Set Anew[, m-1] = 0 on the right
      x_[(j_ == mesh$m) & (x_ > 1)] <- 1.0
      x_[(j_ == mesh$m - 1L) & (x_ < 0)] <- 0.0
    }
  } else if (mesh$degree == 2) {
    if (mesh$cyclic) {
      knots <- mesh$loc - mesh$loc[1]
      loc <- loc - mesh$loc[1]
      inter <- c(0, diff(mesh$interval))
    } else {
      knots <- mesh$loc - mesh$loc[1]
      loc <- loc - mesh$loc[1]
      inter <- range(knots)
    }

    info <-
      fm_bary(
        fm_mesh_1d(
          knots,
          interval = inter,
          boundary = if (isTRUE(mesh$cyclic)) "cyclic" else "free",
          degree = 1
        ),
        loc = loc,
        method = "linear"
      )

    if (mesh$cyclic) {
      d <-
        (knots[c(seq_len(length(knots) - 1L) + 1L, 1)] - knots) %%
        diff(mesh$interval)
      d2 <- (knots[c(seq_len(length(knots) - 2L) + 2L, seq_len(2))] -
               knots) %% diff(mesh$interval)
      d2[d2 == 0] <- diff(mesh$interval)
      d <- d[c(length(d), seq_len(length(d) - 1L))]
      d2 <- d2[c(length(d2), seq_len(length(d2) - 1L))]
    } else {
      d <- knots[-1] - knots[-length(knots)]
      d2 <- knots[-seq_len(2)] - knots[seq_len(length(knots) - 2L)]
    }

    if (mesh$cyclic) {
      ## Left intervals for each basis function:
      i.l <- seq_along(info$t[, 1])
      j.l <- info$t[, 1] + 2L
      x.l <- (info$bary[, 2] * d[info$t[, 2]] / d2[info$t[, 2]] * info$bary[, 2])
      if (derivatives) {
        x.d1.l <- (2 / d2[info$t[, 2]] * info$bary[, 2])
        x.d2.l <- (2 / d2[info$t[, 2]] / d[info$t[, 2]])
      }
      ## Right intervals for each basis function:
      i.r <- seq_along(info$t[, 1])
      j.r <- info$t[, 1]
      x.r <- (info$bary[, 1] * d[info$t[, 2]] / d2[info$t[, 1]] * info$bary[, 1])
      if (derivatives) {
        x.d1.r <- -(2 / d2[info$t[, 2]] * info$bary[, 1])
        x.d2.r <- (2 / d2[info$t[, 1]] / d[info$t[, 2]])
      }
      ## Middle intervals for each basis function:
      i.m <- seq_along(info$t[, 1])
      j.m <- info$t[, 1] + 1L
      x.m <- (1 - (info$bary[, 1] * d[info$t[, 2]] / d2[info$t[, 1]] * info$bary[, 1] +
                     info$bary[, 2] * d[info$t[, 2]] / d2[info$t[, 2]] * info$bary[, 2]
      ))
      if (derivatives) {
        x.d1.m <- (2 / d2[info$t[, 1]] * info$bary[, 1]) -
          (2 / d2[info$t[, 2]] * info$bary[, 2])
        x.d2.m <- -(2 / d2[info$t[, 1]] / info$t[, 2]) -
          (2 / d2[info$t[, 2]] / info$t[, 2])
      }
    } else {
      d2 <- c(2 * d[1], 2 * d[1], d2, 2 * d[length(d)], 2 * d[length(d)])
      d <- c(d[1], d[2], d, d[length(d)], d[length(d)])
      ok <- (info$t[, 1] >= 1L) & (info$t[, 2] <= length(knots))
      index <- info$t[ok, , drop = FALSE] + 1L
      bary <- info$bary[ok, , drop = FALSE]
      ## Left intervals for each basis function:
      i.l <- seq_along(loc)[ok]
      j.l <- index[, 1] + 1L
      x.l <- (bary[, 2] * d[index[, 2]] / d2[index[, 2]] * bary[, 2])
      if (derivatives) {
        x.d1.l <- (2 / d2[index[, 2]] * bary[, 2])
        x.d2.l <- (2 / d2[index[, 2]] / d[index[, 2]])
      }
      ## Right intervals for each basis function:
      i.r <- seq_along(loc)[ok]
      j.r <- index[, 1] - 1L
      x.r <- (bary[, 1] * d[index[, 2]] / d2[index[, 1]] * bary[, 1])
      if (derivatives) {
        x.d1.r <- -(2 / d2[index[, 1]] * bary[, 1])
        x.d2.r <- (2 / d2[index[, 1]] / d[index[, 2]])
      }
      ## Middle intervals for each basis function:
      i.m <- seq_along(loc)[ok]
      j.m <- index[, 1]
      x.m <- (1 - (bary[, 1] * d[index[, 2]] / d2[index[, 1]] * bary[, 1] +
                     bary[, 2] * d[index[, 2]] / d2[index[, 2]] * bary[, 2]
      ))
      if (derivatives) {
        x.d1.m <- (2 / d2[index[, 1]] * bary[, 1]) - (2 / d2[index[, 2]] * bary[, 2])
        x.d2.m <- -(2 / d2[index[, 1]] / d[index[, 2]]) -
          (2 / d2[index[, 2]] / d[index[, 2]])
      }
    }

    i_ <- c(i.l, i.r, i.m)
    j_ <- c(j.l, j.r, j.m)
    x_ <- c(x.l, x.r, x.m)
    if (derivatives) {
      x_d1 <- c(x.d1.l, x.d1.r, x.d1.m)
      x_d2 <- c(x.d2.l, x.d2.r, x.d2.m)
    }

    if (!mesh$cyclic) {
      # Convert boundary basis functions to linear
      # First remove anything from above outside the interval, then add back in
      # the appropriate values
      ok <- (loc >= inter[1]) & (loc <= inter[2])
      i_ <- i_[ok]
      j_ <- j_[ok]
      x_ <- x_[ok]
      if (derivatives) {
        x_d1 <- x_d1[ok]
        x_d2 <- x_d2[ok]
      }

      # left
      ok <- (loc < 0) & (info$t[, 1] == 1L)
      i_l <- c(seq_along(loc)[ok], seq_along(loc)[ok])
      j_l <- c(info$t[ok, 1], info$t[ok, 2])
      x_l <- c(
        0.5 + (info$bary[ok, 1] - 1),
        0.5 - (info$bary[ok, 1] - 1)
      )
      if (derivatives) {
        x_d1_l <- rep(c(-1, 1) / d[1], each = sum(ok))
        x_d2_l <- rep(c(0, 0), each = sum(ok))
      }

      # right
      ok <- (loc > inter[2]) & (info$t[, 2] == length(knots))
      i_r <- c(seq_along(loc)[ok], seq_along(loc)[ok])
      j_r <- c(info$t[ok, 2], info$t[ok, 1]) + 1L
      x_r <- c(
        0.5 + (info$bary[ok, 2] - 1),
        0.5 - (info$bary[ok, 2] - 1)
      )
      if (derivatives) {
        x_d1_r <- rep(c(1, -1) / d[length(d)], each = sum(ok))
        x_d2_r <- rep(c(0, 0), each = sum(ok))
      }

      i_ <- c(i_, i_l, i_r)
      j_ <- c(j_, j_l, j_r)
      x_ <- c(x_, x_l, x_r)
      if (derivatives) {
        x_d1 <- c(x_d1, x_d1_l, x_d1_r)
        x_d2 <- c(x_d2, x_d2_l, x_d2_r)
      }

      if (mesh$boundary[1] == "dirichlet") {
        ok <- j_ > 1L
        j_[ok] <- j_[ok] - 1L
        x_[!ok] <- -x_[!ok]
        if (derivatives) {
          x_d1[!ok] <- -x_d1[!ok]
          x_d2[!ok] <- -x_d2[!ok]
        }
      } else if (mesh$boundary[1] == "neumann") {
        ok <- j_ > 1L
        j_[ok] <- j_[ok] - 1L
      } else if ((mesh$boundary[1] == "free") &&
                 (mesh$free.clamped[1])) {
        # new1 <- 2 * basis1
        # new2 <- basis2 - basis1
        ok1 <- j_ == 1L
        i2 <- i_[ok1]
        j2 <- j_[ok1] + 1L
        x2 <- -x_[ok1]
        x_[ok1] <- 2 * x_[ok1]
        if (derivatives) {
          x2_d1 <- -x_d1[ok1]
          x_d1[ok1] <- 2 * x_d1[ok1]
          x2_d2 <- -x_d2[ok1]
          x_d2[ok1] <- 2 * x_d2[ok1]
        }
        i_ <- c(i_, i2)
        j_ <- c(j_, j2)
        x_ <- c(x_, x2)
        if (derivatives) {
          x_d1 <- c(x_d1, x2_d1)
          x_d2 <- c(x_d2, x2_d1)
        }
      }
      if (mesh$boundary[2] == "dirichlet") {
        ok <- j_ > mesh$m
        j_[ok] <- j_[ok] - 1L
        x_[ok] <- -x_[ok]
        if (derivatives) {
          x_d1[ok] <- -x_d1[ok]
          x_d2[ok] <- -x_d2[ok]
        }
      } else if (mesh$boundary[2] == "neumann") {
        ok <- j_ > mesh$m
        j_[ok] <- mesh$m
      } else if ((mesh$boundary[2] == "free") &&
                 (mesh$free.clamped[2])) {
        # new_m <- m + {m-1};     m = 1, m - 1 = 2
        # new_{m-1} <- {m-1} - m; m = 1, m - 1 = 2
        # new1 <- 2 * basis1
        # new2 <- basis2 - basis1
        ok1 <- j_ == mesh$m
        i2 <- i_[ok1]
        j2 <- j_[ok1] - 1L
        x2 <- -x_[ok1]
        x_[ok1] <- 2 * x_[ok1]
        if (derivatives) {
          x2_d1 <- -x_d1[ok1]
          x_d1[ok1] <- 2 * x_d1[ok1]
          x2_d2 <- -x_d2[ok1]
          x_d2[ok1] <- 2 * x_d2[ok1]
        }
        i_ <- c(i_, i2)
        j_ <- c(j_, j2)
        x_ <- c(x_, x2)
        if (derivatives) {
          x_d1 <- c(x_d1, x2_d1)
          x_d2 <- c(x_d2, x2_d1)
        }
      }
    }

    if (mesh$cyclic) {
      j_ <- (j_ - 1L - 1L) %% mesh$m + 1L
    }
  } else {
    stop("Unsupported B-spline degree = ", mesh$degree)
  }

  info$A <- Matrix::sparseMatrix(
    i = i_,
    j = j_,
    x = (weights[i_] * x_),
    dims = c(length(loc), mesh$m)
  )
  if (derivatives) {
    info$dA <- Matrix::sparseMatrix(
      i = i_,
      j = j_,
      x = weights[i_] * x_d1,
      dims = c(length(loc), mesh$m)
    )
    info$d2A <- Matrix::sparseMatrix(
      i = i_,
      j = j_,
      x = weights[i_] * x_d2,
      dims = c(length(loc), mesh$m)
    )
  }

  info[["ok"]] <- rep(TRUE, length(loc))

  structure(
    info,
    class = "fm_basis"
  )
}







#' @export
#' @describeIn fm_evaluate The `...` arguments are passed on to `fm_evaluator_lattice()`
#' if no `loc` or `lattice` is provided.
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
      structure(list(
        x = x,
        y = y,
        lattice = lattice,
        loc = NULL,
        proj = proj,
        crs = crs
      ),
      class = "fm_evaluator")
  } else {
    proj <- fm_basis_mesh_2d(mesh, loc = loc, crs = crs)
    projector <-
      structure(list(
        x = NULL,
        y = NULL,
        lattice = NULL,
        loc = loc,
        proj = proj,
        crs = crs
      ),
      class = "fm_evaluator")
  }

  return(projector)
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

  return(projector)
}



#' @describeIn fm_evaluate
#' Creates an [fm_lattice_2d()] object, by default covering the input mesh.
#' @export
fm_evaluator_lattice <- function(mesh,
                                 xlim = NULL,
                                 ylim = NULL,
                                 dims = c(100, 100),
                                 projection = NULL,
                                 crs = NULL,
                                 ...) {
  stopifnot(inherits(mesh, c("fm_mesh_2d", "inla.mesh")))
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


#' @export
#' @describeIn fm_evaluate Converts legacy `inla.mesh` to `fm_mesh_2d` and calls
#' the `fm_evaluator` method again.
fm_evaluator.inla.mesh <- function(mesh, ...) {
  fm_evaluator(fm_as_mesh_2d(mesh), ...)
}
#' @export
#' @describeIn fm_evaluate Converts legacy `inla.mesh` to `fm_mesh_1d` and calls
#' the `fm_evaluator` method again.
fm_evaluator.inla.mesh.1d <- function(mesh, ...) {
  fm_evaluator(fm_as_mesh_1d(mesh), ...)
}





# fm_contains ####

#' Check which mesh triangles are inside a polygon
#'
#' Wrapper for the [sf::st_contains()] (previously `sp::over()`) method to find triangle centroids
#' or vertices inside `sf` or `sp` polygon objects
#'
#' @param x geometry (typically an `sf` or `sp::SpatialPolygons` object) for the queries
#' @param y an [fm_mesh_2d()] or `inla.mesh` object
#' @param \dots Passed on to other methods
#' @param type the query type; either `'centroid'` (default, for triangle centroids),
#' or `'vertex'` (for mesh vertices)
#'
#' @return List of vectors of triangle indices (when `type` is `'centroid'`) or
#' vertex indices (when `type` is `'vertex'`). The list has one entry per row of the `sf` object.
#' Use `unlist(fm_contains(...))` if the combined union is needed.
#'
#' @author Haakon Bakka, \email{bakka@@r-inla.org}, and Finn Lindgren \email{finn.lindgren@@gmail.com}
#'
#' @examples
#' if (TRUE &&
#'   fm_safe_sp()) {
#'   # Create a polygon and a mesh
#'   obj <- sp::SpatialPolygons(
#'     list(sp::Polygons(
#'       list(sp::Polygon(rbind(
#'         c(0, 0),
#'         c(50, 0),
#'         c(50, 50),
#'         c(0, 50)
#'       ))),
#'       ID = 1
#'     )),
#'     proj4string = fm_CRS("longlat_globe")
#'   )
#'   mesh <- fm_rcdt_2d_inla(globe = 2, crs = fm_crs("sphere"))
#'
#'   ## 3 vertices found in the polygon
#'   fm_contains(obj, mesh, type = "vertex")
#'
#'   ## 3 triangles found in the polygon
#'   fm_contains(obj, mesh)
#'
#'   ## Multiple transformations can lead to slightly different results due to edge cases
#'   ## 4 triangles found in the polygon
#'   fm_contains(
#'     obj,
#'     fm_transform(mesh, crs = fm_crs("mollweide_norm"))
#'   )
#' }
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
  if (!inherits(y, c("fm_mesh_2d", "inla.mesh"))) {
    stop(paste0(
      "'y' must be an 'fm_mesh_2d' or 'inla.mesh' object, not '",
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




# fm_is_within ####

#' @title Query if points are inside a mesh
#'
#' @description
#'  Queries whether each input point is within a mesh or not.
#'
#' @param x A set of points of a class supported by `fm_evaluator(y, loc = x)`
#' @param y An `inla.mesh`
#' @param \dots Currently unused
#' @returns A logical vector
#' @examples
#' all(fm_is_within(fmexample$loc, fmexample$mesh))
#' @export
fm_is_within <- function(x, y, ...) {
  UseMethod("fm_is_within")
}

#' @rdname fm_is_within
#' @export
fm_is_within.default <- function(x, y, ...) {
  fm_evaluator(y, loc = x)$proj$ok
}



# fm_basis ####

#' @title Compute mapping matrix between mesh function space and points
#'
#' @description
#'  Computes the basis mapping matrix between a function space on a mesh, and locations.
#'
#' @param x An object supported by the [fm_evaluator()] class
#' @param loc A set of points of a class supported by `fm_evaluator(x, loc = loc)`
#' @param full logical; if `TRUE`, return a `fm_basis` object, containing at least
#' a projection matrix `A` and logical vector `ok` indicating which evaluations
#' are valid. If `FALSE`, return only the projection matrix `A`. Default is `FALSE`.
#' @param \dots Currently unused
#' @returns A `sparseMatrix` object (if `full = FALSE`),
#' or a `fm_basis` object (if `full = TRUE` or `isTRUE(derivatives)`).
#' @seealso [fm_raw_basis()]
#' @examples
#' # Compute basis mapping matrix
#' str(fm_basis(fmexample$mesh, fmexample$loc))
#' print(fm_basis(fmexample$mesh, fmexample$loc), full = TRUE)
#' @export
fm_basis <- function(x, ..., full = FALSE) {
  UseMethod("fm_basis")
}

#' Print method for `fm_basis()`
#'
#' Prints information for an [fm_basis()] object.
#'
#' @param x [fm_basis()] object
#' @param \dots Unused
#' @returns `invisible(x)`
#' @seealso [fm_basis()]
#' @export
#' @examples
#' print(fm_basis(fmexample$mesh, fmexample$loc, full = TRUE))
print.fm_basis <- function(x, ...) {
  cat("fm_basis object\n")
  cat("  Projection matrix (A): ", paste0(dim(x$A), collapse = "-by-"), "\n", sep = "")
  cat("  Valid evaluations (ok): ", sum(x$ok), " out of ", length(x$ok), "\n", sep = "")
  cat("  Additional information: ",
      paste(names(x)[!names(x) %in% c("A", "ok")], collapse = ", "),
      "\n",
      sep = "")

  invisible(x)
}

#' @rdname fm_basis
#' @export
fm_basis.default <- function(x, ..., full = FALSE) {
  lifecycle::deprecate_stop(
    "0.1.7.9002",
    "fm_basis.default()",
    details = "Each mesh class needs its own `fm_basis()` method.")
}

#' @param weights Optional weight vector to apply (from the left, one
#' weight for each row of the basis matrix)
#' @param derivatives If non-NULL and logical, include derivative matrices
#' in the output. Forces `full = TRUE`.
#' \item{A }{The projection matrix, `u(loc_i)=sum_j A_ij w_i`}
#' \item{d1A, d2A }{Derivative weight matrices,
#' `du/dx(loc_i)=sum_j dx_ij w_i`, etc.}
#' @rdname fm_basis
#' @export
fm_basis.fm_mesh_1d <- function(x, loc, weights = NULL, derivatives = NULL, ..., full = FALSE) {
  result <- fm_basis_mesh_1d(
    x,
    loc = loc,
    weights = weights,
    derivatives = derivatives,
    ...
  )
  if (isTRUE(derivatives) && !full) {
    full <- TRUE
  }
  fm_basis(result, full = full)
}

#' @return For `fm_mesh_2d`, a matrix, or if `derivatives` is `TRUE`,
#' a list with elements
#' \item{A }{The projection matrix, `u(loc_i)=sum_j A_ij w_i`}
#' \item{dx, dy, dz }{Derivative weight matrices, `du/dx(loc_i)=sum_j
#' dx_ij w_i`, etc.}
#' @rdname fm_basis
#' @export
fm_basis.fm_mesh_2d <- function(x, loc, weights = NULL, derivatives = NULL, ...,
                                full = FALSE) {
  result <- fm_basis_mesh_2d(
    x,
    loc = loc,
    weights = weights,
    derivatives = derivatives,
    ...
  )
  if (isTRUE(derivatives) && !full) {
    full <- TRUE
  }
  fm_basis(result, full = full)
}

#' @rdname fm_basis
#' @export
#' @method fm_basis inla.mesh.1d
fm_basis.inla.mesh.1d <- function(x, ...) {
  fm_basis.fm_mesh_1d(fm_as_mesh_1d(x), ...)
}

#' @rdname fm_basis
#' @export
#' @method fm_basis inla.mesh
fm_basis.inla.mesh <- function(x, ...) {
  fm_basis.fm_mesh_2d(fm_as_mesh_2d(x), ...)
}

#' @rdname fm_basis
#' @export
fm_basis.fm_evaluator <- function(x, ..., full = FALSE) {
  fm_basis(
    structure(
      x$proj,
      class = "fm_basis"
      ),
    full = full)
}

#' @rdname fm_basis
#' @export
fm_basis.fm_basis <- function(x, ..., full = FALSE) {
  if (full) {
    x
  } else {
    x[["A"]]
  }
}

#' @param x [fm_tensor()] object
#' @export
#' @rdname fm_basis
fm_basis.fm_tensor <- function(x,
                               loc,
                               weights = NULL,
                               ...,
                               full = FALSE) {
  if (length(loc) != length(x[["fun_spaces"]])) {
    stop(
      paste0(
        "Length of location list (",
        length(loc), ") doesn't match the number of function spaces (",
        length(x[["fun_spaces"]]),
        ")"
      )
    )
  }
  if (is.null(names(loc))) {
    names(loc) <- names(x[["fun_spaces"]])
  } else if (!setequal(names(x[["fun_spaces"]]), names(loc))) {
    stop("Name mismatch between location list names and function space names.")
  }
  idx <- names(x[["fun_spaces"]])
  if (is.null(idx)) {
    idx <- seq_along(x[["fun_spaces"]])
  }
  proj <- lapply(
    idx,
    function(k) {
      fm_basis(x[["fun_spaces"]][[k]], loc = loc[[k]], full = TRUE)
    }
  )
  names(proj) <- names(x[["fun_spaces"]])

  # Combine the matrices
  # (A1, A2, A3) -> rowkron(A3, rowkron(A2, A1))
  A <- proj[[1]][["A"]]
  if (!is.null(weights)) {
    A <- Matrix::Diagonal(nrow(A), x = weights) %*% A
  }
  ok <- proj[[1]][["ok"]]
  for (k in seq_len(length(x[["fun_spaces"]]) - 1)) {
    A <- fm_row_kron(proj[[k + 1]][["A"]], A)
    ok <- proj[[k + 1]][["ok"]] & ok
  }

  out <- structure(
    list(A = A, ok = ok),
    class = "fm_basis"
  )
  fm_basis(out, full = full)
}










internal_spline_mesh_1d <- function(interval, m, degree, boundary, free.clamped) {
  boundary <-
    match.arg(
      boundary,
      c("neumann", "dirichlet", "free", "cyclic")
    )
  if (degree <= 1) {
    n <- (switch(boundary,
      neumann = m,
      dirichlet = m + 2,
      free = m,
      cyclic = m + 1
    ))
    if (n < 2) {
      n <- 2
      degree <- 0
      boundary <- "c"
    }
  } else {
    stopifnot(degree == 2)
    n <- (switch(boundary,
      neumann = m + 1,
      dirichlet = m + 1,
      free = m - 1,
      cyclic = m
    ))
    if (boundary == "free") {
      if (m <= 1) {
        n <- 2
        degree <- 0
        boundary <- "c"
      } else if (m == 2) {
        n <- 2
        degree <- 1
      }
    } else if (boundary == "cyclic") {
      if (m <= 1) {
        n <- 2
        degree <- 0
      }
    }
  }
  return(fm_mesh_1d(seq(interval[1], interval[2], length.out = n),
    degree = degree,
    boundary = boundary,
    free.clamped = free.clamped
  ))
}


#' Basis functions for mesh manifolds
#'
#' Calculate basis functions on [fm_mesh_1d()] or [fm_mesh_2d()],
#' without necessarily matching the default function space of the given mesh
#' object.
#'
#' @param mesh An [fm_mesh_1d()] or [fm_mesh_2d()] object.
#' @param type `b.spline` (default) for B-spline basis functions,
#' `sph.harm` for spherical harmonics (available only for meshes on the
#' sphere)
#' @param n For B-splines, the number of basis functions in each direction (for
#' 1d meshes `n` must be a scalar, and for planar 2d meshes a 2-vector).
#' For spherical harmonics, `n` is the maximal harmonic order.
#' @param degree Degree of B-spline polynomials.  See
#' [fm_mesh_1d()].
#' @param knot.placement For B-splines on the sphere, controls the latitudinal
#' placements of knots. `"uniform.area"` (default) gives uniform spacing
#' in `sin(latitude)`, `"uniform.latitude"` gives uniform spacing in
#' latitudes.
#' @param rot.inv For spherical harmonics on a sphere, `rot.inv=TRUE`
#' gives the rotationally invariant subset of basis functions.
#' @param boundary Boundary specification, default is free boundaries.  See
#' [fm_mesh_1d()] for more information.
#' @param free.clamped If `TRUE` and `boundary` is `"free"`, the
#' boundary basis functions are clamped to 0/1 at the interval boundary by
#' repeating the boundary knots. See
#' [fm_mesh_1d()] for more information.
#' @param ... Unused
#' @returns A matrix with evaluated basis function
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @seealso [fm_mesh_1d()], [fm_mesh_2d()], [fm_basis()]
#' @examples
#'
#' loc <- rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1))
#' mesh <- fm_mesh_2d(loc, max.edge = 0.15)
#' basis <- fm_raw_basis(mesh, n = c(4, 5))
#'
#' proj <- fm_evaluator(mesh, dims = c(10, 10))
#' image(proj$x, proj$y, fm_evaluate(proj, basis[, 7]), asp = 1)
#' \donttest{
#' if (interactive() && require("rgl")) {
#'   plot_rgl(mesh, col = basis[, 7], draw.edges = FALSE, draw.vertices = FALSE)
#' }
#' }
#'
#' @export
fm_raw_basis <- function(mesh,
                         type = "b.spline",
                         n = 3,
                         degree = 2,
                         knot.placement = "uniform.area",
                         rot.inv = TRUE,
                         boundary = "free",
                         free.clamped = TRUE,
                         ...) {
  type <- match.arg(type, c("b.spline", "sph.harm"))
  knot.placement <- (match.arg(
    knot.placement,
    c(
      "uniform.area",
      "uniform.latitude"
    )
  ))

  if (identical(type, "b.spline")) {
    if (fm_manifold(mesh, c("R1", "S1"))) {
      mesh1 <-
        internal_spline_mesh_1d(
          mesh$interval, n, degree,
          boundary, free.clamped
        )
      basis <- fm_basis(mesh1, mesh$loc)
    } else if (identical(mesh$manifold, "R2")) {
      if (length(n) == 1) {
        n <- rep(n, 2)
      }
      if (length(degree) == 1) {
        degree <- rep(degree, 2)
      }
      if (length(boundary) == 1) {
        boundary <- rep(boundary, 2)
      }
      if (length(free.clamped) == 1) {
        free.clamped <- rep(free.clamped, 2)
      }
      mesh1x <-
        internal_spline_mesh_1d(
          range(mesh$loc[, 1]),
          n[1], degree[1],
          boundary[1], free.clamped[1]
        )
      mesh1y <-
        internal_spline_mesh_1d(
          range(mesh$loc[, 2]),
          n[2], degree[2],
          boundary[2], free.clamped[2]
        )
      basis <-
        fm_row_kron(
          fm_basis(mesh1y, mesh$loc[, 2]),
          fm_basis(mesh1x, mesh$loc[, 1])
        )
    } else if (identical(mesh$manifold, "S2")) {
      loc <- mesh$loc
      uniform.lat <- identical(knot.placement, "uniform.latitude")
      degree <- max(0L, min(n - 1L, degree))
      basis <- fmesher_spherical_bsplines1(
        loc[, 3],
        n = n,
        degree = degree,
        uniform = uniform.lat
      )
      if (!rot.inv) {
        warning("Currently only 'rot.inv=TRUE' is supported for B-splines.")
      }
    } else {
      stop("Only know how to make B-splines on R2 and S2.")
    }
  } else if (identical(type, "sph.harm")) {
    if (!identical(mesh$manifold, "S2")) {
      stop("Only know how to make spherical harmonics on S2.")
    }
    # With GSL activated:
    #        if (rot.inv) {
    #            basis <- (inla.fmesher.smorg(
    #                mesh$loc,
    #                mesh$graph$tv,
    #                sph0 = n
    #            )$sph0)
    #        } else {
    #            basis <- (inla.fmesher.smorg(
    #                mesh$loc,
    #                mesh$graph$tv,
    #                sph = n
    #            )$sph)
    #        }

    fm_require_stop(
      "gsl",
      "The 'gsl' R package is needed for spherical harmonics."
    )

    # Make sure we have radius-1 coordinates
    loc <- mesh$loc / rowSums(mesh$loc^2)^0.5
    if (rot.inv) {
      basis <- matrix(0, nrow(loc), n + 1)
      for (l in seq(0, n)) {
        basis[, l + 1] <- sqrt(2 * l + 1) *
          gsl::legendre_Pl(l = l, x = loc[, 3])
      }
    } else {
      angle <- atan2(loc[, 2], loc[, 1])
      basis <- matrix(0, nrow(loc), (n + 1)^2)
      for (l in seq(0, n)) {
        basis[, 1 + l * (l + 1)] <-
          sqrt(2 * l + 1) *
            gsl::legendre_Pl(l = l, x = loc[, 3])
        for (m in seq_len(l)) {
          scaling <- sqrt(2 * (2 * l + 1) * exp(lgamma(l - m + 1) - lgamma(l + m + 1)))
          poly <- gsl::legendre_Plm(l = l, m = m, x = loc[, 3])
          basis[, 1 + l * (l + 1) - m] <-
            scaling * sin(-m * angle) * poly
          basis[, 1 + l * (l + 1) + m] <-
            scaling * cos(m * angle) * poly
        }
      }
    }
  }

  return(basis)
}







# Block methods ####

#' Blockwise aggregation matrices
#'
#' Creates an aggregation matrix for blockwise aggregation, with optional
#' weighting.
#'
#' @param block integer vector; block information. If `NULL`,
#' `rep(1L, block_len)` is used, where `block_len` is determined by
#' `length(log_weights)))` or `length(weights)))`.
#' A single scalar is also repeated
#' to a vector of corresponding length to the weights.
#' @param weights Optional weight vector
#' @param log_weights Optional `log(weights)` vector. Overrides `weights` when
#' non-NULL.
#' @param rescale logical; If `TRUE`, normalise the weights by `sum(weights)`
#' or `sum(exp(log_weights))` within each block.
#' Default: `FALSE`
#' @param n_block integer; The number of conceptual blocks. Only needs to be
#' specified if it's larger than `max(block)`, or to keep the output of
#' consistent size for different inputs.
#'
#' @returns A (sparse) matrix
#' @export
#' @describeIn fm_block A (sparse) matrix of size `n_block` times `length(block)`.
#' @examples
#' block <- rep(1:2, 3:2)
#' fm_block(block)
#' fm_block(block, rescale = TRUE)
#' fm_block(block, log_weights = -2:2, rescale = TRUE)
#' fm_block_eval(
#'   block,
#'   weights = 1:5,
#'   rescale = TRUE,
#'   values = 11:15
#' )
#' fm_block_logsumexp_eval(
#'   block,
#'   weights = 1:5,
#'   rescale = TRUE,
#'   values = log(11:15),
#'   log = FALSE
#' )
fm_block <- function(block = NULL,
                     weights = NULL,
                     log_weights = NULL,
                     rescale = FALSE,
                     n_block = NULL) {
  info <-
    fm_block_prep(
      block = block,
      n_block = n_block,
      log_weights = log_weights,
      weights = weights
    )
  if (length(info$block) == 0) {
    return(
      Matrix::sparseMatrix(
        i = integer(0),
        j = integer(0),
        x = numeric(0),
        dims = c(0L, 0L)
      )
    )
  }
  weights <-
    fm_block_weights(
      block = info$block,
      weights = info$weights,
      log_weights = info$log_weights,
      n_block = info$n_block,
      rescale = rescale
    )

  Matrix::sparseMatrix(
    i = info$block,
    j = seq_len(length(info$block)),
    x = as.numeric(weights),
    dims = c(info$n_block, length(info$block))
  )
}

#' @param values Vector to be blockwise aggregated
#' @describeIn fm_block Evaluate aggregation. More efficient alternative to to
#' `as.vector(fm_block(...) %*% values)`.
#' @export
fm_block_eval <- function(block = NULL,
                          weights = NULL,
                          log_weights = NULL,
                          rescale = FALSE,
                          n_block = NULL,
                          values = NULL) {
  info <-
    fm_block_prep(
      block = block,
      n_block = n_block,
      log_weights = log_weights,
      weights = weights,
      values = values
    )
  if (length(info$block) == 0) {
    return(numeric(0L))
  }
  weights <-
    fm_block_weights(
      block = info$block,
      weights = info$weights,
      log_weights = info$log_weights,
      n_block = info$n_block,
      rescale = rescale
    )

  val <-
    Matrix::sparseMatrix(
      i = info$block,
      j = rep(1L, length(info$block)),
      x = as.numeric(values * weights),
      dims = c(info$n_block, 1)
    )
  as.vector(val)
}


#' @param log If `TRUE` (default), return log-sum-exp. If `FALSE`,
#' return sum-exp.
#' @describeIn fm_block Evaluate log-sum-exp aggregation.
#' More efficient and numerically stable alternative to to
#' `log(as.vector(fm_block(...) %*% exp(values)))`.
#' @export
fm_block_logsumexp_eval <- function(block = NULL,
                                    weights = NULL,
                                    log_weights = NULL,
                                    rescale = FALSE,
                                    n_block = NULL,
                                    values = NULL,
                                    log = TRUE) {
  info <-
    fm_block_prep(
      block = block,
      n_block = n_block,
      log_weights = log_weights,
      weights = weights,
      values = values
    )
  if (length(info$block) == 0) {
    return(numeric(0L))
  }
  log_weights <-
    fm_block_log_weights(
      block = info$block,
      weights = info$weights,
      log_weights = info$log_weights,
      n_block = info$n_block,
      rescale = rescale
    )

  # Compute shift for stable log-sum-exp
  w_values <- values + log_weights
  shift <- fm_block_log_shift(
    log_weights = w_values,
    block = info$block,
    n_block = info$n_block
  )

  val <-
    Matrix::sparseMatrix(
      i = info$block,
      j = rep(1L, length(info$block)),
      x = as.numeric(exp(w_values - shift[info$block])),
      dims = c(info$n_block, 1)
    )

  if (log) {
    val <- log(as.vector(val)) + shift
  } else {
    val <- as.vector(val) * exp(shift)
  }

  val
}


#' @describeIn fm_block Computes (optionally) blockwise renormalised weights
#' @export
fm_block_weights <-
  function(block = NULL,
           weights = NULL,
           log_weights = NULL,
           rescale = FALSE,
           n_block = NULL) {
    info <-
      fm_block_prep(
        block = block,
        n_block = n_block,
        log_weights = log_weights,
        weights = weights
      )
    if (length(info$block) == 0) {
      return(numeric(0))
    }
    if (rescale) {
      # Compute blockwise normalised weights
      if (!is.null(info$log_weights)) {
        info$log_weights <- fm_block_log_weights(
          log_weights = info$log_weights,
          block = info$block,
          n_block = info$n_block,
          rescale = rescale
        )
        info$weights <- exp(info$log_weights)
      } else {
        if (is.null(info$weights)) {
          info$weights <- rep(1.0, length(info$block))
        }
        scale <- as.vector(Matrix::sparseMatrix(
          i = info$block,
          j = rep(1L, length(info$block)),
          x = as.numeric(info$weights),
          dims = c(info$n_block, 1L)
        ))
        info$weights <- info$weights / scale[block]
      }
    } else if (!is.null(info$log_weights)) {
      info$weights <- exp(info$log_weights)
    } else if (is.null(info$weights)) {
      info$weights <- rep(1.0, length(info$block))
    }
    info$weights
  }

#' @describeIn fm_block Computes (optionally) blockwise renormalised log-weights
#' @export
fm_block_log_weights <- function(block = NULL,
                                 weights = NULL,
                                 log_weights = NULL,
                                 rescale = FALSE,
                                 n_block = NULL) {
  info <-
    fm_block_prep(
      block = block,
      n_block = n_block,
      log_weights = log_weights,
      weights = weights,
      force_log = TRUE
    )
  if (length(info$block) == 0) {
    return(numeric(0))
  }
  if (is.null(info$log_weights)) {
    if (is.null(info$weights)) {
      info$log_weights <- rep(0.0, length(info$block))
    } else {
      info$log_weights <- log(info$weights)
    }
  }
  if (rescale) {
    shift <- fm_block_log_shift(
      block = info$block, log_weights = info$log_weights,
      n_block = info$n_block
    )
    log_rescale <- as.vector(
      Matrix::sparseMatrix(
        i = info$block,
        j = rep(1L, length(info$block)),
        x = as.numeric(exp(info$log_weights - shift[info$block])),
        dims = c(info$n_block, 1L)
      )
    )
    log_rescale <- (log(log_rescale) + shift)[block]
    info$log_weights <- info$log_weights - log_rescale
  }

  info$log_weights
}


#' @describeIn fm_block Computes shifts for stable blocked log-sum-exp.
#' To compute \eqn{\log(\sum_{i; \textrm{block}_i=k} \exp(v_i) w_i)}{
#' log(sum_(i;block_i=k) exp(v_i) w_i)
#' } for
#' each block `k`, first compute combined values and weights, and a shift:
#' ```
#' w_values <- values + fm_block_log_weights(block, log_weights = log_weights)
#' shift <- fm_block_log_shift(block, log_weights = w_values)
#' ```
#' Then aggregate the values within each block:
#' ```
#' agg <- aggregate(exp(w_values - shift[block]),
#'                  by = list(block = block),
#'                  \(x) log(sum(x)))
#' agg$x <- agg$x + shift[agg$block]
#' ```
#' The implementation uses a faster method:
#' ```
#' as.vector(
#'   Matrix::sparseMatrix(
#'     i = block,
#'     j = rep(1L, length(block)),
#'     x = exp(w_values - shift[block]),
#'     dims = c(n_block, 1))
#' ) + shift
#' ```
#' @export
fm_block_log_shift <- function(block = NULL,
                               log_weights = NULL,
                               n_block = NULL) {
  info <-
    fm_block_prep(
      block = block,
      n_block = n_block,
      log_weights = log_weights,
      force_log = TRUE
    )
  if (length(info$block) == 0) {
    return(0.0)
  }
  block_k <- sort(unique(info$block))
  shift <- numeric(info$n_block)
  if (!is.null(info$log_weights)) {
    if (info$n_block <= 200) {
      # Fast for small n_block
      shift[block_k] <-
        vapply(
          block_k,
          function(k) {
            max(info$log_weights[info$block == k])
          },
          0.0
        )
    } else {
      # Fast for large n_block
      shift[block_k] <-
        stats::aggregate(
          info[["log_weights"]],
          by = list(block = info[["block"]]),
          FUN = max,
          simplify = TRUE
        )$x
    }
  }

  shift
}


#' @describeIn fm_block Helper function for preparing `block`, `weights`, and
#' `log_weights`, `n_block` inputs.
#' @export
#' @param n_values When supplied, used instead of `length(values)` to determine
#' the value vector input length.
#' @param force_log When `FALSE` (default),
#' passes either `weights` and `log_weights` on, if provided, with `log_weights`
#' taking precedence. If `TRUE`, forces the computation of `log_weights`,
#' whether given in the input or not.
fm_block_prep <- function(block = NULL,
                          log_weights = NULL,
                          weights = NULL,
                          n_block = NULL,
                          values = NULL,
                          n_values = NULL,
                          force_log = FALSE) {
  if (is.null(n_values)) {
    if (is.null(values)) {
      n_values <- max(length(block), length(weights), length(log_weights))
    } else {
      n_values <- length(values)
    }
  }
  if (is.null(block) || (length(block) == 0)) {
    block <- rep(1L, n_values)
  } else if (length(block) == 1L) {
    block <- rep(block, n_values)
  }
  if (min(block) < 1L) {
    warning(paste0(
      "min(block) = ", min(block),
      " < 1L. Setting too small values to 1L."
    ))
    block <- pmax(1L, block)
  }
  if (is.null(n_block)) {
    n_block <- max(block)
  } else if (max(block) > n_block) {
    warning(paste0(
      max(block), " = max(block) > n_block = ",
      n_block,
      ". Setting too large values to n_block."
    ))
    block <- pmin(n_block, block)
  }

  if (force_log) {
    if (is.null(log_weights)) {
      if (is.null(weights)) {
        log_weights <- rep(0.0, n_values)
      } else {
        log_weights <- log(weights)
        weights <- NULL
      }
    } else if (!is.null(weights)) {
      warning("Both weights and log_weights supplied. Using log_weights.",
        immediate. = TRUE
      )
      weights <- NULL
    }
  } else {
    if (is.null(log_weights) && is.null(weights)) {
      weights <- rep(1.0, n_values)
    } else if (!is.null(weights)) {
      # log_weights is non-NULL
      if (!is.null(log_weights)) {
        warning("Both weights and log_weights supplied. Using log_weights.",
          immediate. = TRUE
        )
        weights <- NULL
      }
    }
  }
  if (length(weights) == 1L) {
    weights <- rep(weights, n_values)
  }
  if (length(log_weights) == 1L) {
    log_weights <- rep(log_weights, n_values)
  }

  list(
    block = block, weights = weights, log_weights = log_weights,
    n_block = n_block
  )
}
