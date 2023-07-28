# Point/mesh connection methods ####

# fm_evaluate ####

#' @title Methods for projecting to/from an inla.mesh
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
#' if (fm_safe_inla()) {
#'   n <- 20
#'   loc <- matrix(runif(n * 2), n, 2)
#'   mesh <- fm_rcdt_2d_inla(loc, refine = list(max.edge = 0.05))
#'   proj <- fm_evaluator(mesh)
#'   field <- cos(mesh$loc[, 1] * 2 * pi * 3) * sin(mesh$loc[, 2] * 2 * pi * 7)
#'   image(proj$x, proj$y, fm_evaluate(proj, field))
#' }
#' \donttest{
#' if (fm_safe_inla() &&
#'   require("ggplot2") &&
#'   fm_safe_sp() &&
#'   require("ggpolypath") &&
#'   require("inlabru")) {
#'   ggplot() +
#'     gg(mesh, col = field)
#' }
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
fm_evaluate <- function(...) {
  UseMethod("fm_evaluate")
}

#' @export
#' @rdname fm_evaluate
fm_evaluate.inla.mesh <- function(mesh, field, ...) {
  fm_evaluate.fm_mesh_2d(fm_as_mesh_2d(mesh), field = field, ...)
}
#' @export
#' @rdname fm_evaluate
fm_evaluate.fm_mesh_2d <- function(mesh, field, ...) {
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
fm_evaluate.fm_tensor <- function(mesh, field, ...) {
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
fm_evaluate.inla.mesh.1d <- function(mesh, field, ...) {
  fm_evaluate.fm_mesh_1d(fm_as_mesh_1d(mesh), field = field, ...)
}
#' @export
#' @rdname fm_evaluate
fm_evaluate.fm_mesh_1d <- function(mesh, field, ...) {
  if (missing(field) || is.null(field)) {
    lifecycle::deprecate_stop(
      "0.0.1",
      "fm_evaluate(field = ' must not be missing or NULL.')",
      "fm_evaluator()"
    )
  }

  proj <- fm_evaluator(mesh, ...)
  fm_evaluate(proj, field)
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


# fm_evaluator ####

#' @describeIn fm_evaluate
#' Returns an `fm_evaluator` list object with evaluation information.
#' The `proj` element contains a mapping matrix `A` and a logical vector `ok`,
#' that indicates which locations were mappable to the input mesh.
#' For `fm_mesh_2d` and `inla.mesh`
#' input, `proj` also contains a matrix `bary` and vector `t`, with the
#' barycentric coordinates within the triangle each input location falls in.
#' @export
fm_evaluator <- function(...) {
  UseMethod("fm_evaluator")
}


#' @param weights Optional weight vector
#' @param derivatives logical; If true, also return matrices `dA` and `d2A`
#' for `fm_mesh_1d` objects, and `dx`, `dy`, `dz` for `fm_mesh_2d`.
#' @export
#' @rdname fm_evaluate
fm_evaluator_mesh_2d <- function(mesh,
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
    x = as.vector(b[ii, ]) * weights[rep(ii, 3)]
  ))

  mesh_deriv <- function(mesh, info, weights) {
    n.mesh <- mesh$n

    ii <- which(info$ok)
    n.ok <- sum(info$ok)
    tv <- mesh$graph$tv[info$t[ii, 1L], ]
    e1 <- mesh$loc[tv[, 3], ] - mesh$loc[tv[, 2], ]
    e2 <- mesh$loc[tv[, 1], ] - mesh$loc[tv[, 3], ]
    e3 <- mesh$loc[tv[, 2], ] - mesh$loc[tv[, 1], ]
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

  info
}


#' @param method character; with "default", uses the object definition of the
#' function space. Otherwise overrides the object definition.
#' @export
#' @rdname fm_evaluate
fm_evaluator_mesh_1d <- function(mesh,
                                 loc,
                                 weights = NULL,
                                 derivatives = NULL,
                                 method = c(
                                   "default",
                                   "nearest",
                                   "linear",
                                   "quadratic"
                                 ),
                                 ...) {
  method <- match.arg(method)
  ok <- (loc >= mesh$interval[1]) & (loc <= mesh$interval[2])
  if (identical(method, "default")) {
    ## Compute basis based on mesh$degree and mesh$boundary
    if (mesh$degree == 0) {
      info <- fm_basis(
        mesh,
        loc,
        method = "nearest",
        weights = weights,
        derivatives = ifelse(is.null(derivatives), FALSE, TRUE)
      )
      if (mesh$boundary[1] == "dirichlet") {
        info$A <- info$A[, -1, drop = FALSE]
        if (!is.null(derivatives) && derivatives) {
          info$dA <- info$dA[, -1, drop = FALSE]
        }
      }
      if (mesh$boundary[2] == "dirichlet") {
        info$A <- info$A[, -(mesh$m + 1), drop = FALSE]
        if (!is.null(derivatives) && derivatives) {
          info$dA <- info$dA[, -(mesh$m + 1), drop = FALSE]
        }
      }
      info[["ok"]] <- ok
      return(info)
    } else if (mesh$degree == 1) {
      info <- fm_basis(
        mesh,
        loc,
        method = "linear",
        weights = weights,
        derivatives =
          ifelse(is.null(derivatives),
            FALSE, TRUE
          )
      )
      if (mesh$boundary[1] == "dirichlet") {
        info$A <- info$A[, -1, drop = FALSE]
        if (!is.null(derivatives) && derivatives) {
          info$dA <- info$dA[, -1, drop = FALSE]
        }
      }
      if (mesh$boundary[2] == "dirichlet") {
        info$A <- info$A[, -(mesh$m + 1), drop = FALSE]
        if (!is.null(derivatives) && derivatives) {
          info$dA <- info$dA[, -(mesh$m + 1), drop = FALSE]
        }
      }
      info[["ok"]] <- ok
      return(info)
    } else if (mesh$degree == 2) {
      info <-
        fm_basis(
          mesh,
          loc,
          method = "quadratic",
          weights = weights,
          derivatives =
            ifelse(is.null(derivatives),
              FALSE, TRUE
            )
        )
      adjust.matrix <- function(mesh, A) {
        A <- fm_as_dgTMatrix(A)
        i <- A@i + 1L
        j <- A@j + 1L
        x <- A@x
        if (mesh$cyclic) {
          j[j == 1L] <- mesh$m + 1L
          j[j == mesh$m + 2L] <- 2L
          j <- j - 1L
        } else {
          if (mesh$boundary[1] == "neumann") {
            j[j == 1L] <- 2L
            j <- j - 1L
          } else if (mesh$boundary[1] == "dirichlet") {
            x[j == 1L] <- -x[j == 1L]
            j[j == 1L] <- 2L
            j <- j - 1L
          } else if ((mesh$boundary[1] == "free") &&
            mesh$free.clamped[1]) {
            j1 <- which(j == 1L)
            i <- c(i, i[j1])
            j <- c(j, j[j1] + 1L)
            x <- c(x, -x[j1])
            x[j1] <- 2 * x[j1]
          }
          if (mesh$boundary[2] == "neumann") {
            j[j == (mesh$m + 1L)] <- mesh$m
          } else if (mesh$boundary[2] == "dirichlet") {
            x[j == (mesh$m + 1L)] <- -x[j == (mesh$m + 1L)]
            j[j == (mesh$m + 1L)] <- mesh$m
          } else if ((mesh$boundary[2] == "free") &&
            mesh$free.clamped[2]) {
            j1 <- which(j == (mesh$m))
            i <- c(i, i[j1])
            j <- c(j, j[j1] - 1L)
            x <- c(x, -x[j1])
            x[j1] <- 2 * x[j1]
          }
        }
        return(Matrix::sparseMatrix(
          i = i, j = j, x = x,
          dims = c(length(loc), mesh$m)
        ))
      }

      if (!is.null(derivatives) && derivatives) {
        return(list(
          A = adjust.matrix(mesh, info$A),
          dA = adjust.matrix(mesh, info$dA),
          d2A = adjust.matrix(mesh, info$d2A),
          ok = ok
        ))
      } else {
        info$A <- adjust.matrix(mesh, info$A)
        info$ok <- ok
        return(info)
      }
    } else {
      stop(paste("'degree' must be 0, 1, or 2.  'degree=",
        mesh$degree,
        "' is not supported.",
        sep = ""
      ))
    }
  } else {
    if (is.null(weights)) {
      weights <- rep(1, length(loc))
    }

    if (!is.na(pmatch(method, c("linear", "nearest")))) {
      idx <- fm_bary(mesh, loc, method = method)
      dA <- NULL

      if (method == "linear") {
        ## Compute the n 1st order B-splines.
        i <- rep(seq_along(loc), times = 2)
        weights.i <- weights[i]
        A <- (Matrix::sparseMatrix(
          i = i,
          j = as.vector(idx$index),
          x = weights.i * as.vector(idx$bary),
          dims = c(length(loc), mesh$n)
        ))
        if (!is.null(derivatives) && derivatives) {
          if (mesh$cyclic) {
            d <- (c(mesh$loc[-1], mesh$interval[2]) -
              mesh$loc)
          } else {
            d <- (mesh$loc[-1] - mesh$loc[-mesh$n])
          }
          dA <- (Matrix::sparseMatrix(
            i = i,
            j = as.vector(idx$index),
            x = (weights.i * c(
              -1 / d[idx$index[, 1]],
              1 / d[idx$index[, 1]]
            )),
            dims = c(length(loc), mesh$n)
          ))
        }
      } else {
        ## Nearest neighbours.
        A <- (Matrix::sparseMatrix(
          i = seq_along(loc),
          j = idx$index[, 1],
          x = weights * idx$bary[, 1],
          dims = c(length(loc), mesh$n)
        ))
      }

      if (!is.null(derivatives) && derivatives) {
        return(list(A = A, dA = dA, ok = ok))
      } else {
        return(list(A = A, ok = ok))
      }
    } else { ## Quadratic
      ## Compute the n+1 2nd order B-splines.

      if (mesh$cyclic) {
        Boundary.knots <-
          (c(mesh$loc, mesh$interval[2])[c(mesh$n, 2)] +
            diff(mesh$interval) * c(-1, 1))
        knots <- c(mesh$loc, mesh$interval[2])
      } else {
        Boundary.knots <-
          c(
            mesh$loc[1] - (mesh$loc[2] - mesh$loc[1]),
            mesh$loc[mesh$n] + (mesh$loc[mesh$n] - mesh$loc[mesh$n - 1])
          )
        knots <- mesh$loc
      }

      if (FALSE) { ## Only for debugging.
        ## Using bs():
        ## Note: Intermediate step constructs dense matrix.
        bsobj <-
          splines::bs(
            x = loc,
            knots = knots,
            degree = 2,
            Boundary.knots = Boundary.knots
          )
        bsobj <-
          Matrix::Matrix(
            as.vector(bsobj[, 1:(mesh$n + mesh$cyclic + 1)]),
            nrow(bsobj), mesh$n + mesh$cyclic + 1
          )
      }

      ## Direct calculation:
      knots <- c(Boundary.knots[1], knots)
      if (mesh$cyclic) {
        idx <-
          fm_bary(
            fm_mesh_1d(knots, boundary = "free"),
            (loc - mesh$interval[1]) %%
              diff(mesh$interval) + mesh$interval[1],
            method = "linear"
          )
      } else {
        idx <-
          fm_bary(fm_mesh_1d(knots, boundary = "free"),
            loc,
            method = "linear"
          )
      }
      knots <- c(knots, Boundary.knots[2])
      idx$index <- idx$index[, 1] - 1L ## Indices into mesh intervals.

      d <- knots[2:length(knots)] - knots[1:(length(knots) - 1)]
      d2 <- knots[3:length(knots)] - knots[1:(length(knots) - 2)]

      ## Left intervals for each basis function:
      i.l <- seq_along(idx$index)
      j.l <- idx$index + 2L
      x.l <- (idx$bary[, 2] * d[idx$index + 1] / d2[idx$index + 1] * idx$bary[, 2])
      ## Right intervals for each basis function:
      i.r <- seq_along(idx$index)
      j.r <- idx$index
      x.r <- (idx$bary[, 1] * d[idx$index + 1] / d2[idx$index] * idx$bary[, 1])
      ## Middle intervals for each basis function:
      i.m <- seq_along(idx$index)
      j.m <- idx$index + 1L
      x.m <- (1 - (idx$bary[, 1] * d[idx$index + 1] / d2[idx$index] * idx$bary[, 1] +
        idx$bary[, 2] * d[idx$index + 1] / d2[idx$index + 1] * idx$bary[, 2]
      ))

      i <- c(i.l, i.r, i.m)
      j <- c(j.l, j.r, j.m)
      weights.i <- weights[i]

      A <- (Matrix::sparseMatrix(
        i = i, j = j,
        x = weights.i * c(x.l, x.r, x.m),
        dims = c(length(idx$index), mesh$n + mesh$cyclic + 1L)
      ))

      dA <- NULL
      d2A <- NULL
      if (!is.null(derivatives) && derivatives) {
        ## dA:
        ## Left, right, middle intervals for each basis function:
        x.l <- (2 / d2[idx$index + 1] * idx$bary[, 2])
        x.r <- (-2 / d2[idx$index] * idx$bary[, 1])
        x.m <- (-(-2 / d2[idx$index] * idx$bary[, 1] +
          2 / d2[idx$index + 1] * idx$bary[, 2]
        ))
        dA <- (Matrix::sparseMatrix(
          i = i, j = j,
          x = weights.i * c(x.l, x.r, x.m),
          dims = (c(
            length(idx$index),
            mesh$n + mesh$cyclic + 1L
          ))
        ))

        ## d2A:
        ## Left, right, middle intervals for each basis function:
        x.l <- (2 / d[idx$index + 1] / d2[idx$index + 1])
        x.r <- (2 / d[idx$index + 1] / d2[idx$index])
        x.m <- (-(2 / d[idx$index + 1] / d2[idx$index] +
          2 / d[idx$index + 1] / d2[idx$index + 1]))
        d2A <- (Matrix::sparseMatrix(
          i = i, j = j,
          x = weights.i * c(x.l, x.r, x.m),
          dims = (c(
            length(idx$index),
            mesh$n + mesh$cyclic + 1L
          ))
        ))
      }

      if (is.null(derivatives) && !derivatives) {
        return(list(A = A, ok = ok))
      } else {
        return(list(A = A, dA = dA, d2A = d2A, ok = ok))
      }
    }
  }
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
#' @describeIn fm_evaluate The `...` arguments are passed on to `fm_evaluator_lattice()`
#' if no `loc` or `lattice` is provided.
fm_evaluator.inla.mesh <- function(mesh,
                                   loc = NULL,
                                   lattice = NULL,
                                   crs = NULL,
                                   ...) {
  fm_evaluator.fm_mesh_2d(
    fm_as_mesh_2d(mesh),
    loc = loc,
    lattice = lattice,
    crs = crs,
    ...
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
      proj <- fm_evaluator_mesh_2d(mesh, lattice$loc)
    } else {
      proj <- fm_evaluator_mesh_2d(mesh,
        loc = lattice$loc,
        crs = lattice$crs
      )
    }
    projector <- list(x = x, y = y, lattice = lattice, loc = NULL, proj = proj, crs = crs)
    class(projector) <- "fm_evaluator"
  } else {
    proj <- fm_evaluator_mesh_2d(mesh, loc = loc, crs = crs)
    projector <- list(x = NULL, y = NULL, lattice = NULL, loc = loc, proj = proj, crs = crs)
    class(projector) <- "fm_evaluator"
  }

  return(projector)
}


#' @export
#' @rdname fm_evaluate
fm_evaluator.inla.mesh.1d <- function(mesh,
                                      loc = NULL,
                                      xlim = mesh$interval,
                                      dims = 100,
                                      ...) {
  fm_evaluator.fm_mesh_1d(fm_as_mesh_1d(mesh), loc = loc, xlim = xlim, dims = dims, ...)
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

  proj <- fm_evaluator_mesh_1d(mesh, loc)
  projector <- list(x = loc, lattice = NULL, loc = loc, proj = proj)
  class(projector) <- "fm_evaluator"

  return(projector)
}


#' @param x [fm_tensor()] object
#' @export
#' @rdname fm_evaluate
fm_evaluator.fm_tensor <- function(x,
                                   loc,
                                   ...) {
  if (length(loc) != length(x[["fun_spaces"]])) {
    stop("")
  }
  if (is.null(names(loc))) {
    names(loc) <- names(x[["fun_spaces"]])
  } else if (!setequal(names(x[["fun_spaces"]]), names(loc))) {
    stop("")
  }
  proj <- lapply(
    names(x[["fun_spaces"]]),
    function(k) {
      fm_evaluator(x[["fun_spaces"]][[k]], loc = loc[[k]])$proj
    }
  )

  # Combine the matrices
  # (A1, A2, A3) -> rowkron(A3, rowkron(A2, A1))
  A <- proj[[1]][["A"]]
  ok <- proj[[1]][["ok"]]
  for (k in seq_len(length(x[["fun_spaces"]]) - 1)) {
    A <- fm_row_kron(proj[[k + 1]][["A"]], A)
    ok <- proj[[k + 1]][["ok"]] & ok
  }

  structure(
    list(proj = list(A = A, ok = ok)),
    class = "fm_evaluator"
  )
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
#' if (fm_safe_inla() &&
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
#' \dontrun{
#' if (fm_safe_inla(quietly = TRUE)) {
#'   # Load Gorilla data
#'
#'   data("gorillas", package = "inlabru")
#'
#'   # Check if all Gorilla nests are inside the mesh
#'
#'   all(fm_is_within(gorillas$nests, gorillas$mesh))
#'
#'   # Also works for locations not stored as SpatialPoints object
#'
#'   loc <- coordinates(gorillas$nests)
#'   all(fm_is_within(loc, gorillas$mesh))
#' }
#' }
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
#' @param \dots Currently unused
#' @returns A `sparseMatrix`
#' @examples
#' \dontrun{
#' if (fm_safe_inla(quietly = TRUE)) {
#'   # Load Gorilla data
#'
#'   data("gorillas", package = "inlabru")
#'
#'   # Compute basis mapping matrix
#'   str(fm_basis(gorillas$mesh, gorillas$nests))
#' }
#' }
#' @export
fm_basis <- function(x, ...) {
  UseMethod("fm_basis")
}

#' @rdname fm_basis
#' @export
fm_basis.default <- function(x, loc, ...) {
  fm_evaluator(x, loc = loc)$proj$A
}

#' @param weights Optional weight matrix to apply (from the left)
#' @param derivatives If non-NULL and logical, return a list, optionally
#' including derivative matrices.
#' @returns For `fm_mesh_1d`, a list with elements
#' \item{A }{The projection matrix, `u(loc_i)=sum_j A_ij w_i`}
#' \item{d1A, d2A }{Derivative weight matrices,
#' `du/dx(loc_i)=sum_j dx_ij w_i`, etc.}
#' @rdname fm_basis
#' @export
fm_basis.fm_mesh_1d <- function(x, loc, weights = NULL, derivatives = NULL, ...) {
  result <- fm_evaluator_mesh_1d(
    x,
    loc = loc,
    weights = weights,
    derivatives = derivatives,
    ...
  )
  if (is.null(derivatives)) {
    return(result$A)
  }
  return(result)
}

#' @return For `fm_mesh_2d`, a list with elements
#' \item{A }{The projection matrix, `u(loc_i)=sum_j A_ij w_i`}
#' \item{dx, dy, dz }{Derivative weight matrices, `du/dx(loc_i)=sum_j
#' dx_ij w_i`, etc.}
#' @rdname fm_basis
#' @export
fm_basis.fm_mesh_2d <- function(x, loc, weights = NULL, derivatives = NULL, ...) {
  result <- fm_evaluator_mesh_2d(
    x,
    loc = loc,
    weights = weights,
    derivatives = derivatives,
    ...
  )
  if (is.null(derivatives)) {
    return(result$A)
  }
  return(result)
}

#' @rdname fm_basis
#' @export
#' @method fm_basis inla.mesh.1d
fm_basis.inla.mesh.1d <- function(x, loc, ...) {
  fm_basis.fm_mesh_1d(fm_as_mesh_1d(x), loc = loc, ...)
}

#' @rdname fm_basis
#' @export
#' @method fm_basis inla.mesh
fm_basis.inla.mesh <- function(x, loc, ...) {
  fm_basis.fm_mesh_2d(fm_as_mesh_2d(x), loc = loc, ...)
}

#' @rdname fm_basis
#' @export
fm_basis.fm_evaluator <- function(x, ...) {
  x$proj$A
}











# Block methods ####

block_prep <- function(block = NULL, log_weights = NULL, weights = NULL) {
  block_len <- NULL
  if (is.null(block)) {
    if (!is.null(log_weights)) {
      block_len <- length(log_weights)
    } else {
      block_len <- length(weights)
    }
    if (block_len == 0) {
      return(integer(0))
    }
    block <- rep(1L, block_len)
  } else if (length(block) == 1L) {
    if (!is.null(log_weights)) {
      block_len <- length(log_weights)
    } else {
      block_len <- length(weights)
    }
    if (block_len > 1L) {
      block <- rep(block, block_len)
    }
  }
  block
}

#' Blockwise aggregation matrices
#'
#' Creates an aggregation matrix for blockwise aggregation, with optional
#' weighting.
#'
#' @param block integer vector; block information. If `NULL`,
#' `rep(1L, block_len)` is used, where `block_len` is determined by
#' `length(log_weights)))` or `length(weights)))`.A single scalar is also repeated
#' to a vector of corresponding length to the weights.
#' @param weights Optional weight vector
#' @param log_weights Optional `log(weights)` vector. Overrides `weights` when
#' non-NULL.
#' @param rescale logical; If `TRUE`, normalise the weights by `sum(weights)`
#' or `sum(exp(log_weights))` within each block.
#' Default: `FALSE`
#' @param n_block integer; The number of conceptual blocks. Only needs to be
#' specified if it's larger than `max(block)`.
#'
#' @export
#'
#' @describeIn fm_block A (sparse) matrix of size `n_block` times `length(block)`.
#' @examples
#' fm_block(c(1, 1, 1, 2, 2))
#' fm_block(c(1, 1, 1, 2, 2), rescale = TRUE)
#' fm_block(c(1, 1, 1, 2, 2), log_weights = -2:2, rescale = TRUE)
fm_block <- function(block = NULL,
                     weights = NULL,
                     log_weights = NULL,
                     rescale = FALSE,
                     n_block = NULL) {
  block <- block_prep(block = block, log_weights = log_weights, weights = weights)
  if (length(block) == 0) {
    return(
      Matrix::sparseMatrix(
        i = integer(0),
        j = integer(0),
        x = numeric(0),
        dims = c(0L, 0L)
      )
    )
  }
  if (is.null(n_block)) {
    n_block <- max(block)
  }
  weights <-
    fm_block_weights(
      block = block,
      weights = weights,
      log_weights = log_weights,
      n_block = n_block,
      rescale = rescale
    )

  Matrix::sparseMatrix(
    i = block,
    j = seq_len(length(block)),
    x = weights,
    dims = c(n_block, length(block))
  )
}


#' @describeIn fm_block Computes (optionally) blockwise renormalised weights
#' @export
fm_block_weights <- function(block = NULL, weights = NULL, log_weights = NULL,
                             n_block = NULL, rescale = FALSE) {
  block <- block_prep(block = block, log_weights = log_weights, weights = weights)
  if (length(block) == 0) {
    return(numeric(0))
  }
  if (is.null(n_block)) {
    n_block <- max(block)
  }
  if (rescale) {
    # Compute blockwise normalised weights
    if (!is.null(log_weights)) {
      log_weights <- fm_block_log_weights(
        log_weights = log_weights,
        block = block,
        n_block = n_block,
        rescale = rescale
      )
      weights <- exp(log_weights)
    } else {
      if (is.null(weights)) {
        weights <- rep(1.0, length(block))
      }
      scale <- as.vector(
        Matrix::sparseMatrix(
          i = block,
          j = rep(1L, length(block)),
          x = weights,
          dims = c(n_block, 1)
        )
      )
      weights <- weights / scale[block]
    }
  } else if (!is.null(log_weights)) {
    weights <- exp(log_weights)
  } else if (is.null(weights)) {
    weights <- rep(1.0, length(block))
  }
  weights
}

#' @describeIn fm_block Computes (optionally) blockwise renormalised log-weights
#' @export
fm_block_log_weights <- function(block = NULL, weights = NULL, log_weights = NULL,
                                 n_block = NULL,
                                 rescale = FALSE) {
  block <- block_prep(block = block, log_weights = log_weights, weights = weights)
  if (length(block) == 0) {
    return(numeric(0))
  }
  if (is.null(n_block)) {
    n_block <- max(block)
  }
  if (is.null(log_weights)) {
    if (is.null(weights)) {
      log_weights <- rep(0.0, length(block))
    } else {
      log_weights <- log(weights)
    }
  }
  if (rescale) {
    shift <- fm_block_log_shift(
      block = block, log_weights = log_weights,
      n_block = n_block
    )
    log_rescale <- as.vector(
      Matrix::sparseMatrix(
        i = block,
        j = rep(1L, length(block)),
        x = exp(log_weights - shift[block]),
        dims = c(n_block, 1)
      )
    )
    log_rescale <- (log(log_rescale) + shift)[block]
    log_weights <- log_weights - log_rescale
  }
  log_weights
}


#' @describeIn fm_block Computes shifts for stable blocked log-sum-exp.
#' To compute \eqn{\log(\sum_{i; \text{block}_i=k} \exp(v_i) w_i)}{
#' log(sum_(i;block_i=k) exp(v_i) w_i)
#' } for
#' each block `k`:
#' ```
#' w_values <- values + fm_block_log_weights(block, log_weights)
#' shift <- fm_block_log_shift(block, w_values)
#' agg <- aggregate(exp(w_values - shift[block]),
#'                  by = list(block = block),
#'                  \(x) log(sum(x)))
#' agg$x <- agg$x + shift[agg$block]
#' ```
#' @export
fm_block_log_shift <- function(block = NULL, log_weights = NULL, n_block = NULL) {
  block <- block_prep(block = block, log_weights = log_weights, weights = NULL)
  if (length(block) == 0) {
    return(0.0)
  }
  if (is.null(n_block)) {
    n_block <- max(block)
  }
  block_k <- sort(unique(block))
  shift <- numeric(n_block)
  if (!is.null(log_weights)) {
    shift[block_k] <-
      vapply(
        block_k,
        function(k) {
          max(log_weights[block == k])
        },
        0.0
      )
  }
  shift
}
