# fm_is_within ####

#' @title Query if points are inside a mesh
#'
#' @description
#'  Queries whether each input point is within a mesh or not.
#'
#' @param x A set of points/locations of a class supported by `fm_basis(y, loc =
#'   x, ..., full = TRUE)`
#' @param y An [fm_mesh_2d] or other class supported by
#' `fm_basis(y, loc = x, ..., full = TRUE)`
#' @param \dots Passed on to [fm_basis()]
#' @returns A logical vector
#' @examples
#' all(fm_is_within(fmexample$loc, fmexample$mesh))
#' @export
fm_is_within <- function(x, y, ...) {
  fm_basis(y, loc = x, ..., full = TRUE)$ok
}


# fm_basis ####

#' @title Compute mapping matrix between mesh function space and points
#'
#' @description Computes the basis mapping matrix between a function space on a
#' mesh, and locations.
#'
#' @param x An function space object, or other supported object
#'   (`matrix`, `Matrix`, `list`)
#' @param loc A location/value information object (`numeric`, `matrix`, `sf`,
#' `fm_bary`, etc, depending on the class of `x`)
#' @param full logical; if `TRUE`, return a `fm_basis` object, containing at
#'   least a projection matrix `A` and logical vector `ok` indicating which
#'   evaluations are valid. If `FALSE`, return only the projection matrix `A`.
#'   Default is `FALSE`.
#' @param \dots Passed on to submethods
#' @returns A `sparseMatrix` object (if `full = FALSE`), or a `fm_basis` object
#'   (if `full = TRUE` or `isTRUE(derivatives)`). The `fm_basis` object contains
#'   at least the projection matrix `A` and logical vector `ok`; If `x_j`
#'   denotes the latent basis coefficient for basis function `j`, the field is
#'   defined as `u(loc_i)=sum_j A_ij x_j` for all `i` where `ok[i]` is `TRUE`,
#'   and `u(loc_i)=0.0` where `ok[i]` is `FALSE`.
#' @seealso [fm_raw_basis()]
#' @examples
#' # Compute basis mapping matrix
#' dim(fm_basis(fmexample$mesh, fmexample$loc))
#' print(fm_basis(fmexample$mesh, fmexample$loc, full = TRUE))
#'
#' # From precomputed `fm_bary` information:
#' bary <- fm_bary(fmexample$mesh, fmexample$loc)
#' print(fm_basis(fmexample$mesh, bary, full = TRUE))
#' @export
fm_basis <- function(x, ..., full = FALSE) {
  UseMethod("fm_basis")
}

#' @rdname fm_basis
#' @export
fm_basis.default <- function(x, ..., full = FALSE) {
  lifecycle::deprecate_stop(
    "0.1.7.9002",
    "fm_basis.default()",
    details = "Each mesh class needs its own `fm_basis()` method."
  )
}

#' @param derivatives If non-NULL and logical, include derivative matrices
#' in the output. Forces `full = TRUE`.
#' @describeIn fm_basis If `derivatives=TRUE`, the `fm_basis` object contains
#'   additional derivative weight matrices, `d1A` and `d2A`, `du/dx(loc_i)=sum_j
#'   dx_ij w_i`.
#' @export
fm_basis.fm_mesh_1d <- function(x,
                                loc,
                                weights = NULL,
                                derivatives = NULL,
                                ...,
                                full = FALSE) {
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

#' @describeIn fm_basis If `derivatives=TRUE`, additional derivative weight
#'   matrices are included in the `full=TRUE` output: Derivative weight matrices
#' `dx`, `dy`, `dz`; `du/dx(loc_i)=sum_j dx_ij w_i`, etc.
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

#' @describeIn fm_basis `fm_mesh_3d` basis functions.
#' @export
fm_basis.fm_mesh_3d <- function(x, loc, weights = NULL, ...,
                                full = FALSE) {
  bary <- fm_bary(x, loc, ...)
  n_loc <- NROW(bary)
  ok <- !is.na(bary$index)
  simplex <- fm_bary_simplex(x, bary[ok, , drop = FALSE])
  if (is.null(weights)) {
    weights <- rep(1.0, n_loc)
  } else if (length(weights) == 1) {
    weights <- rep(weights, n_loc)
  }
  A <- sparseMatrix_nonzero(
    i = rep(which(ok), 4),
    j = as.vector(simplex),
    x = as.numeric(as.vector(bary$where[ok, ]) * weights[rep(which(ok), 4)]),
    dims = c(n_loc, fm_dof(x))
  )

  fm_basis(list(A = A, ok = ok, bary = bary), full = full)
}

#' @describeIn fm_basis `fm_lattice_2d` bilinear basis functions.
#' @export
fm_basis.fm_lattice_2d <- function(x, loc, weights = NULL, ...,
                                   full = FALSE) {
  bary <- fm_bary(x, loc, ...)
  n_loc <- NROW(bary)
  ok <- !is.na(bary$index)
  simplex <- fm_bary_simplex(x, bary[ok, , drop = FALSE])
  if (is.null(weights)) {
    weights <- rep(1.0, n_loc)
  } else if (length(weights) == 1) {
    weights <- rep(weights, n_loc)
  }
  A <- sparseMatrix_nonzero(
    i = rep(which(ok), 4),
    j = as.vector(simplex),
    x = as.numeric(as.vector(bary$where[ok, ]) * weights[rep(which(ok), 4)]),
    dims = c(n_loc, fm_dof(x))
  )

  fm_basis(list(A = A, ok = ok), full = full)
}

#' @describeIn fm_basis `fm_lattice_Nd` multilinear basis functions.
#' @export
fm_basis.fm_lattice_Nd <- function(x, loc, weights = NULL, ...,
                                   full = FALSE) {
  bary <- fm_bary(x, loc, ...)
  n_loc <- NROW(bary)
  ok <- !is.na(bary$index)
  simplex <- fm_bary_simplex(x, bary[ok, , drop = FALSE])
  if (is.null(weights)) {
    weights <- rep(1.0, n_loc)
  } else if (length(weights) == 1) {
    weights <- rep(weights, n_loc)
  }
  A <- sparseMatrix_nonzero(
    i = rep(which(ok), ncol(bary$where)),
    j = as.vector(simplex),
    x = as.numeric(as.vector(bary$where[ok, ]) *
      weights[rep(which(ok), ncol(bary$where))]),
    dims = c(n_loc, fm_dof(x))
  )

  fm_basis(list(A = A, ok = ok), full = full)
}

#' @export
#' @describeIn fm_basis Evaluates a basis matrix for a `fm_tensor` function
#'   space.
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

  fm_basis(
    list(A = A, ok = ok),
    full = full
  )
}


#' @export
#' @describeIn fm_basis Evaluates a basis matrix for a `fm_collect` function
#'   space. The `loc` argument must be a `list` or `tibble` with elements
#'   `loc` (the locations) and `index` (the indices into the function space
#'   collection).
#' @importFrom rlang .env
fm_basis.fm_collect <- function(x,
                                loc,
                                weights = NULL,
                                ...,
                                full = FALSE) {
  loc_names <- names(loc)
  if (!is.null(loc_names) &&
    (!("loc" %in% loc_names) || !("index" %in% loc_names))) {
    stop(
      paste0(
        "Location data for fm_collect must have elements `loc` and ",
        "`index`.\n",
        "Found: ", paste0(names(loc), collapse = ", ")
      )
    )
  }

  if (!tibble::is_tibble(loc)) {
    if (is.null(loc_names)) {
      # The .env construction is needed to avoid `loc` name clash effects
      loc <- tibble::tibble(
        loc = .env$loc[[1]],
        index = .env$loc[[2]]
      )
    } else {
      loc <- tibble::tibble(
        loc = .env$loc[["loc"]],
        index = .env$loc[["index"]]
      )
    }
  }

  if (is.numeric(loc[["index"]]) && !is.integer(loc[["index"]])) {
    loc[["index"]] <- as.integer(loc[["index"]])
  }
  if (!is.null(names(x[["fun_spaces"]]))) {
    if (is.factor(loc[["index"]])) {
      loc[["index"]] <- as.character(loc[["index"]])
    }
    if (is.character(loc[["index"]])) {
      # Convert character indices to integer
      loc[["index"]] <- match(loc[["index"]], names(x[["fun_spaces"]]))
    }
  }

  idx <- seq_along(x[["fun_spaces"]])
  valid <- loc[["index"]] %in% idx

  proj <- lapply(
    idx,
    function(k) {
      fm_basis(
        x[["fun_spaces"]][[k]],
        loc = loc[loc[["index"]] == k, , drop = FALSE][["loc"]],
        full = TRUE
      )
    }
  )

  # Combine the matrices
  A <- Matrix::.bdiag(lapply(proj, fm_basis))
  ok <- do.call(c, lapply(proj, function(xx) xx[["ok"]]))

  # Reorder to original order and fill in invalid rows
  block_order <- order(loc[["index"]][valid])
  reorder <- order(block_order)
  A_ <- A[reorder, , drop = FALSE]
  ok_ <- ok[reorder]

  A <- Matrix::sparseMatrix(
    i = integer(0),
    j = integer(0),
    x = numeric(0),
    dims = c(nrow(loc), ncol(A_))
  )
  A[valid, ] <- A_
  ok <- logical(nrow(loc))
  ok[valid] <- ok_

  if (!is.null(weights)) {
    A <- Matrix::Diagonal(n = nrow(A), x = weights) %*% A
  }

  fm_basis(
    list(A = A, ok = ok),
    full = full
  )
}


#' @describeIn fm_basis Creates a new `fm_basis` object with elements `A` and
#'   `ok`, from a pre-evaluated basis matrix, including optional additional
#'   elements in the `...` arguments. If a `ok` is `NULL`, it is inferred as
#'   `rep(TRUE, NROW(x))`, indicating that all rows correspond to successful
#'   basis evaluations. If `full = FALSE`,
#'   returns the matrix unchanged.
#' @param ok numerical of length `NROW(x)`, indicating which rows of `x` are
#'   valid/successful basis evaluations. If `NULL`, inferred as
#'   `rep(TRUE, NROW(x))`.
#' @param weights Optional weight vector to apply (from the left, one
#' weight for each row of the basis matrix)
#' @export
fm_basis.matrix <- function(x, ok = NULL, weights = NULL, ..., full = FALSE) {
  if (!full && is.null(weights)) {
    return(x)
  }
  fm_basis(list(A = x, ok = ok, ...), weights = weights, full = full)
}

#' @describeIn fm_basis Creates a new `fm_basis` object with elements `A` and
#'   `ok`, from a pre-evaluated basis matrix, including optional additional
#'   elements in the `...` arguments. If a `ok` is `NULL`, it is inferred as
#'   `rep(TRUE, NROW(x))`, indicating that all rows correspond to successful
#'   basis evaluations. If `full = FALSE`,
#'   returns the matrix unchanged.
#' @export
fm_basis.Matrix <- function(x, ok = NULL, weights = NULL, ..., full = FALSE) {
  if (!full && is.null(weights)) {
    return(x)
  }
  fm_basis(list(A = x, ok = ok, ...), weights = weights, full = full)
}

#' @describeIn fm_basis Creates a new `fm_basis` object from a plain list
#'   containing at least an element `A`. If an `ok` element is missing,
#'   it is inferred as `rep(TRUE, NROW(x$A))`. If `full = FALSE`,
#'   extracts the `A` matrix.
#' @export
fm_basis.list <- function(x, weights = NULL, ..., full = FALSE) {
  stopifnot("A" %in% names(x))
  if (!is.null(weights)) {
    x[["A"]] <- Matrix::Diagonal(nrow(x[["A"]]), x = weights) %*% x[["A"]]
  }
  if (!full) {
    return(x[["A"]])
  }
  if (is.null(x[["ok"]])) {
    x[["ok"]] <- rep(TRUE, NROW(x[["A"]]))
  } else if (!is.logical(x[["ok"]]) ||
    (length(x[["ok"]]) != NROW(x[["A"]]))) {
    stop(
      "Invalid 'ok' element in 'x'; should be a logical vector of length ",
      NROW(x[["A"]])
    )
  }
  structure(x, class = "fm_basis")
}

#' @describeIn fm_basis If `full` is `TRUE`, returns `x` unchanged, otherwise
#'   returns the `A` matrix contained in `x`.
#' @export
fm_basis.fm_basis <- function(x, ..., full = FALSE) {
  if (full) {
    x
  } else {
    x[["A"]]
  }
}

#' @describeIn fm_basis Extract `fm_basis` information from an `fm_evaluator`
#'   object. If `full = FALSE`, returns the `A` matrix contained in the
#'   `fm_basis` object.
#' @export
fm_basis.fm_evaluator <- function(x, ..., full = FALSE) {
  fm_basis(x$proj, full = full)
}


internal_spline_mesh_1d <- function(interval,
                                    m,
                                    degree,
                                    boundary,
                                    free.clamped) {
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
  fm_mesh_1d(seq(interval[1], interval[2], length.out = n),
    degree = degree,
    boundary = boundary,
    free.clamped = free.clamped
  )
}


# fmesher_spherical_harmonics_gsl <- function(loc,
#                                             max_order,
#                                             rot_inv) {
#   n <- max_order
#   loc <- loc / rowSums(loc^2)^0.5
#   if (rot_inv) {
#     basis <- matrix(0, nrow(loc), n + 1)
#     for (l in seq(0, n)) {
#       basis[, l + 1] <- sqrt(2 * l + 1) *
#         gsl::legendre_Pl(l = l, x = loc[, 3])
#     }
#   } else {
#     angle <- atan2(loc[, 2], loc[, 1])
#     basis <- matrix(0, nrow(loc), (n + 1)^2)
#     for (l in seq(0, n)) {
#       basis[, 1 + l * (l + 1)] <-
#         sqrt(2 * l + 1) *
#         gsl::legendre_Pl(l = l, x = loc[, 3])
#       for (m in seq_len(l)) {
#         scaling <- sqrt(2 * (2 * l + 1) * exp(lgamma(l - m + 1) -
#                                                 lgamma(l + m + 1)))
#         poly <- gsl::legendre_Plm(l = l, m = m, x = loc[, 3])
#         basis[, 1 + l * (l + 1) - m] <-
#           scaling * sin(-m * angle) * poly
#         basis[, 1 + l * (l + 1) + m] <-
#           scaling * cos(m * angle) * poly
#       }
#     }
#   }
#   basis
# }


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
#' @author Finn Lindgren <Finn.Lindgren@@gmail.com>
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
    # Make sure we have radius-1 coordinates
    loc <- mesh$loc / rowSums(mesh$loc^2)^0.5
    basis <- fmesher_spherical_harmonics(
      loc,
      max_order = as.integer(n),
      rot_inv = isTRUE(rot.inv)
    )
  }

  basis
}

# Create sparse matrix with no explicit zeros
sparseMatrix_nonzero <- function(i, j, x, dims) {
  nonzero <- (x != 0)
  i <- i[nonzero]
  j <- j[nonzero]
  x <- x[nonzero]
  Matrix::sparseMatrix(i = i, j = j, x = x, dims = dims)
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
  if (!inherits(loc, "fm_bary")) {
    loc <- fm_bary(mesh, loc = loc, crs = crs, ...)
  }
  n_loc <- NROW(loc)

  ok <- !is.na(loc$index)

  if (is.null(weights)) {
    weights <- rep(1.0, n_loc)
  } else if (length(weights) == 1) {
    weights <- rep(weights, n_loc)
  }

  ii <- which(ok)
  A <- (sparseMatrix_nonzero(
    dims = c(n_loc, mesh$n),
    i = rep(ii, 3),
    j = as.vector(mesh$graph$tv[loc$index[ii], ]),
    x = as.numeric(as.vector(loc$where[ii, ]) * weights[rep(ii, 3)])
  ))

  mesh_deriv <- function(mesh, bary, ok, weights) {
    n.mesh <- mesh$n

    ii <- which(ok)
    n.ok <- sum(ok)
    tv <- mesh$graph$tv[bary$index[ii], , drop = FALSE]
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
    dx <- (sparseMatrix_nonzero(
      dims = c(n_loc, n.mesh),
      i = rep(ii, 3),
      j = as.vector(tv),
      x = as.vector(x) * weights[rep(ii, 3)]
    ))
    dy <- (sparseMatrix_nonzero(
      dims = c(n_loc, n.mesh),
      i = rep(ii, 3),
      j = as.vector(tv),
      x = as.vector(y) * weights[rep(ii, 3)]
    ))
    dz <- (sparseMatrix_nonzero(
      dims = c(n_loc, n.mesh),
      i = rep(ii, 3),
      j = as.vector(tv),
      x = as.vector(z) * weights[rep(ii, 3)]
    ))

    list(dx = dx, dy = dy, dz = dz)
  }

  info <- list(bary = loc, A = A, ok = ok)

  if (!is.null(derivatives) && derivatives) {
    info <-
      c(
        info,
        mesh_deriv(
          mesh = mesh,
          bary = info$bary,
          ok = info$ok,
          weights = weights
        )
      )
  }

  fm_basis(info, full = TRUE)
}


#' @export
#' @rdname fm_basis_helpers
fm_basis_mesh_1d <- function(mesh,
                             loc,
                             weights = NULL,
                             derivatives = NULL,
                             ...) {
  if (is.null(weights)) {
    weights <- rep(1.0, NROW(loc))
  } else if (length(weights) == 1L) {
    weights <- rep(weights, NROW(loc))
  }

  derivatives <- !is.null(derivatives) && derivatives
  info_ <- list()

  ## Compute basis based on mesh$degree and mesh$boundary
  if (mesh$degree == 0) {
    info <- fm_bary(mesh, loc = loc, method = "nearest")
    info_ <- list(bary = info)
    bary_ok <- !is.na(info$index)
    i_ <- seq_along(loc)[bary_ok]
    j_ <- info$index[bary_ok]
    x_ <- info$where[bary_ok, 1]
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
    info_ <- list(bary = info)
    bary_ok <- !is.na(info$index)
    info <- info[bary_ok, ]
    i_ <- c(which(bary_ok), which(bary_ok))
    simplex <- fm_bary_simplex(mesh, info)
    j_ <- as.vector(simplex)
    x_ <- as.vector(info$where)
    if (derivatives) {
      j_curr <- simplex[, 1]
      j_next <- simplex[, 2]
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
      x_d <- rep(c(-1, 1), each = sum(bary_ok)) / rep(dist, times = 2)
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
      inter <- c(0, diff(mesh$interval))
    } else {
      knots <- mesh$loc - mesh$loc[1]
      inter <- range(knots)
    }
    if (!inherits(loc, "fm_bary")) {
      loc <- loc - mesh$loc[1]
    }

    # Note: If loc is `fm_bary`, it's also valid for this local fm_mesh_1d.
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
    info_ <- list(bary = info)
    bary_ok <- !is.na(info$index)
    info <- info[bary_ok, ]

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
      simplex <- fm_bary_simplex(mesh, info)
      i.l <- which(bary_ok)
      j.l <- simplex[, 1] + 2L
      x.l <- (info$where[, 2] * d[simplex[, 2]] / d2[simplex[, 2]] *
        info$where[, 2])
      if (derivatives) {
        x.d1.l <- (2 / d2[simplex[, 2]] * info$where[, 2])
        x.d2.l <- (2 / d2[simplex[, 2]] / d[simplex[, 2]])
      }
      ## Right intervals for each basis function:
      i.r <- seq_along(simplex[, 1])
      j.r <- simplex[, 1]
      x.r <- (info$where[, 1] * d[simplex[, 2]] / d2[simplex[, 1]] *
        info$where[, 1])
      if (derivatives) {
        x.d1.r <- -(2 / d2[simplex[, 2]] * info$where[, 1])
        x.d2.r <- (2 / d2[simplex[, 1]] / d[simplex[, 2]])
      }
      ## Middle intervals for each basis function:
      i.m <- seq_along(simplex[, 1])
      j.m <- simplex[, 1] + 1L
      x.m <- (1 - (info$where[, 1] * d[simplex[, 2]] / d2[simplex[, 1]] *
        info$where[, 1] +
        info$where[, 2] * d[simplex[, 2]] / d2[simplex[, 2]] *
          info$where[, 2]))
      if (derivatives) {
        x.d1.m <- (2 / d2[simplex[, 1]] * info$where[, 1]) -
          (2 / d2[simplex[, 2]] * info$where[, 2])
        x.d2.m <- -(2 / d2[simplex[, 1]] / d[simplex[, 2]]) -
          (2 / d2[simplex[, 2]] / d[simplex[, 2]])
      }
    } else {
      d2 <- c(2 * d[1], 2 * d[1], d2, 2 * d[length(d)], 2 * d[length(d)])
      d <- c(d[1], d[2], d, d[length(d)], d[length(d)])
      simplex <- fm_bary_simplex(mesh, info)
      ok <- (simplex[, 1] >= 1L) & (simplex[, 2] <= length(knots))
      index <- simplex[ok, , drop = FALSE] + 1L
      bary <- info$where[ok, , drop = FALSE]
      ## Left intervals for each basis function:
      i.l <- which(bary_ok)[ok]
      j.l <- index[, 2]
      x.l <- (bary[, 2] * d[index[, 2]] / d2[index[, 2]] * bary[, 2])
      if (derivatives) {
        x.d1.l <- (2 / d2[index[, 2]] * bary[, 2])
        x.d2.l <- (2 / d2[index[, 2]] / d[index[, 2]])
      }
      ## Right intervals for each basis function:
      i.r <- which(bary_ok)[ok]
      j.r <- index[, 1] - 1L
      x.r <- (bary[, 1] * d[index[, 2]] / d2[index[, 1]] * bary[, 1])
      if (derivatives) {
        x.d1.r <- -(2 / d2[index[, 1]] * bary[, 1])
        x.d2.r <- (2 / d2[index[, 1]] / d[index[, 2]])
      }
      ## Middle intervals for each basis function:
      i.m <- which(bary_ok)[ok]
      j.m <- index[, 1]
      x.m <- (1 - (bary[, 1] * d[index[, 2]] / d2[index[, 1]] * bary[, 1] +
        bary[, 2] * d[index[, 2]] / d2[index[, 2]] * bary[, 2]
      ))
      if (derivatives) {
        x.d1.m <- (2 / d2[index[, 1]] * bary[, 1]) -
          (2 / d2[index[, 2]] * bary[, 2])
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
      simplex <- fm_bary_simplex(mesh, info)

      # Convert boundary basis functions to linear
      # First remove anything from above outside the interval, then add back in
      # the appropriate values
      if (inherits(loc, "fm_bary")) {
        ok <- (loc$where[bary_ok, 1] >= 0) & (loc$where[bary_ok, 2] >= 0)
      } else {
        ok <- (loc[bary_ok] >= inter[1]) & (loc[bary_ok] <= inter[2])
      }
      i_ <- i_[ok]
      j_ <- j_[ok]
      x_ <- x_[ok]
      if (derivatives) {
        x_d1 <- x_d1[ok]
        x_d2 <- x_d2[ok]
      }

      # left
      if (inherits(loc, "fm_bary")) {
        ok <- (loc$where[bary_ok, 2] < 0) & (simplex[, 1] == 1L)
      } else {
        ok <- (loc[bary_ok] < 0) & (simplex[, 1] == 1L)
      }
      i_l <- c(which(bary_ok)[ok], which(bary_ok)[ok])
      j_l <- c(simplex[ok, 1], simplex[ok, 2])
      x_l <- c(
        0.5 + (info$where[ok, 1] - 1),
        0.5 - (info$where[ok, 1] - 1)
      )
      if (derivatives) {
        x_d1_l <- rep(c(-1, 1) / d[1], each = sum(ok))
        x_d2_l <- rep(c(0, 0), each = sum(ok))
      }

      # right
      if (inherits(loc, "fm_bary")) {
        ok <- (loc$where[bary_ok, 1] < 0) & (simplex[, 2] == length(knots))
      } else {
        ok <- (loc[bary_ok] > inter[2]) & (simplex[, 2] == length(knots))
      }
      i_r <- c(which(bary_ok)[ok], which(bary_ok)[ok])
      j_r <- c(simplex[ok, 2], simplex[ok, 1]) + 1L
      x_r <- c(
        0.5 + (info$where[ok, 2] - 1),
        0.5 - (info$where[ok, 2] - 1)
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

  info_$A <- sparseMatrix_nonzero(
    i_,
    j_,
    weights[i_] * x_,
    dims = c(NROW(loc), mesh$m)
  )
  if (derivatives) {
    if (mesh$degree <= 1) {
      info_$dA <- sparseMatrix_nonzero(
        i = i_d, j = j_d, x = weights[i_d] * x_d,
        dims = c(NROW(loc), mesh$m)
      )
    } else {
      # degree is 2
      info_$dA <- sparseMatrix_nonzero(
        i = i_,
        j = j_,
        x = weights[i_] * x_d1,
        dims = c(NROW(loc), mesh$m)
      )
      info_$d2A <- sparseMatrix_nonzero(
        i = i_,
        j = j_,
        x = weights[i_] * x_d2,
        dims = c(NROW(loc), mesh$m)
      )
    }
  }

  info_[["ok"]] <- bary_ok

  fm_basis(info_, full = TRUE)
}


# Plain B-spline basis evaluation by Farin eq 10.13-10.14,
# building the basis function matrices recursively via index vectors
internal_bspline <- function(x, knots, degree = 1, deriv = 0) {
  if (min(x) < min(knots)) {
    stop("Some x out of range (too small)")
  }
  if (max(x) > max(knots)) {
    stop("Some x out of range (too large)")
  }
  k <- findInterval(x, knots, all.inside = TRUE)
  basis <- list(
    i = seq_along(x),
    j = k,
    values = numeric(length(x)) + 1.0
  )
  #  message("knots: ", knots)
  #  message("unique j: ", unique(basis$j))
  if (degree == 0) {
    return(basis)
  }
  knots <- c(rep(min(knots), degree - 1), knots, rep(max(knots), degree - 1))
  basis$j <- basis$j + degree
  #  message("knots: ", knots)
  #  message("unique j: ", unique(basis$j))
  l_range <- unique(sort(basis$j)) - 1L
  for (deg in seq_len(degree)) {
    basis_prev <- basis
    basis <- list(
      i = integer(0),
      j = integer(0),
      x = numeric(0),
      values = numeric(0)
    )
    l_range <- unique(sort(c(l_range - 1L, l_range)))
    for (l in l_range) {
      # +1L since j is base-1 but l is base-0
      left <- basis_prev$j == l + 1L
      right <- basis_prev$j == l + 2L
      if (any(left)) {
        basis$i <- c(basis$i, basis_prev$i[left])
        basis$j <- c(basis$j, basis_prev$j[left])
        basis$values <- c(
          basis$values,
          basis_prev$values[left] *
            (x[basis_prev$i[left]] - knots[l]) /
            (knots[l + deg] - knots[l])
        )
      }
      if (any(right)) {
        basis$i <- c(basis$i, basis_prev$i[right])
        basis$j <- c(basis$j, basis_prev$j[right] - 1L)
        basis$values <- c(
          basis$values,
          basis_prev$values[right] *
            (knots[l + deg + 1L] - x[basis_prev$i[right]]) /
            (knots[l + deg + 1L] - knots[l + 1L])
        )
      }
    }
    #    message("knots: ", knots)
    #    message("unique j: ", unique(basis$j))
  }
  basis
}


# Plain B-spline basis evaluation by Farin eq 10.13-10.14,
# building the basis function matrices recursively via index vectors
internal_bspline2 <- function(x, knots, degree = 1, deriv = 0) {
  if (min(x) < min(knots)) {
    stop("Some x out of range (too small)")
  }
  if (max(x) > max(knots)) {
    stop("Some x out of range (too large)")
  }

  if (deriv > 0) {
    basis_lower <-
      internal_bspline2(x, knots, degree = degree - 1, deriv = deriv - 1)
    m <- length(knots) + degree - 1L
    m_lower <- m - 1L
    A <- sparseMatrix_nonzero(
      i = basis_lower$i,
      j = basis_lower$j,
      x = basis_lower$values,
      dims = c(length(x), m_lower)
    )
    basis_diff <- sparseMatrix_nonzero(
      i = c(seq_len(m_lower), seq_len(m_lower)),
      j = c(seq_len(m_lower), seq_len(m_lower) - 1L),
      x = rep(c(1, -1), c(length(knots) - 1, length(knots) - 1)),
      dims = c(m_lower, m)
    )
    A <- A %*% basis_diff
    return(A)
  }

  k <- findInterval(x, knots, all.inside = TRUE)
  basis <- list(
    i = seq_along(x),
    j = k,
    values = numeric(length(x)) + 1.0
  )
  if (degree == 0) {
    return(basis)
  }
  knots <- c(rep(min(knots), degree - 1), knots, rep(max(knots), degree - 1))
  basis$j <- basis$j + degree
  l_range <- unique(sort(basis$j)) - 1L
  for (deg in seq_len(degree)) {
    basis_prev <- basis
    l_range <- unique(sort(c(l_range - 1L, l_range)))
    # +1L since j is base-1 but l is base-0
    left <- basis_prev$j %in% (l_range + 1L)
    right <- basis_prev$j %in% (l_range + 2L)
    sz <- c(sum(left), sum(right))
    basis <- list(
      i = integer(sum(sz)),
      j = integer(sum(sz)),
      values = numeric(sum(sz))
    )
    if (sz[1] > 0L) {
      idx_left <- seq_len(sz[1])
      i <- basis_prev$i[left]
      j <- basis_prev$j[left]
      l <- j - 1L
      basis$i[idx_left] <- i
      basis$j[idx_left] <- j
      basis$values[idx_left] <- basis_prev$values[left] *
        (x[i] - knots[l]) /
        (knots[l + deg] - knots[l])
    }
    if (sz[2] > 0L) {
      idx_right <- seq_len(sz[2]) + sz[1]
      i <- basis_prev$i[right]
      j <- basis_prev$j[right]
      l <- j - 2L
      basis$i[idx_right] <- i
      basis$j[idx_right] <- j - 1L
      basis$values[idx_right] <- basis_prev$values[right] *
        (knots[l + deg + 1L] - x[i]) /
        (knots[l + deg + 1L] - knots[l + 1L])
    }
  }
  basis
}


# Block methods ####

#' Blockwise aggregation matrices
#'
#' Creates an aggregation matrix for blockwise aggregation, with optional
#' weighting.
#'
#' @param block integer vector; block information. If `NULL`,
#'   `rep(1L, block_len)` is used, where `block_len` is determined by
#'   `length(log_weights)))` or `length(weights)))`. A single scalar is also
#'   repeated to a vector of corresponding length to the weights.
#'
#'   Note: from version `0.2.0.9017` to `0.4.0.9005`, 'character'
#'   input was converted to integer with `as.integer(factor(block))`. As this
#'   could lead to unintended ordering of the output, this is no longer allowed.
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
#' @describeIn fm_block A (sparse) matrix of size `n_block` times
#'   `length(block)`.
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

  sparseMatrix_nonzero(
    i = info$block,
    j = seq_along(info$block),
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

  if (FALSE) {
    val <-
      Matrix::sparseMatrix(
        i = info$block,
        j = rep(1L, length(info$block)),
        x = as.numeric(values * weights),
        dims = c(info$n_block, 1)
      )
    as.vector(val)
  } else {
    agg <- stats::aggregate(
      data.frame(x = values * weights),
      by = list(block = info[["block"]]),
      FUN = sum,
      simplify = TRUE,
      drop = TRUE
    )
    val <- numeric(info$n_block)
    val[agg[["block"]]] <- agg[["x"]]
    val
  }
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
  if (is.character(block)) {
    lifecycle::deprecate_stop(
      "0.4.0.9006",
      "fm_block_prep(block = 'as `character` is no longer supported')",
      details =
        c(
          "Converting character block information to integer",
          "with `as.integer(factor(block))` is no longer supported,",
          "as it may lead to incorrect ordering of the results."
        )
    )
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
