#' @include deprecated.R

# fm_fem ####

#' @title Compute finite element matrices
#'
#' @description (...)
#'
#' @param mesh `inla.mesh` object
#' @param order integer
#' @param ... Currently unused
#'
#' @export
fm_fem <- function(mesh, order = 2, ...) {
  UseMethod("fm_fem")
}

#' @rdname fm_fem
#' @return `fm_fem.fm_mesh_1d`: A list with elements `c0`, `c1`, `g1`, `g2`.
#' When `mesh$degree == 2`, also `g01`, `g02`, and `g12`.
#' @export
fm_fem.fm_mesh_1d <- function(mesh, order = 2, ...) {
  if (order > 2) {
    warning("Only fem order <= 2 implemented for fm_mesh_1d")
    order <- 2
  }

  ## Use the same matrices for degree 0 as for degree 1
  if ((mesh$degree == 0) || (mesh$degree == 1)) {
    if (mesh$cyclic) {
      loc <-
        c(
          mesh$loc[mesh$n] - diff(mesh$interval),
          mesh$loc,
          mesh$loc[1] + diff(mesh$interval)
        )
      c0 <- (loc[3:length(loc)] - loc[1:(length(loc) - 2)]) / 2
      c1.l <- (loc[2:(length(loc) - 1)] - loc[1:(length(loc) - 2)]) / 6
      c1.r <- (loc[3:length(loc)] - loc[2:(length(loc) - 1)]) / 6
      c1.0 <- (c1.l + c1.r) * 2
      g1.l <- -1 / (loc[2:(length(loc) - 1)] - loc[1:(length(loc) - 2)])
      g1.r <- -1 / (loc[3:length(loc)] - loc[2:(length(loc) - 1)])
      g1.0 <- -g1.l - g1.r
      i.l <- seq_len(mesh$n)
      i.r <- seq_len(mesh$n)
      i.0 <- seq_len(mesh$n)
      if (mesh$n > 1) {
        j.l <- c(mesh$n, 1:(mesh$n - 1))
        j.r <- c(2:mesh$n, 1)
        j.0 <- 1:mesh$n
      } else {
        j.l <- 1L
        j.r <- 1L
        j.0 <- 1L
      }
    } else {
      c0 <-
        c(
          (mesh$loc[2] - mesh$loc[1]) / 2,
          (mesh$loc[mesh$n] - mesh$loc[mesh$n - 1]) / 2
        )
      if (mesh$n > 2) {
        c0 <-
          c(
            c0[1],
            (mesh$loc[3:mesh$n] - mesh$loc[1:(mesh$n - 2)]) / 2,
            c0[2]
          )
      }
      c1.l <- (mesh$loc[2:mesh$n] - mesh$loc[1:(mesh$n - 1)]) / 6
      c1.r <- c1.l
      c1.0 <- (c(0, c1.l) + c(c1.r, 0)) * 2
      g1.l <- -1 / (mesh$loc[2:mesh$n] - mesh$loc[1:(mesh$n - 1)])
      g1.r <- g1.l
      g1.0 <- -c(0, g1.l) - c(g1.r, 0)
      i.l <- 2:mesh$n
      i.r <- 1:(mesh$n - 1)
      i.0 <- 1:mesh$n
      j.l <- 1:(mesh$n - 1)
      j.r <- 2:mesh$n
      j.0 <- 1:mesh$n

      if (mesh$boundary[1] == "dirichlet") {
        g1.0 <- g1.0[-1]
        g1.l <- g1.l[-1]
        g1.r <- g1.r[-1]
        c1.0 <- c1.0[-1]
        c1.l <- c1.l[-1]
        c1.r <- c1.r[-1]
        c0 <- c0[-1]
        i.l <- i.l[-1] - 1
        i.r <- i.r[-1] - 1
        i.0 <- i.0[-1] - 1
        j.l <- j.l[-1] - 1
        j.r <- j.r[-1] - 1
        j.0 <- j.0[-1] - 1
      } else if (mesh$boundary[1] == "free") {
        g1.0[1] <- 0
        g1.r[1] <- 0
      }
      if (mesh$boundary[2] == "dirichlet") {
        m <- mesh$m
        g1.0 <- g1.0[-(m + 1)]
        g1.l <- g1.l[-m]
        g1.r <- g1.r[-m]
        c1.0 <- c1.0[-(m + 1)]
        c1.l <- c1.l[-m]
        c1.r <- c1.r[-m]
        c0 <- c0[-(m + 1)]
        i.l <- i.l[-m]
        i.r <- i.r[-m]
        i.0 <- i.0[-(m + 1)]
        j.l <- j.l[-m]
        j.r <- j.r[-m]
        j.0 <- j.0[-(m + 1)]
      } else if (mesh$boundary[2] == "free") {
        g1.0[mesh$m] <- 0
        g1.l[mesh$m - 1] <- 0
      }
    }

    g1 <-
      Matrix::sparseMatrix(
        i = c(i.l, i.r, i.0),
        j = c(j.l, j.r, j.0),
        x = c(g1.l, g1.r, g1.0),
        dims = c(mesh$m, mesh$m)
      )
    c1 <-
      Matrix::sparseMatrix(
        i = c(i.l, i.r, i.0),
        j = c(j.l, j.r, j.0),
        x = c(c1.l, c1.r, c1.0),
        dims = c(mesh$m, mesh$m)
      )
    g2 <- Matrix::t(g1) %*% Matrix::Diagonal(mesh$m, 1 / c0) %*% g1
    c0 <- Matrix::Diagonal(mesh$m, c0)
  } else if (mesh$degree == 2) {
    if (mesh$cyclic) {
      knots1 <- mesh$loc
      knots2 <- c(mesh$loc[-1], mesh$interval[2])
    } else {
      knots1 <- mesh$loc[-mesh$n]
      knots2 <- mesh$loc[-1]
    }
    knots.m <- (knots1 + knots2) / 2
    knots.d <- (knots2 - knots1) / 2
    ## 3-point Gaussian quadrature
    info <-
      fm_basis(mesh,
        loc = (c(
          knots.m,
          knots.m - knots.d * sqrt(3 / 5),
          knots.m + knots.d * sqrt(3 / 5)
        )),
        weights =
          c(knots.d * 8 / 9, knots.d * 5 / 9, knots.d * 5 / 9)^0.5,
        derivatives = TRUE
      )
    c1 <- Matrix::t(info$A) %*% info$A
    g1 <- Matrix::t(info$dA) %*% info$dA
    g2 <- Matrix::t(info$d2A) %*% info$d2A

    g01 <- Matrix::t(info$A) %*% info$dA
    g02 <- Matrix::t(info$A) %*% info$d2A
    g12 <- Matrix::t(info$dA) %*% info$d2A

    c0 <- Matrix::Diagonal(nrow(c1), Matrix::rowSums(c1))

    return(list(c0 = c0, c1 = c1, g1 = g1, g2 = g2, g01 = g01, g02 = g02, g12 = g12))
  } else {
    stop(paste("Mesh basis degree=", mesh$degree,
      " is not supported by fm_fem.fm_mesh_1d.",
      sep = ""
    ))
  }

  return(list(c0 = c0, c1 = c1, g1 = g1, g2 = g2))
}

#' @rdname fm_fem
#' @param aniso If non-NULL, a `list(gamma, v)`. Calculates anisotropic structure
#' matrices (in addition to the regular) for \eqn{\gamma}{gamma} and \eqn{v}{v} for
#' an anisotropic operator \eqn{\nabla\cdot H \nabla}{div H grad}, where
#' \eqn{H=\gamma I + v v^\top}{H = gamma I + v v'}.
#' Currently (2023-08-05) the fields need to be given per vertex.
#' @return `fm_fem.fm_mesh_2d`: A list with elements `c0`, `c1`, `g1`, `va`, `ta`,
#' and more if `order > 1`. When `aniso` is non-NULL, also `g1aniso` matrices, etc.
#'
#' @export
fm_fem.fm_mesh_2d <- function(mesh, order = 2,
                              aniso = NULL,
                              ...) {
  if (length(order) != 1) {
    stop("'order' must have length 1.")
  }
  if (!is.null(aniso)) {
    if (!is.list(aniso) || length(aniso) != 2) {
      stop("'aniso' must be NULL or a list of length 2.")
    }
  }
  result <- fmesher_fem(
    mesh_loc = mesh$loc,
    mesh_tv = mesh$graph$tv - 1L,
    fem_order_max = order,
    aniso = aniso,
    options = list()
  )
  result
}

#' @rdname fm_fem
#' @export
#' @method fm_fem inla.mesh.1d
fm_fem.inla.mesh.1d <- function(mesh, order = 2, ...) {
  fm_fem(fm_as_fm(mesh), order = order, ...)
}

#' @rdname fm_fem
#' @export
#' @method fm_fem inla.mesh
fm_fem.inla.mesh <- function(mesh, order = 2, ...) {
  fm_fem(fm_as_fm(mesh), order = order, ...)
}
