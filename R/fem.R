#' @include deprecated.R

# fm_fem ####

#' @title Compute finite element matrices
#'
#' @description Compute finite element mass and structure matrices
#'
#' @param mesh [fm_mesh_1d()], [fm_mesh_2d()], or other supported mesh class
#' object
#' @param order integer; the maximum operator order
#' @param ... Currently unused
#'
#' @export
#' @examples
#' names(fm_fem(fm_mesh_1d(1:4), order = 3))
#' names(fm_fem(fmexample$mesh, order = 3))
#'
fm_fem <- function(mesh, order = 2, ...) {
  UseMethod("fm_fem")
}

make_symmetric <- function(x) {
  (x + Matrix::t(x)) / 2
}

#' @rdname fm_fem
#' @returns `fm_fem.fm_mesh_1d`: A list with elements `c0`, `c1`, `g1`, `g2`,
#' etc.
#' When `mesh$degree == 2`, also `g01`, `g02`, and `g12`.
#' @export
fm_fem.fm_mesh_1d <- function(mesh, order = 2, ...) {
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

    c0_ <- c0
    c0 <- Matrix::Diagonal(mesh$m, c0_)

    g_list <- list(
      g1 = g1,
      g2 = make_symmetric(
        (Matrix::t(g1) %*% Matrix::Diagonal(mesh$m, 1 / c0_)) %*% g1
      )
    )
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
        derivatives = TRUE,
        full = TRUE
      )
    c1 <- Matrix::t(info$A) %*% info$A
    g1 <- Matrix::t(info$dA) %*% info$dA
    g2 <- Matrix::t(info$d2A) %*% info$d2A

    g01 <- Matrix::t(info$A) %*% info$dA
    g02 <- Matrix::t(info$A) %*% info$d2A
    g12 <- Matrix::t(info$dA) %*% info$d2A

    c0_ <- Matrix::rowSums(c1)
    c0 <- Matrix::Diagonal(nrow(c1), c0_)

    g_list <- list(g1 = g1, g2 = g2, g01 = g01, g02 = g02, g12 = g12)
  } else {
    stop(paste("Mesh basis degree=", mesh$degree,
      " is not supported by fm_fem.fm_mesh_1d.",
      sep = ""
    ))
  }

  if (order > 2) {
    tmp <- Matrix::t(g_list[["g2"]]) %*% Matrix::Diagonal(mesh$m, 1 / c0_)
    for (k in seq_len(order - 2) + 2) {
      g_list[[paste0("g", k)]] <-
        make_symmetric(tmp %*% g_list[[paste0("g", k - 2)]])
    }
  }

  return(c(list(c0 = c0, c1 = c1), g_list))
}

#' @rdname fm_fem
#' @param aniso If non-NULL, a `list(gamma, v)`. Calculates anisotropic
#'   structure matrices (in addition to the regular) for \eqn{\gamma}{gamma} and
#'   \eqn{v}{v} for an anisotropic operator \eqn{\nabla\cdot H \nabla}{div H
#'   grad}, where \eqn{H=\gamma I + v v^\top}{H = gamma I + v v'}. Currently
#'   (2023-08-05) the fields need to be given per vertex.
#' @returns `fm_fem.fm_mesh_2d`: A list with elements `c0`, `c1`, `g1`, `va`,
#'   `ta`, and more if `order > 1`. When `aniso` is non-NULL, also `g1aniso`
#'   matrices, etc.
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
#' @returns `fm_fem.fm_tensor`: A list with elements `cc`, `g1`, `g2`.
#' @export
fm_fem.fm_tensor <- function(mesh, order = 2, ...) {
  if (order > 2) {
    warning("Only fem order <= 2 implemented for fm_tensor")
    order <- 2
  }

  fem_list <- lapply(mesh$fun_spaces, fm_fem, order = order)
  cc_list <- lapply(seq_along(mesh$fun_spaces), function(i) {
    if (inherits(mesh$fun_spaces[[i]], "fm_mesh_1d") &&
      mesh$fun_spaces[[i]]$degree == 2) {
      return(fem_list[[i]]$c1)
    }
    fem_list[[i]]$c0
  })

  kron_multi <- function(x) {
    if (length(x) == 1) {
      return(x[[1]])
    }
    result <- x[[1]]
    for (k in seq_len(length(x) - 1)) {
      result <- kronecker(x[[k + 1]], result)
    }
    result
  }

  cc <- kron_multi(cc_list)
  g1 <- cc * 0.0
  g2 <- cc * 0.0
  mat_list <- cc_list
  for (i in seq_along(fem_list)) {
    mat_list[[i]] <- fem_list[[i]]$g2
    g2 <- g2 + kron_multi(mat_list)
    mat_list[[i]] <- fem_list[[i]]$g1
    g1 <- g1 + kron_multi(mat_list)
    for (j in seq_len(i - 1L)) {
      mat_list[[j]] <- fem_list[[j]]$g1
      g2 <- g2 + 2.0 * kron_multi(mat_list)
      mat_list[[j]] <- cc_list[[j]]
    }
    mat_list[[i]] <- cc_list[[i]]
  }

  return(list(cc = cc, g1 = g1, g2 = g2))
}



#' @rdname fm_fem
#' @returns `fm_fem.fm_collect`: A list with elements `c0`, `c1`,
#' `g1`, `g2`, etc, and `cc` (`c0` for every model except `fm_mesh_1d` with
#' `degree=2`, for which it is `c1`). If the base type for the collection
#' provides `va` and `ta` values, those are also returned.
#' @export
fm_fem.fm_collect <- function(mesh, order = 2, ...) {
  fem_list <- lapply(mesh$fun_spaces, fm_fem, order = order)
  cc_list <- lapply(seq_along(mesh$fun_spaces), function(i) {
    if (inherits(mesh$fun_spaces[[i]], "fm_mesh_1d") &&
      mesh$fun_spaces[[i]]$degree == 2) {
      return(fem_list[[i]]$c1)
    }
    fem_list[[i]]$c0
  })

  result <- list(
    cc = Matrix::.bdiag(cc_list),
    c0 = Matrix::.bdiag(lapply(fem_list, function(x) x$c0)),
    c1 = Matrix::.bdiag(lapply(fem_list, function(x) x$c1))
  )
  if ("va" %in% names(fem_list[[1]])) {
    result$va <- unlist(lapply(fem_list, function(x) x$va))
  }
  if ("ta" %in% names(fem_list[[1]])) {
    result$ta <- unlist(lapply(fem_list, function(x) x$ta))
  }
  for (k in seq_len(order)) {
    result[[paste0("g", k)]] <-
      Matrix::.bdiag(lapply(fem_list, function(x) x[[paste0("g", k)]]))
  }

  result
}




row_cross_product <- function(e1, e2) {
  cbind(
    e1[, 2] * e2[, 3] - e1[, 3] * e2[, 2],
    e1[, 3] * e2[, 1] - e1[, 1] * e2[, 3],
    e1[, 1] * e2[, 2] - e1[, 2] * e2[, 1]
  )
}

row_volume_product <- function(e1, e2, e3) {
  rowSums(row_cross_product(e1, e2) * e3)
}

#' @rdname fm_fem
#' @returns `fm_fem.fm_mesh_3d`: A list with elements `c0`, `c1`, `g1`, `g2`,
#'   `va`, `ta`, and more if `order > 2`.
#'
#' @export
fm_fem.fm_mesh_3d <- function(mesh, order = 2, ...) {
  if (order > 2) {
    warning("Only fem order <= 2 implemented for fm_mesh_3d")
  }
  v1 <- mesh$loc[mesh$graph$tv[, 1], , drop = FALSE]
  v2 <- mesh$loc[mesh$graph$tv[, 2], , drop = FALSE]
  v3 <- mesh$loc[mesh$graph$tv[, 3], , drop = FALSE]
  v4 <- mesh$loc[mesh$graph$tv[, 4], , drop = FALSE]
  e1 <- v2 - v1
  e2 <- v3 - v2
  e3 <- v4 - v3
  e4 <- v1 - v4
  vols_t <- abs(row_volume_product(e1, e2, e3)) / 6

  c0 <- Matrix::sparseMatrix(
    i = as.vector(mesh$graph$tv),
    j = as.vector(mesh$graph$tv),
    x = rep(vols_t / 4, times = 4),
    dims = c(mesh$n, mesh$n)
  )
  vols_v <- Matrix::diag(c0)

  # Sign changes for b2 and b4 for consistent in/out vector orientation
  b1 <- row_cross_product(e2, e3)
  b2 <- -row_cross_product(e3, e4)
  b3 <- row_cross_product(e4, e1)
  b4 <- -row_cross_product(e1, e2)

  g_i <- g_j <- g_x <- numeric(nrow(mesh$graph$tv) * 16)
  for (tt in seq_len(nrow(mesh$graph$tv))) {
    GG <- rbind(
      b1[tt, , drop = FALSE],
      b2[tt, , drop = FALSE],
      b3[tt, , drop = FALSE],
      b4[tt, , drop = FALSE]
    )
    ii <- (tt - 1) * 16 + seq_len(16)
    g_i[ii] <- rep(mesh$graph$tv[tt, ], each = 4)
    g_j[ii] <- rep(mesh$graph$tv[tt, ], times = 4)
    g_x[ii] <- as.vector((GG %*% t(GG)) / vols_t[tt] / 36)
  }
  g1 <- Matrix::sparseMatrix(
    i = g_i,
    j = g_j,
    x = g_x,
    dims = c(mesh$n, mesh$n)
  )

  list(
    c0 = c0,
    g1 = g1,
    g2 = g1 %*% Matrix::Diagonal(mesh$n, 1 / vols_v) %*% g1,
    va = vols_v,
    ta = vols_t
  )
}


#' @title fm_sizes
#' @description `r lifecycle::badge("experimental")`
#'   Compute effective sizes of faces/cells and vertices in a mesh
#' @param ... Passed on to submethods
#' @returns A `list` with elements `face` and `vertex` for 2D meshes, or `cell`
#'   and `vertex` for 3D meshes. The elements are vectors of effective sizes of
#'   the faces/cells and vertices, respectively.
#' @export
#' @examples
#' str(fm_sizes(fmexample$mesh))
#'
fm_sizes <- function(...) {
  UseMethod("fm_sizes")
}

#' @rdname fm_sizes
#' @param mesh object of a supported mesh class
#' @export
fm_sizes.fm_mesh_2d <- function(mesh, ...) {
  if (fm_manifold(mesh, "S")) {
    warning("`fm_sizes()` does not handle spherical triangles.")
  }
  v1 <- mesh$loc[mesh$graph$tv[, 1], , drop = FALSE]
  v2 <- mesh$loc[mesh$graph$tv[, 2], , drop = FALSE]
  v3 <- mesh$loc[mesh$graph$tv[, 3], , drop = FALSE]
  e1 <- v2 - v1
  e2 <- v3 - v2
  e3 <- v3 - v1
  areas_t <- rowSums((row_cross_product(e1, e2) +
    row_cross_product(e2, e3) +
    row_cross_product(e3, e1))^2)^0.5 / 6

  c0 <- Matrix::sparseMatrix(
    i = as.vector(mesh$graph$tv),
    j = as.vector(mesh$graph$tv),
    x = rep(areas_t / 3, times = 3),
    dims = c(mesh$n, mesh$n)
  )
  areas_v <- Matrix::diag(c0)

  list(face = areas_t, vertex = areas_v)
}

#' @rdname fm_sizes
#' @export
fm_sizes.fm_mesh_3d <- function(mesh, ...) {
  v1 <- mesh$loc[mesh$graph$tv[, 1], , drop = FALSE]
  v2 <- mesh$loc[mesh$graph$tv[, 2], , drop = FALSE]
  v3 <- mesh$loc[mesh$graph$tv[, 3], , drop = FALSE]
  v4 <- mesh$loc[mesh$graph$tv[, 4], , drop = FALSE]
  e1 <- v2 - v1
  e2 <- v3 - v2
  e3 <- v4 - v3
  e4 <- v1 - v4
  vols_t <- abs(row_volume_product(e1, e2, e3)) / 6

  c0 <- Matrix::sparseMatrix(
    i = as.vector(mesh$graph$tv),
    j = as.vector(mesh$graph$tv),
    x = rep(vols_t / 4, times = 4),
    dims = c(mesh$n, mesh$n)
  )
  vols_v <- Matrix::diag(c0)

  list(cell = vols_t, vertex = vols_v)
}
