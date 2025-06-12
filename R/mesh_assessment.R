## meshassessment.R
##
## Copyright (C) 2016, 2018, 2025, Finn Lindgren

#' Interactive mesh building and diagnostics
#'
#' Assess the finite element approximation errors in a mesh for interactive R
#' sessions.
#'
#' @param mesh An [fm_mesh_2d] object
#' @param spatial.range numeric; the spatial range parameter to use for the
#'   assessment
#' @param alpha numeric; A valid [fm_matern_precision()] `alpha` parameter
#' @param dims 2-numeric; the grid size
#' @returns An `sf` object with gridded mesh assessment information
#' @author Finn Lindgren <Finn.Lindgren@@gmail.com>
#' @seealso [fm_mesh_2d()], [fm_rcdt_2d]
#' @examples
#'
#' bnd <- fm_segm(cbind(
#'   c(0, 10, 10, 0, 0),
#'   c(0, 0, 10, 10, 0)
#' ), is.bnd = TRUE)
#' mesh <- fm_rcdt_2d_inla(boundary = bnd, max.edge = 1)
#' out <- fm_assess(mesh, spatial.range = 3, alpha = 2)
#'
#' @export
fm_assess <- function(mesh, spatial.range, alpha = 2,
                      dims = NULL) {
  mesh.edgelengths <- function(mesh, proj) {
    i <- c()
    val <- c()
    num <- c()
    tri_idx <- rbind(c(1, 2), c(2, 3), c(3, 1))
    for (k in 1:3) {
      ti1 <- tri_idx[k, 1]
      ti2 <- tri_idx[k, 2]
      len <- Matrix::rowSums((mesh$loc[mesh$graph$tv[, ti2], ] -
        mesh$loc[mesh$graph$tv[, ti1], ])^2)^0.5
      i <- c(i, mesh$graph$tv[, ti1], mesh$graph$tv[, ti2])
      val <- c(val, rep(as.vector(len), times = 2))
      num <- c(num, rep(1, 2 * length(as.vector(len))))
    }
    avg_len <-
      as.vector(Matrix::sparseMatrix(i = i, j = rep(1, length(i)), x = val)) /
        as.vector(Matrix::sparseMatrix(i = i, j = rep(1, length(i)), x = num))

    b <- fm_basis(proj, full = TRUE)
    proj_len <- as.vector(b$A %*% avg_len)
    proj_len[!b$ok] <- NA
    proj_len
  }
  mesh.proj <- function(mesh, dims) {
    fm_evaluator(mesh, dims = dims)
  }
  mesh.spde <- function(mesh, alpha) {
    list(
      mesh = mesh,
      alpha = alpha,
      prior.range = c(1, 0.5),
      prior.sigma = c(1, 0.5)
    )
  }
  mesh.Q <- function(spde, spatial.range) {
    fm_matern_precision(spde$mesh,
      alpha = spde$alpha,
      rho = spatial.range,
      sigma = 1
    )
  }
  mesh.S <- function(Q) {
    fm_qinv(Q)
  }
  mesh.sd <- function(proj, S) {
    b <- fm_basis(proj, full = TRUE)
    v <- Matrix::rowSums(b$A * (b$A %*% S))
    v[!b$ok] <- NA
    array(v^0.5, dim = proj$lattice$dims)
  }
  mesh.sd.deviation.approx <- function(proj, S, sd0) {
    b <- fm_basis(proj, full = TRUE)
    val <- b$A %*% (
      as.vector(Matrix::t(b$A[b$ok, , drop = FALSE]) %*%
        as.vector(sd0)[b$ok]) /
        Matrix::colSums(b$A[b$ok, , drop = FALSE]))
    val[!b$ok] <- NA
    array(
      1 + (as.vector(sd0) - val),
      dim = proj$lattice$dims
    )
  }
  mesh.sd.bound <- function(proj, S) {
    fm_evaluate(proj, field = Matrix::diag(S)^0.5)
  }

  if (!fm_manifold(mesh, "R2")) {
    warning("fm_assess has only been tested on flat 2D manifolds.")
  }

  d <- fm_manifold_dim(mesh)
  if (is.null(dims)) {
    dims <- rep(ceiling(1e5^(1 / d)), d)
  }

  spde <- mesh.spde(mesh = mesh, alpha = alpha)
  Q <- mesh.Q(spde = spde, spatial.range = spatial.range)
  S <- mesh.S(Q = Q)
  proj <- mesh.proj(mesh = mesh, dims = dims)
  sd0 <- mesh.sd(proj = proj, S = S)
  sd.deviation <- mesh.sd.deviation.approx(proj = proj, S = S, sd0 = sd0)
  edgelengths <- mesh.edgelengths(mesh, proj)
  sd.bound <- mesh.sd.bound(proj, S)

  if (inherits(mesh, "fm_mesh_2d")) {
    out <- tibble::tibble(
      x = proj$lattice$loc[, 1],
      y = proj$lattice$loc[, 2],
      sd = as.vector(sd0),
      sd.dev = as.vector(sd.deviation),
      sd.bound = as.vector(sd.bound),
      edge.len = edgelengths
    )
    out <- sf::st_as_sf(out, coords = c("x", "y"), crs = fm_crs(mesh))
  } else {
    out <- tibble::tibble(
      loc = proj$lattice$loc,
      sd = as.vector(sd0),
      sd.dev = as.vector(sd.deviation),
      sd.bound = as.vector(sd.bound),
      edge.len = edgelengths
    )
  }
  out
}
