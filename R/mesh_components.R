#' @title Compute connected mesh subsets
#'
#' @description Compute subsets of vertices and triangles/tetrahedrons in an
#'   [fm_mesh_2d] or [fm_mesh_3d] object that are connected by edges/triangles.
#'
#' @return A list with elements `vertex` and `triangle`/`tetra`, vectors of
#' integer labels for which connected component they belong, and `info`, a
#' `data.frame` with columns
#' \item{component}{Connected component integer label.}
#' \item{nV}{The number of vertices in the component.}
#' \item{nT}{The number of triangles/tetrahedrons in the component.}
#' \item{area/volume}{The surface area or volume associated with the component.
#' Component labels are not comparable across
#' different meshes, but some ordering stability is guaranteed by initiating
#' each component from the lowest numbered triangle whenever a new component is
#' initiated.}
#'
#' @param mesh An [fm_mesh_2d] or [fm_mesh_3d] object
#' @author Finn Lindgren <Finn.Lindgren@@gmail.com>
#' @seealso [fm_mesh_2d()], [fm_rcdt_2d()], [fm_mesh_3d()]
#' @export
fm_mesh_components <- function(mesh) {
  UseMethod("fm_mesh_components")
}

#' @rdname fm_mesh_components
#' @export
#' @examples
#'
#' # Construct two simple meshes:
#' loc <- matrix(c(0, 1, 0, 1), 2, 2)
#' mesh1 <- fm_mesh_2d(loc = loc, max.edge = 0.1)
#' bnd <- fm_nonconvex_hull(loc, 0.3)
#' mesh2 <- fm_mesh_2d(boundary = bnd, max.edge = 0.1)
#'
#' # Compute connectivity information:
#' conn1 <- fm_mesh_components(mesh1)
#' conn2 <- fm_mesh_components(mesh2)
#' # One component, simply connected mesh
#' conn1$info
#' # Two disconnected components
#' conn2$info
#'
#' # Extract the subset mesh for each component:
#' # (Note: some information is lost, such as fixed segments,
#' # and boundary edge labels.)
#' mesh3_1 <- fm_rcdt_2d_inla(
#'   loc = mesh2$loc,
#'   tv = mesh2$graph$tv[conn2$triangle == 1, , drop = FALSE],
#'   delaunay = FALSE
#' )
#' mesh3_2 <- fm_rcdt_2d_inla(
#'   loc = mesh2$loc,
#'   tv = mesh2$graph$tv[conn2$triangle == 2, , drop = FALSE],
#'   delaunay = FALSE
#' )
#'
#' if (require("ggplot2")) {
#'   ggplot() +
#'     geom_fm(data = mesh3_1, fill = "red", alpha = 0.5) +
#'     geom_fm(data = mesh3_2, fill = "blue", alpha = 0.5)
#' }
fm_mesh_components.fm_mesh_2d <- function(mesh) {
  vertex <- integer(mesh$n)
  Nt <- nrow(mesh$graph$tv)
  triangle <- integer(Nt)
  ok <- !is.na(mesh$graph$tt)
  tt <- Matrix::sparseMatrix(
    i = c(
      which(ok[, 1]),
      which(ok[, 2]),
      which(ok[, 3])
    ),
    j = c(
      mesh$graph$tt[ok[, 1], 1],
      mesh$graph$tt[ok[, 2], 2],
      mesh$graph$tt[ok[, 3], 3]
    ),
    x = rep(1, sum(ok)),
    dims = c(Nt, Nt)
  )
  component <- 0
  while (any(triangle == 0)) {
    component <- component + 1
    tri <- integer(Nt)
    tri[min(which(triangle == 0))] <- 1
    tri_prev <- integer(Nt)
    while (any(tri_prev != tri)) {
      tri_prev <- tri
      tri <- tri | as.vector(tt %*% tri > 0)
    }
    triangle[tri > 0] <- component
    vtx <- sort(unique(as.vector(mesh$graph$tv[tri > 0, ])))
    if (any(vertex[vtx] > 0)) {
      warning(paste0(
        "Corner-only connected triangles detected.\n",
        "  Vertices = ",
        paste0(vtx[vertex[vtx] > 0], collapse = ", "), "\n",
        "  Components = ",
        paste0(vertex[vtx[vertex[vtx] > 0]], collapse = ", "), "\n",
        "  New component = ", component, "\n",
        "  Vertex component information will be inconsistent.", "\n",
        "  Triangle component information will ignore corner-only connections."
      ))
    }
    vertex[vtx] <- component
  }

  if (component == 0) {
    info <- data.frame(
      component = integer(0),
      nV = integer(0),
      nT = integer(0),
      area = numeric(0)
    )
  } else {
    fem <- fm_fem(mesh, order = 1)
    info <- do.call(
      rbind,
      lapply(
        sort(unique(triangle)),
        function(x) {
          data.frame(
            component = x,
            nV = sum(vertex == x),
            nT = sum(triangle == x),
            area = sum(fem$ta[triangle == x])
          )
        }
      )
    )
  }

  list(
    vertex = vertex,
    triangle = triangle,
    info = info
  )
}

#' @rdname fm_mesh_components
#' @export
#' @examples
#'
#' (m <- fm_mesh_3d(
#'   matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0), 4, 3, byrow = TRUE),
#'   matrix(c(1, 2, 3, 4), 1, 4, byrow = TRUE)
#' ))
#' # Compute connectivity information:
#' (conn <- fm_mesh_components(m))
fm_mesh_components.fm_mesh_3d <- function(mesh) {
  vertex <- integer(mesh$n)
  Nt <- nrow(mesh$graph$tv)
  tetra <- integer(Nt)
  ok <- !is.na(mesh$graph$tt)
  tt <- Matrix::sparseMatrix(
    i = c(
      which(ok[, 1]),
      which(ok[, 2]),
      which(ok[, 3]),
      which(ok[, 4])
    ),
    j = c(
      mesh$graph$tt[ok[, 1], 1],
      mesh$graph$tt[ok[, 2], 2],
      mesh$graph$tt[ok[, 3], 3],
      mesh$graph$tt[ok[, 4], 4]
    ),
    x = rep(1, sum(ok)),
    dims = c(Nt, Nt)
  )
  component <- 0
  while (any(tetra == 0)) {
    component <- component + 1
    tri <- integer(Nt)
    tri[min(which(tetra == 0))] <- 1
    tri_prev <- integer(Nt)
    while (any(tri_prev != tri)) {
      tri_prev <- tri
      tri <- tri | as.vector(tt %*% tri > 0)
    }
    tetra[tri > 0] <- component
    vtx <- sort(unique(as.vector(mesh$graph$tv[tri > 0, ])))
    if (any(vertex[vtx] > 0)) {
      warning(paste0(
        "Corner- or edge-only connected tetrahedrons detected.\n",
        "  Vertices = ",
        paste0(vtx[vertex[vtx] > 0], collapse = ", "), "\n",
        "  Components = ",
        paste0(vertex[vtx[vertex[vtx] > 0]], collapse = ", "), "\n",
        "  New component = ", component, "\n",
        "  Vertex component information will be inconsistent.", "\n",
        "  Tetrahedron component information will ignore corner-only ",
        "connections."
      ))
    }
    vertex[vtx] <- component
  }

  if (component == 0) {
    info <- data.frame(
      component = integer(0),
      nV = integer(0),
      nT = integer(0),
      area = numeric(0)
    )
  } else {
    fem <- fm_fem(mesh, order = 1)
    info <- do.call(
      rbind,
      lapply(
        sort(unique(tetra)),
        function(x) {
          data.frame(
            component = x,
            nV = sum(vertex == x),
            nT = sum(tetra == x),
            volume = sum(fem$ta[tetra == x])
          )
        }
      )
    )
  }

  list(
    vertex = vertex,
    tetra = tetra,
    info = info
  )
}
