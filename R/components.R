#' @title Compute connected mesh subsets
#'
#' @description Compute subsets of vertices and triangles/tetrahedrons in an
#'   [fm_mesh_2d] or [fm_mesh_3d] object that are connected by edges/triangles,
#'   and split [fm_segm] objects into connected components.
#'
#' @param x An object to extract components from
#' @param ... Additional arguments passed to methods
#' @author Finn Lindgren <Finn.Lindgren@@gmail.com>
#' @seealso [fm_mesh_2d()], [fm_rcdt_2d()], [fm_mesh_3d()], [fm_segm()]
#' @export
fm_components <- function(x, ...) {
  UseMethod("fm_components")
}

#' @describeIn fm_components Backwards compatibility for version `0.4.0`
#' @export
fm_mesh_components <- function(...) {
  lifecycle::deprecate_warn(
    "0.4.0.9001",
    "fm_mesh_components()",
    "fm_components()"
  )
  fm_components(...)
}

#' @rdname fm_components
#' @return For `fm_mesh_2d` and `fm_mesh_3d`, returns a list with elements
#' `vertex` and `triangle`/`tetra`, vectors of
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
#' conn1 <- fm_components(mesh1)
#' conn2 <- fm_components(mesh2)
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
fm_components.fm_mesh_2d <- function(x, ...) {
  mesh <- x
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

#' @rdname fm_components
#' @export
#' @examples
#'
#' (m <- fm_mesh_3d(
#'   matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0), 4, 3, byrow = TRUE),
#'   matrix(c(1, 2, 3, 4), 1, 4, byrow = TRUE)
#' ))
#' # Compute connectivity information:
#' (conn <- fm_components(m))
fm_components.fm_mesh_3d <- function(x, ...) {
  mesh <- x
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




#' @rdname fm_components
#' @returns For `fm_segm`, returns a list of segments, each with component
#'   either a single closed loop of segments, or an open segment chain.
#' @export
#' @examples
#'
#' (segm <- c(
#'   fm_segm(
#'     matrix(c(0, 0, 1, 0, 1, 1, 0, 1), 4, 2, byrow = TRUE),
#'     matrix(c(1, 2, 2, 3, 3, 4, 4, 1), 4, 2, byrow = TRUE)
#'   ),
#'   fm_segm(
#'     matrix(c(0, 0, 1, 0, 1, 1, 0, 1), 4, 2, byrow = TRUE),
#'     matrix(c(3, 4, 1, 2, 2, 3), 3, 2, byrow = TRUE),
#'     is.bnd = FALSE
#'   )
#' ))
#' # Compute connectivity information:
#' (conn <- lapply(segm, fm_components))
#' (conn2 <- fm_components(segm))
fm_components.fm_segm <- function(x, ...) {
  segm <- x
  bnd_seg <- which(segm$is.bnd)
  int_seg <- which(!segm$is.bnd)
  used_seg <- integer(0)
  comp_segments <- list()
  comp_is_closed_loop <- logical(0)
  comp_is_bnd <- logical(0)

  comp <- integer(nrow(segm$idx))
  active_comp <- 0L
  while (any(comp == 0L)) {
    active_comp <- active_comp + 1L
    comp_is_closed_loop <- c(comp_is_closed_loop, FALSE)
    curr_seg <- which.min(comp)
    curr_next_idx <- 2L
    comp_segments[[active_comp]] <- curr_seg
    comp_is_bnd <- c(comp_is_bnd, segm$is.bnd[curr_seg])
    forward <- TRUE
    local_forward <- TRUE
    repeat {
      if (comp_is_bnd[active_comp]) {
        if (forward) {
          next_seg <- which(segm$idx[, 1] == segm$idx[curr_seg, 2])
        } else {
          next_seg <- which(segm$idx[, 2] == segm$idx[curr_seg, 1])
        }
        next_seg <- intersect(
          setdiff(next_seg, used_seg),
          bnd_seg
        )
        if (length(next_seg) == 0) {
          next_seg <- NULL
        } else {
          next_seg <- min(next_seg)
        }
      } else {
        v0 <- segm$idx[curr_seg, curr_next_idx]
        curr_next_idx <- 2L
        next_seg <- intersect(
          setdiff(
            setdiff(which(segm$idx[, 1] == v0), curr_seg),
            used_seg
          ),
          int_seg
        )
        if (length(next_seg) == 0) {
          curr_next_idx <- 1L
          next_seg <- intersect(
            setdiff(
              setdiff(which(segm$idx[, 2] == v0), curr_seg),
              used_seg
            ),
            int_seg
          )
        }
        if (length(next_seg) == 0) {
          next_seg <- NULL
        } else {
          next_seg <- min(next_seg)
        }
      }
      if (length(next_seg) == 0) {
        if (forward) {
          forward <- FALSE
          curr_next_idx <- 1L
          curr_seg <- comp_segments[[active_comp]][1L]
          next
        } else {
          break
        }
      }
      if (any(next_seg %in% comp_segments[[active_comp]])) {
        comp_is_closed_loop[active_comp] <- TRUE
        if (forward) {
          sub <-
            seq(
              which.max(next_seg == comp_segments[[active_comp]]),
              length(comp_segments[[active_comp]])
            )
        } else {
          sub <-
            seq(
              1L,
              which.min(next_seg == comp_segments[[active_comp]])
            )
        }
        comp_segments[[active_comp]] <- comp_segments[[active_comp]][sub]
        break
      }
      if (forward) {
        comp_segments[[active_comp]] <-
          c(comp_segments[[active_comp]], next_seg)
      } else {
        comp_segments[[active_comp]] <-
          c(next_seg, comp_segments[[active_comp]])
      }
      curr_seg <- next_seg
    }
    comp[comp_segments[[active_comp]]] <- active_comp
    used_seg <- c(used_seg, comp_segments[[active_comp]])
  }

  if (any(comp_is_bnd & !comp_is_closed_loop)) {
    warning(paste0(
      "Some boundary segments are not closed loops.\n",
      "  Components = ",
      paste0(
        which(comp_is_bnd & !comp_is_closed_loop),
        collapse = ", "
      ), "\n",
      "  These components will be treated as interior segments."
    ))
    comp_is_bnd[comp_is_bnd & !comp_is_closed_loop] <- FALSE
  }

  fm_as_segm_list(
    lapply(seq_along(comp_segments), function(x) {
      fm_segm(
        loc = segm$loc,
        idx = segm$idx[comp_segments[[x]], , drop = FALSE],
        is.bnd = comp_is_bnd[x],
        grp = segm$grp[comp_segments[[x]]],
        crs = fm_crs(segm)
      )
    })
  )
}

#' @rdname fm_components
#' @export
fm_components.fm_segm_list <- function(x, ...) {
  do.call(c, lapply(x, fm_components))
}
