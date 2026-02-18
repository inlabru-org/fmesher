#' @include deprecated.R

#' @title Generate lattice points covering a mesh
#'
#' @description Generate `terra`, `sf`, or `sp` lattice locations
#'
#' @export
#'
#' @author Finn Lindgren <Finn.Lindgren@@gmail.com>
#'
#' @param mesh An `fm_mesh_2d` object
#' @param dims A length 2 integer vector giving the dimensions of
#' the target lattice.
#' @param xlim,ylim Length 2 numeric vectors of x- and y- axis limits.
#' Defaults taken from the range of the mesh or mask; see `minimal`.
#' @param mask If logical and TRUE, remove pixels that are outside the mesh.
#' If `mask` is an `sf` or `Spatial` object, only return pixels covered by this
#' object.
#' @param format character; "sf", "terra" or "sp"
#' @param minimal logical; if `TRUE` (default), the default range is determined
#' by the minimum of the ranges of the mesh and mask, otherwise only the mesh.
#' @returns `sf`, `SpatRaster`, or `SpatialPixelsDataFrame` covering the mesh or
#' mask.
#'
#' @examples
#' if (require("ggplot2", quietly = TRUE)) {
#'   dims <- c(50, 50)
#'   pxl <- fm_pixels(
#'     fmexample$mesh,
#'     dims = dims,
#'     mask = fmexample$boundary_sf[[1]],
#'     minimal = TRUE
#'   )
#'   pxl$val <- rnorm(NROW(pxl)) +
#'     fm_evaluate(fmexample$mesh, pxl, field = 2 * fmexample$mesh$loc[, 1])
#'   ggplot() +
#'     geom_tile(
#'       data = pxl,
#'       aes(geometry = geometry, fill = val),
#'       stat = "sf_coordinates"
#'     ) +
#'     geom_sf(data = fm_as_sfc(fmexample$mesh), alpha = 0.2)
#' }
#'
#' \donttest{
#' if (require("ggplot2", quietly = TRUE) &&
#'   require("terra", quietly = TRUE) &&
#'   require("tidyterra", quietly = TRUE)) {
#'   pxl <- fm_pixels(fmexample$mesh,
#'     dims = c(50, 50), mask = fmexample$boundary_sf[[1]],
#'     format = "terra"
#'   )
#'   pxl$val <- rnorm(NROW(pxl) * NCOL(pxl))
#'   pxl <-
#'     terra::mask(
#'       pxl,
#'       mask = pxl$.mask,
#'       maskvalues = c(FALSE, NA),
#'       updatevalue = NA
#'     )
#'   ggplot() +
#'     geom_spatraster(data = pxl, aes(fill = val)) +
#'     geom_sf(data = fm_as_sfc(fmexample$mesh), alpha = 0.2)
#' }
#' }
fm_pixels <- function(mesh,
                      dims = c(150, 150),
                      xlim = NULL,
                      ylim = NULL,
                      mask = TRUE,
                      format = "sf",
                      minimal = TRUE) {
  format <- match.arg(format, c("sf", "terra", "sp"))
  if (!fm_manifold(mesh, "R2")) {
    stop("fmesher::fm_pixels() currently works for R2 meshes only.")
  }
  if (is.null(mask)) {
    mask <- FALSE
  }

  if (!is.logical(mask)) {
    if (inherits(mask, "SpatialPolygonsDataFrame")) {
      mask <- as(mask, "SpatialPolygons")
    }
    mask <- sf::st_as_sf(mask)
    mask_bbox <- sf::st_bbox(mask)
  }

  if (is.null(xlim)) {
    xlim <- range(mesh$loc[, 1])
    if (!is.logical(mask) && minimal) {
      xlim <- c(max(xlim[1], mask_bbox[1]), min(xlim[2], mask_bbox[3]))
    }
  }
  x <- seq(xlim[1], xlim[2], length.out = dims[1])
  if (is.null(ylim)) {
    ylim <- range(mesh$loc[, 2])
    if (!is.logical(mask) && minimal) {
      ylim <- c(max(ylim[1], mask_bbox[2]), min(ylim[2], mask_bbox[4]))
    }
  }
  y <- seq(ylim[1], ylim[2], length.out = dims[2])

  pixels <- expand.grid(x = x, y = y)
  pixels <- sf::st_as_sf(pixels, coords = c("x", "y"), crs = fm_crs(mesh))

  pixels_within <- rep(TRUE, NROW(pixels))
  if (is.null(mask)) {
    mask <- FALSE
  }
  if (is.logical(mask)) {
    if (mask) {
      pixels_within <- fm_is_within(pixels, mesh)
      pixels <- pixels[pixels_within, , drop = FALSE]
    }
  } else {
    if (inherits(mask, "SpatialPolygonsDataFrame")) {
      mask <- as(mask, "SpatialPolygons")
    }
    mask <- sf::st_as_sf(mask)
    pixels_within <- sf::st_covered_by(pixels, mask)
    pixels_within <- lengths(pixels_within) > 0
    pixels <- pixels[pixels_within, , drop = FALSE]
  }

  if (identical(format, "sp")) {
    pixels <- as(pixels, "Spatial")
    pixels <- as(pixels, "SpatialPixelsDataFrame")
  } else if (identical(format, "terra")) {
    fm_require_stop("terra")
    pixels <- as(pixels, "Spatial")
    pixels <- as(pixels, "SpatialPixelsDataFrame")
    pixels$.mask <- TRUE
    pixels <- terra::rast(pixels)
  }

  pixels
}


#' @title Refine a 2d mesh
#'
#' @description Refine an existing mesh
#'
#' @keywords internal
#'
#' @param mesh An [fm_mesh_2d()] object
#' @param refine A list of refinement options passed on to
#' [fm_rcdt_2d_inla]
#' @returns A refined `fm_mesh_2d` object
#' @author Finn Lindgren <Finn.Lindgren@@gmail.com>
#' @export
#' @examples
#' fm_dof(fmexample$mesh)
#' fm_dof(fm_refine(fmexample$mesh, refine = list(max.edge = 1)))
#'
fm_refine <- function(mesh, refine = list(max.edge = 1)) {
  rmesh <- fm_rcdt_2d_inla(
    loc = mesh$loc,
    boundary = fm_segm(mesh, boundary = TRUE),
    interior = fm_segm(mesh, boundary = FALSE),
    crs = fm_crs(mesh),
    refine = refine
  )
  rmesh
}


#' Split triangles of a mesh into subtriangles
#'
#' `r lifecycle::badge("experimental")`
#' Splits each mesh triangle into `(n + 1)^2` subtriangles.
#' The current version drops any edge constraint information from the mesh.
#'
#' @param mesh an [fm_mesh_2d] object
#' @param n number of added points along each edge. Default is 1.
#' @param delaunay logical; if `TRUE`, the subdivided mesh is forced into a
#'   Delaunay triangle structure. If `FALSE` (default), the triangles are
#'   subdivided uniformly instead.
#' @returns A refined [fm_mesh_2d] object, with added `bary` information
#'   (an [fm_bary()] object), that can be used for interpolating functions from
#'   the original mesh to the new mesh (from version `0.5.0.9002`).
#' @author Finn Lindgren <Finn.Lindgren@@gmail.com>
#' @export
#' @examples
#' mesh <- fm_rcdt_2d_inla(
#'   loc = rbind(c(0, 0), c(1, 0), c(0, 1)),
#'   tv = rbind(c(1, 2, 3))
#' )
#' mesh_sub <- fm_subdivide(mesh, 3)
#' mesh
#' mesh_sub
#'
#' # Difference should be zero for flat triangle meshes:
#' sum((mesh_sub$loc - fm_basis(mesh, mesh_sub$bary) %*% mesh$loc)^2)
#'
#' plot(mesh_sub, edge.color = 2)
#'
#' plot(fm_subdivide(fmexample$mesh, 3), edge.color = 2)
#' plot(fmexample$mesh, add = TRUE, edge.color = 1)
fm_subdivide <- function(mesh, n = 1, delaunay = FALSE) {
  # In case of old stored meshes, unify the graph storage:
  mesh$graph <- fm_graph(mesh)
  if (n < 1) {
    bary_index <- vapply(mesh$graph$vt, function(x) x[1, 1], integer(1))
    bary_where <- as.matrix(Matrix::sparseMatrix(
      i = seq_len(mesh$n),
      j = vapply(mesh$graph$vt, function(x) x[1, 2], integer(1)),
      x = rep(1, mesh$n),
      dims = c(mesh$n, 3)
    ))
    mesh$bary <- fm_bary(list(bary_index, bary_where))

    # Map original points, to be consistent with n >= 1:
    mesh$idx$loc <- seq_len(nrow(mesh$loc))

    return(mesh)
  }

  sub <- fmesher_subdivide(
    mesh_loc = fm_unify_coords(mesh$loc),
    mesh_tv = mesh$graph$tv - 1L,
    mesh_boundary = mesh$segm$bnd$idx - 1L,
    mesh_interior = mesh$segm$int$idx - 1L,
    subdivisions = n,
    options = list()
  )

  if (fm_manifold(mesh, "S2")) {
    radius <- mean(rowSums(mesh$loc^2)^0.5)
    renorm <- function(loc) {
      loc * (radius / rowSums(loc^2)^0.5)
    }
    sub$loc <- renorm(sub$loc)
  }

  new_mesh <- fm_rcdt_2d_inla(
    loc = sub$loc,
    tv = sub$tv + 1L,
    crs = fm_crs(mesh),
    delaunay = delaunay
  )

  bary_index <- sub$bary_index[new_mesh$idx$loc] + 1L
  bary_where <- sub$bary_where[new_mesh$idx$loc, , drop = FALSE]
  new_mesh$bary <- fm_bary(list(bary_index, bary_where))

  # Map original points:
  # Relies on the input points always being placed first, which
  # is supposed to be true.
  new_mesh$idx$loc <- seq_len(nrow(mesh$loc))

  new_mesh
}


#' Extract a subset of a mesh
#'
#' `r lifecycle::badge("experimental")` (from version `0.5.0.9003`)
#' Constructs a new mesh based on a subset of the triangles of an existing mesh.
#' The current version drops any edge constraint information from the mesh.
#'
#' @param mesh an mesh to subset
#' @param t_sub triangle or tetrahedron indices.
#' @returns A subset mesh.
#' @author Finn Lindgren <Finn.Lindgren@@gmail.com>
#' @export
#' @examples
#' mesh_sub <- fm_subset(fmexample$mesh, 1:100)
#' mesh_sub
#' plot(mesh_sub)
#'
#' if (requireNamespace("geometry", quietly = TRUE)) {
#'   print(m <- fm_delaunay_3d(matrix(rnorm(30), 10, 3)))
#'   print(fm_subset(m, seq_len(min(5, nrow(m$graph$tv)))))
#' }
fm_subset <- function(mesh, t_sub) {
  tv <- mesh$graph$tv[t_sub, , drop = FALSE]
  v <- sort(unique(as.vector(tv)))
  idx <- rep(as.integer(NA), nrow(mesh$loc))
  idx[v] <- seq_along(v)
  tv <- matrix(idx[tv], nrow(tv), ncol(tv))
  loc <- mesh$loc[v, , drop = FALSE]

  if (inherits(mesh, "fm_mesh_2d")) {
    mesh <- fmesher::fm_rcdt_2d_inla(
      loc = loc,
      tv = tv,
      refine = FALSE,
      crs = fm_crs(mesh),
      delaunay = FALSE
    )
  } else if (inherits(mesh, "fm_mesh_3d")) {
    mesh <- fmesher::fm_mesh_3d(loc = loc, tv = tv)
  } else {
    stop("`fm_subset()` currently only supports fm_mesh_2d and fm_mesh_3d.")
  }

  idx[v] <- mesh$idx$loc
  mesh$idx$loc <- idx

  mesh
}


join_segm <- function(...) {
  segm_list <- list(...)
  loc <- matrix(0, 0, 3)
  idx <- matrix(0, 0, 2)
  for (k in seq_along(segm_list)) {
    idx <- rbind(idx, segm_list[[k]]$idx + nrow(loc))
    loc <- rbind(loc, segm_list[[k]]$loc)
  }

  # Collapse duplicate points
  new_loc <- loc
  new_idx <- seq_len(nrow(loc))
  prev_idx <- 0
  for (k in seq_len(nrow(loc))) {
    if (any(is.na(new_loc[k, ]))) {
      new_idx[k] <- NA
    } else {
      if (prev_idx == 0) {
        prev_dist <- 1
      } else {
        prev_dist <- ((new_loc[seq_len(prev_idx), 1] - new_loc[k, 1])^2 +
          (new_loc[seq_len(prev_idx), 2] - new_loc[k, 2])^2 +
          (new_loc[seq_len(prev_idx), 3] - new_loc[k, 3])^2)^0.5
      }
      if (all(prev_dist > 0)) {
        prev_idx <- prev_idx + 1
        new_idx[k] <- prev_idx
        new_loc[prev_idx, ] <- new_loc[k, ]
      } else {
        new_idx[k] <- which.min(prev_dist)
      }
    }
  }
  idx <- matrix(new_idx[idx], nrow(idx), 2)
  # Remove NA and atomic lines
  ok <-
    !is.na(idx[, 1]) &
      !is.na(idx[, 2]) &
      idx[, 1] != idx[, 2]
  idx <- idx[ok, , drop = FALSE]
  # Set locations
  loc <- new_loc[seq_len(prev_idx), , drop = FALSE]

  fm_segm(
    loc = loc,
    idx = idx,
    is.bnd = FALSE
  )
}


#' Construct the intersection mesh of a mesh and a polygon
#'
#' @description `r lifecycle::badge("experimental")` (from version `0.5.0.9006`)
#'
#' @param mesh `fm_mesh_2d` object to be intersected
#' @param poly `fm_segm` object with a closed polygon to intersect with the
#'   mesh, or a polygon object that can be converted with [fm_as_segm()]
#' @returns An [fm_mesh_2d] object
#' @author Finn Lindgren <Finn.Lindgren@@gmail.com>
#' @keywords internal
#' @export
#' @examples
#' segm <- fm_segm(
#'   rbind(c(-4, -4), c(4, -3), c(0, 4)),
#'   is.bnd = TRUE
#' )
#' (m <- fm_mesh_intersection(fmexample$mesh, segm))
#' plot(fmexample$mesh)
#' lines(segm, col = 4)
#' plot(m, edge.color = 2, add = TRUE)
#'
#' \donttest{
#' # Non-overlapping addition
#' segm2 <- fm_segm(c(
#'   segm,
#'   fm_segm(
#'     rbind(c(-4, 0), c(-3, 0), c(-2, 2)),
#'     is.bnd = TRUE
#'   )
#' ))
#' (m2 <- fm_mesh_intersection(fm_subdivide(fmexample$mesh, 2), segm2))
#' m2_int <- fm_int(m2)
#' plot(m2, edge.color = 2)
#' lines(segm2, col = 4)
#' plot(fmexample$mesh, edge.color = 1, add = TRUE)
#' plot(m2_int$geometry, pch = 20, cex = sqrt(m2_int$weight) * 4, add = TRUE)
#' }
#'
#' \donttest{
#' # Add a hole and restrict to inner part of the original mesh
#' # To avoid issues with intersecting boundary segments, compute
#' # two separate intersection calculations in sequence.
#' # To allow this to be done as a single step, would need to first
#' # cross-intersect the boundary segments.
#' inner_bnd <- fm_segm(fmexample$mesh, boundary = FALSE)
#' fm_is_bnd(inner_bnd) <- TRUE
#' segm3 <- fm_segm(c(
#'   segm2,
#'   fm_segm(
#'     rbind(c(-1.5, 0), c(1, -0.5), c(0, -1.5)),
#'     is.bnd = TRUE
#'   )
#' ))
#' (m3 <- fm_mesh_intersection(
#'   fm_mesh_intersection(
#'     fm_subdivide(fmexample$mesh, 2),
#'     inner_bnd
#'   ),
#'   segm3
#' ))
#' m3_int <- fm_int(m3)
#' plot(fmexample$mesh)
#' plot(m3, edge.color = 2, add = TRUE)
#' lines(segm3, col = 4)
#' plot(m3_int$geometry, pch = 20, cex = sqrt(m3_int$weight) * 4, add = TRUE)
#' }
#'
#' \donttest{
#' # Spherical mesh
#' (m_s2 <- fm_rcdt_2d(globe = 4))
#' segm4 <- fm_segm(
#'   rbind(
#'     c(1, 0, 0.1) / sqrt(1.01),
#'     c(0, 1, 0),
#'     c(-1, -1, 1) / sqrt(3)
#'   ),
#'   is.bnd = TRUE
#' )
#' (m4 <- fm_mesh_intersection(fm_subdivide(m_s2, 1), segm4))
#' m4_int <- fm_int(m4)
#' plot(m_s2)
#' plot(m4, edge.color = 2, add = TRUE)
#' plot(m4_int$geometry, pch = 20, cex = sqrt(m4_int$weight) * 8, add = TRUE)
#' }
fm_mesh_intersection <- function(mesh, poly) {
  poly <- fm_as_segm(poly)
  poly <- fm_transform(poly, fm_crs(mesh), passthrough = TRUE)
  if (ncol(poly$loc) < 3) {
    poly$loc <- cbind(poly$loc, 0)
  }

  all_edges <- fm_segm(
    loc = mesh$loc,
    idx = cbind(
      as.vector(t(mesh$graph$tv)),
      as.vector(t(mesh$graph$tv[, c(2, 3, 1), drop = FALSE]))
    ),
    is.bnd = FALSE
  )

  mesh_cover <- fm_rcdt_2d_inla(
    loc = rbind(mesh$loc, poly$loc),
    interior = all_edges
  )

  split_segm <- fm_split_lines(mesh_cover, segm = poly)

  joint_segm <- join_segm(split_segm, all_edges)

  mesh_joint_cover <- fm_rcdt_2d_inla(
    interior = joint_segm,
    extend = TRUE
  )

  # In case the polygon is a hole, need to ensure the mesh covers the original
  # mesh. Achieved by giving it the joint mesh points as domain
  # points.
  mesh_poly <- fm_rcdt_2d_inla(
    loc = mesh_joint_cover$loc,
    boundary = split_segm
  )

  loc_tri <- fm_centroids(mesh_joint_cover)
  ok_tri <-
    fm_is_within(loc_tri, mesh) &
      fm_is_within(loc_tri, mesh_poly)
  if (any(ok_tri)) {
    mesh_subset <- fm_subset(mesh_joint_cover, which(ok_tri))
  } else {
    mesh_subset <- NULL
  }

  mesh_subset
}

#' @title Planned faster replacement for fm_mesh_intersection
#' @noRd
#' @examples
#' segm <- fm_segm(
#'   rbind(c(-4, -4), c(4, -3), c(0, 4)),
#'   is.bnd = TRUE
#' )
#' (m <- fm_intersect(fmexample$mesh, segm))
#' plot(fmexample$mesh)
#' lines(segm, col = 4)
#' plot(m, edge.color = 2, add = TRUE)
#'
#' \donttest{
#' # Non-overlapping addition
#' segm2 <- fm_segm(c(
#'   segm,
#'   fm_segm(
#'     rbind(c(-4, 0), c(-3, 0), c(-2, 2)),
#'     is.bnd = TRUE
#'   )
#' ))
#' (m2 <- fm_intersect(fm_subdivide(fmexample$mesh, 2), segm2))
#' m2_int <- fm_int(m2)
#' plot(m2, edge.color = 2)
#' lines(segm2, col = 4)
#' plot(fmexample$mesh, edge.color = 1, add = TRUE)
#' plot(m2_int$geometry, pch = 20, cex = sqrt(m2_int$weight) * 4, add = TRUE)
#' }
#'
#' \donttest{
#' # Add a hole and restrict to inner part of the original mesh
#' # To avoid issues with intersecting boundary segments, compute
#' # two separate intersection calculations in sequence.
#' # To allow this to be done as a single step, would need to first
#' # cross-intersect the boundary segments.
#' inner_bnd <- fm_segm(fmexample$mesh, boundary = FALSE)
#' fm_is_bnd(inner_bnd) <- TRUE
#' segm3 <- fm_segm(c(
#'   segm2,
#'   fm_segm(
#'     rbind(c(-1.5, 0), c(1, -0.5), c(0, -1.5)),
#'     is.bnd = TRUE
#'   )
#' ))
#' (m3 <- fm_intersect(
#'   fm_intersect(
#'     fm_subdivide(fmexample$mesh, 2),
#'     inner_bnd
#'   ),
#'   segm3
#' ))
#' mm <- fm_subdivide(fmexample$mesh, 2)
#' bench::mark(
#'   A = {
#'     (m3A <- fm_mesh_intersection(fm_mesh_intersection(mm, inner_bnd), segm3))
#'   },
#'   B = {
#'     (m3B <- fm_intersect(fm_intersect(mm, inner_bnd), segm3))
#'   },
#'   C = {
#'     (m3C <- fm_intersect(fm_intersect(mm, segm3), inner_bnd))
#'   },
#'   check = FALSE, iterations = 5
#' )
#' m3_int <- fm_int(m3)
#' plot(fmexample$mesh)
#' plot(m3, edge.color = 2, add = TRUE)
#' lines(segm3, col = 4)
#' plot(m3_int$geometry, pch = 20, cex = sqrt(m3_int$weight) * 4, add = TRUE)
#' }
#'
#' \donttest{
#' # Spherical mesh
#' (m_s2 <- fm_rcdt_2d(globe = 4))
#' segm4 <- fm_segm(
#'   rbind(
#'     c(1, 0, 0.1) / sqrt(1.01),
#'     c(0, 1, 0),
#'     c(-1, -1, 1) / sqrt(3)
#'   ),
#'   is.bnd = TRUE
#' )
#' (m4 <- fm_intersect(fm_subdivide(m_s2, 1), segm4))
#' m4_int <- fm_int(m4)
#' plot(m_s2)
#' plot(m4, edge.color = 2, add = TRUE)
#' plot(m4_int$geometry, pch = 20, cex = sqrt(m4_int$weight) * 8, add = TRUE)
#' }
fm_intersect <- function(mesh, poly) {
  poly <- fm_as_segm(poly)
  poly <- fm_transform(poly, fm_crs(mesh), passthrough = TRUE)
  if (ncol(poly$loc) < 3) {
    poly$loc <- cbind(poly$loc, 0)
  }
  if (ncol(mesh$loc) < 3) {
    mesh$loc <- cbind(mesh$loc, 0)
  }

  mesh_bnd <- fm_segm(mesh, boundary = TRUE)
  all_edges <- fm_segm(
    loc = mesh$loc,
    idx = cbind(
      as.vector(t(mesh$graph$tv)),
      as.vector(t(mesh$graph$tv[, c(2, 3, 1), drop = FALSE]))
    ),
    is.bnd = FALSE
  )

  mesh_cover <- fm_rcdt_2d_inla(
    loc = rbind(mesh$loc, poly$loc),
    interior = all_edges
  )

  split_segm <- fm_split_lines(mesh_cover, segm = poly)
  fm_is_bnd(split_segm) <- TRUE

  #  joint_segm <- join_segm(split_segm, all_edges)

  new_mesh <- fm_rcdt_2d_inla(
    boundary = join_segm(split_segm, mesh_bnd),
    interior = all_edges,
    extend = FALSE
  )

  # In case the true intersection is empty, need this additional check
  loc_tri <- fm_centroids(new_mesh)
  ok_tri <-
    fm_is_within(loc_tri, mesh) # &
  #    fm_is_within(loc_tri, mesh_poly)
  if (any(ok_tri)) {
    mesh_subset <- fm_subset(new_mesh, which(ok_tri))
  } else {
    mesh_subset <- NULL
  }

  mesh_subset
}


#' @title Store points in different formats
#'
#' @description Convert a matrix of points into different formats.
#'
#' @param loc a coordinate matrix
#' @param crs CRS information to associate with the coordinates
#' @param info An optional data.frame of additional data
#' @param format character; `"sf"`, `"df"`, `"sp"`
#' @return
#' An `sf`, `data.frame`, or `SpatialPointsDataFrame` object, with
#' optional added information.
#' @export
#' @keywords internal
#' @examples
#' fm_store_points(fmexample$loc, format = "sf")
#'
fm_store_points <- function(loc, crs = NULL, info = NULL, format = NULL) {
  format <- match.arg(
    format,
    c("sf", "df", "sp")
  )

  crs <- fm_crs(crs)

  points <- as.data.frame(loc)
  colnames(points) <- c("x", "y", "z")[seq_len(ncol(points))]
  if (!fm_crs_is_null(crs) && !fm_crs_is_geocent(crs)) {
    points <- points[, 1:2, drop = FALSE]
  }

  if (identical(format, "df")) {
    if (!is.null(info)) {
      points <- cbind(points, info)
    }
  } else if (identical(format, "sp")) {
    points <- sp::SpatialPointsDataFrame(
      points,
      data = info,
      proj4string = fm_CRS(crs)
    )
  } else if (identical(format, "sf")) {
    if (is.null(info)) {
      points <- sf::st_as_sf(
        points,
        coords = seq_len(ncol(points)),
        crs = crs
      )
    } else {
      points <- sf::st_as_sf(
        cbind(points, info),
        coords = seq_len(ncol(points)),
        crs = crs
      )
    }
  }

  points # return
}


#' @title Extract vertex locations from an `fm_mesh_2d`
#'
#' @description Extracts the vertices of an `fm_mesh_2d` object.
#'
#' @export
#' @param x An `fm_mesh_2d` object.
#' @param format character; `"sf"`, `"df"`, `"sp"`
#' @return
#' An `sf`, `data.frame`, or `SpatialPointsDataFrame` object, with the vertex
#' coordinates, and a `.vertex` column with the vertex indices.
#'
#' @author Finn Lindgren <Finn.Lindgren@@gmail.com>
#' @seealso [fm_centroids()]
#'
#' @examples
#' if (require("ggplot2", quietly = TRUE)) {
#'   vrt <- fm_vertices(fmexample$mesh, format = "sf")
#'   ggplot() +
#'     geom_sf(data = fm_as_sfc(fmexample$mesh)) +
#'     geom_sf(data = vrt, color = "red")
#' }
#'
fm_vertices <- function(x, format = NULL) {
  fm_store_points(
    loc = x$loc,
    info = data.frame(.vertex = seq_len(nrow(x$loc))),
    crs = fm_crs(x),
    format = format
  )
}

#' @title Extract triangle centroids from an `fm_mesh_2d`
#'
#' @description Computes the centroids of the triangles of an [fm_mesh_2d()]
#' object.
#'
#' @export
#' @param x An `fm_mesh_2d` object.
#' @param format character; `"sf"`, `"df"`, `"sp"`
#' @return
#' An `sf`, `data.frame`, or `SpatialPointsDataFrame` object, with the vertex
#' coordinates, and a `.triangle` column with the triangle indices.
#'
#' @author Finn Lindgren <Finn.Lindgren@@gmail.com>
#' @seealso [fm_vertices()]
#'
#' @examples
#' if (require("ggplot2", quietly = TRUE)) {
#'   vrt <- fm_centroids(fmexample$mesh, format = "sf")
#'   ggplot() +
#'     geom_sf(data = fm_as_sfc(fmexample$mesh)) +
#'     geom_sf(data = vrt, color = "red")
#' }
#'
fm_centroids <- function(x, format = NULL) {
  ## Extract triangle centroids
  loc <- (x$loc[x$graph$tv[, 1], , drop = FALSE] +
    x$loc[x$graph$tv[, 2], , drop = FALSE] +
    x$loc[x$graph$tv[, 3], , drop = FALSE]) / 3

  if (fm_manifold(x, "S2")) {
    loc <- loc / rowSums(loc^2)^0.5 * sum(x$loc[1, ]^2)^0.5
  }

  fm_store_points(
    loc = loc,
    info = data.frame(.triangle = seq_len(nrow(loc))),
    crs = fm_crs(x),
    format = format
  )
}


#' @title Add or remove Z/M information
#' @description `r lifecycle::badge("experimental")`
#' Add and/or remove Z and/or M information from simple feature geometries.
#'
#' @param x An object to modify
#' @param ... Further arguments passed to methods
#' @param add character; one of `NULL`, `"Z"`, `"M"`, or `"ZM"`. Specifies
#'   which dimensions to add.
#' @param remove character; one of `NULL`, `"Z"`, `"M"`, or `"ZM"`. Specifies
#'   which dimensions to remove.
#' @param target character; one of `"XY"`, `"XYZ"`, `"XYM"`, or `"XYZM"`.
#'   Specifies the target dimension format. If provided, overrides `add` and
#'   `remove`. When both `add` and `remove` are `NULL`, the default target is
#'   the smallest format that can hold all the inputs without loss of
#'   information.
#' @param input character or character vector; one of `NULL`, `"XY"`, `"XYZ"`,
#'   `"XYM"`, or `"XYZM"`.
#'   Specifies the input dimension format. If `NULL` (default), the input format
#'   is inferred from the number of columns in `x` (for matrices/numerics) or
#'   from the geometry type (for `sfc` objects).
#' @returns An object of the same class as `x`, with modified Z/M dimensions.
#' @seealso [sf::st_zm()] that supports a subset of these operations.
#' @author Finn Lindgren <Finn.Lindgren@@gmail.com>
#' @export
#' @rdname fm_zm
#' @examples
#' fm_zm(fmexample$loc_sf, add = "Z")
#'
fm_zm <- function(x, ...) {
  UseMethod("fm_zm")
}

#' @export
#' @rdname fm_zm
fm_zm.sf <- function(x, ...) {
  sf::st_geometry(x) <- fm_zm(sf::st_geometry(x), ...)
  x
}

#' @export
#' @rdname fm_zm
fm_zm.sfc <- function(x, ..., add = NULL, remove = NULL, target = NULL) {
  input <- fm_zm_input(x, ...)
  target <- fm_zm_target(input, add = add, remove = remove, target = target)

  sf::st_sfc(lapply(x, fm_zm, ..., target = target), crs = sf::st_crs(x))
}

#' @export
#' @rdname fm_zm
fm_zm.list <- function(x, ..., add = NULL, remove = NULL, target = NULL) {
  input <- fm_zm_input(x, ...)
  target <- fm_zm_target(input, add = add, remove = remove, target = target)

  lapply(x, fm_zm, ..., target = target)
}

#' @export
#' @rdname fm_zm
fm_zm.sfg <- function(x, ..., add = NULL, remove = NULL, target = NULL) {
  input <- fm_zm_input(x)
  target <- fm_zm_target(input, add = add, remove = remove, target = target)
  if (is.list(x)) {
    ret <- lapply(x, fm_zm, input = input, target = target)
  } else {
    ret <- fm_zm(unclass(x), input = input, target = target)
  }

  structure(ret, class = c(target, class(x)[-1]))
}

#' @export
#' @rdname fm_zm
fm_zm.numeric <- function(x,
                          ...,
                          add = NULL,
                          remove = NULL,
                          target = NULL,
                          input = NULL) {
  input <- fm_zm_input(x, ..., input = input)

  fm_zm(
    matrix(x, nrow = 1L),
    input = input,
    add = add,
    remove = remove,
    target = target
  )[1, ]
}

#' @export
#' @rdname fm_zm
fm_zm.matrix <- function(x,
                         ...,
                         add = NULL,
                         remove = NULL,
                         target = NULL,
                         input = NULL) {
  input <- fm_zm_input(x, ..., input = input)

  target <- fm_zm_target(
    input = input,
    add = add,
    remove = remove,
    target = target
  )

  ncol_for_type <- list(
    "XY" = 2L,
    "XYZ" = 3L,
    "XYM" = 3L,
    "XYZM" = 4L
  )[[target]]

  if (input == target) {
    values <- x
  } else if (input == "XY") {
    if (target == "XYZ") {
      values <- cbind(x, 0.0)
    } else if (target == "XYM") {
      values <- cbind(x, 0.0)
    } else if (input == "XY" && target == "XYZM") {
      values <- cbind(x, 0.0, 0.0)
    }
  } else if (input == "XYZ") {
    if (target == "XY") {
      values <- x[, 1:2, drop = FALSE]
    } else if (target == "XYM") {
      values <- cbind(x[, 1:2, drop = FALSE], 0.0)
    } else if (target == "XYZM") {
      values <- cbind(x, 0.0)
    }
  } else if (input == "XYM") {
    if (target == "XY") {
      values <- x[, 1:2, drop = FALSE]
    } else if (target == "XYZ") {
      values <- cbind(x[, 1:2, drop = FALSE], 0.0)
    } else if (target == "XYZM") {
      values <- cbind(x[, 1:2, drop = FALSE], 0.0, x[, 3])
    }
  } else if (input == "XYZM") {
    if (target == "XY") {
      values <- x[, 1:2, drop = FALSE]
    } else if (target == "XYZ") {
      values <- x[, 1:3, drop = FALSE]
    } else if (target == "XYM") {
      values <- x[, c(1, 2, 4), drop = FALSE]
    }
  } else {
    stop("Invalid fm_zm `input` argument")
  }

  values
}

#' @export
#' @describeIn fm_zm Find the set of distinct XY/XYZ/XYM/XYZM types
#' @examples
#' fm_zm_input(fmexample$loc_sf)
#'
fm_zm_input <- function(x, ...) {
  UseMethod("fm_zm_input")
}

#' @export
#' @rdname fm_zm
fm_zm_input.sf <- function(x, ...) {
  fm_zm_input(sf::st_geometry(x), ...)
}

#' @export
#' @rdname fm_zm
fm_zm_input.sfc <- function(x, ...) {
  unique(unlist(lapply(x, function(xx) {
    class(xx)[1]
  })))
}

#' @export
#' @rdname fm_zm
fm_zm_input.list <- function(x, ...) {
  unique(unlist(lapply(x, function(xx) {
    fm_zm_input(xx)
  })))
}

#' @export
#' @rdname fm_zm
fm_zm_input.sfg <- function(x, ...) {
  class(x)[1]
}

#' @export
#' @rdname fm_zm
fm_zm_input.numeric <- function(x, ..., input = NULL) {
  if (is.null(input)) {
    input <- c("", "XY", "XYZ", "XYZM")[length(x)]
  }
  input <- match.arg(input, c("XY", "XYZ", "XYM", "XYZM"))

  input
}

#' @export
#' @rdname fm_zm
fm_zm_input.matrix <- function(x, ..., input = NULL) {
  if (is.null(input)) {
    input <- c("", "XY", "XYZ", "XYZM")[ncol(x)]
  }
  input <- match.arg(input, c("XY", "XYZ", "XYM", "XYZM"))

  input
}

#' @describeIn fm_zm Determines the target XY/XYZ/XYM/XYZM format
#' @export
#' @examples
#' fm_zm_target(c("XY", "XYZ"))
#' fm_zm_target("XY", add = "Z")
#' fm_zm_target(c("XY", "XYZM"), remove = "M")
#'
fm_zm_target <- function(input, add = NULL, remove = NULL, target = NULL) {
  if (!is.null(target)) {
    return(target)
  }

  all_types <- c("XY", "XYZ", "XYM", "XYZM")
  compat <- list(
    "XY" = c("XY", "XYZ", "XYM", "XYZM"),
    "XYZ" = c("XYZ", "XYZM"),
    "XYM" = c("XYM", "XYZM"),
    "XYZM" = c("XYZM")
  )
  output_add <-
    list(
      "Z" = c("XYZ", "XYZM"),
      "M" = c("XYM", "XYZM"),
      "ZM" = c("XYZM")
    )
  output_rm <-
    list(
      "Z" = c("XY", "XYM"),
      "M" = c("XY", "XYZ"),
      "ZM" = c("XY")
    )
  output_compat <- all_types
  for (compat_types in compat[unique(input)]) {
    output_compat <- intersect(output_compat, compat_types)
  }

  if (is.null(target)) {
    output <- output_compat
    if (!is.null(add)) {
      output <- intersect(output, output_add[[add]])
    }
    if (!is.null(remove)) {
      output <- intersect(output, output_rm[[remove]])
      if (length(output) == 0) {
        # Will remove something
        output <- output_rm[[remove]]
        output <- output[length(output)]
      }
    }
    target <- output[1]
  } else {
    target <- target
  }

  target
}

# Convert loc information to raw matrix coordinates for the mesh
fm_onto_mesh <- function(mesh, loc, crs = NULL) {
  if (!is.matrix(loc) && !fm_crs_is_null(crs)) {
    warning("loc is non-matrix but crs specified; will be ignored")
    crs <- NULL
  }
  if (is.null(crs)) {
    crs <- fm_crs(loc)
  }
  mesh_crs <- fm_crs(mesh)

  if (inherits(loc, c("SpatialPoints", "SpatialPointsDataFrame"))) {
    fm_safe_sp(force = TRUE)
    loc <- sp::coordinates(loc)
  } else if (inherits(loc, c("sf", "sfc", "sfg"))) {
    if (inherits(loc, "sf")) {
      loc <- sf::st_geometry(loc)
    } else if (inherits(loc, "sfg")) {
      loc <- sf::st_sfc(loc)
    }
    loc <- fm_zm(loc)
    loc <- sf::st_coordinates(loc)
    c_names <- colnames(loc)
    c_names <- intersect(c_names, c("X", "Y", "Z"))
    loc <- loc[, c_names, drop = FALSE]
  } else if (inherits(loc, c("SpatVector"))) {
    loc <- terra::crds(loc)
  } else if (!is.matrix(loc)) {
    warning(
      paste0(
        "Unclear if the 'loc' class ('",
        paste0(class(loc), collapse = "', '"),
        "') is of a type we know how to handle."
      ),
      immediate. = TRUE
    )
  }

  loc_needs_normalisation <- FALSE
  if (!fm_crs_is_null(crs) && !fm_crs_is_null(mesh_crs)) {
    if (!fm_crs_is_identical(crs, mesh_crs)) {
      loc <- fm_transform(loc,
        crs = mesh_crs,
        crs0 = crs,
        passthrough = FALSE
      )
    }
  } else if (fm_manifold(mesh, "S2")) {
    loc_needs_normalisation <- TRUE
  }

  if (loc_needs_normalisation) {
    loc <- loc / rowSums(loc^2)^0.5 * mean(rowSums(mesh$loc^2)^0.5)
  }

  loc
}


#' @title Function spece degrees of freedom
#'
#' @description
#' Obtain the degrees of freedom of a function space, i.e.
#' the number of basis functions it uses.
#'
#' @param x A function space object, such as [fm_mesh_1d()] or
#' [fm_mesh_2d()]
#' @returns An integer
#' @export
#' @examples
#' fm_dof(fmexample$mesh)
#'
fm_dof <- function(x) {
  UseMethod("fm_dof")
}

#' @rdname fm_dof
#' @export
fm_dof.fm_mesh_1d <- function(x) {
  as.integer(x[["m"]])
}

#' @rdname fm_dof
#' @export
fm_dof.fm_mesh_2d <- function(x) {
  as.integer(x[["n"]])
}

#' @rdname fm_dof
#' @export
fm_dof.fm_mesh_3d <- function(x) {
  as.integer(x[["n"]])
}

#' @rdname fm_dof
#' @export
fm_dof.fm_tensor <- function(x) {
  prod(vapply(x$fun_spaces, fm_dof, 0L))
}

#' @rdname fm_dof
#' @export
fm_dof.fm_collect <- function(x) {
  sum(vapply(x$fun_spaces, fm_dof, 0L))
}

#' @rdname fm_dof
#' @export
fm_dof.fm_lattice_2d <- function(x) {
  length(x$x) * length(x$y)
}

#' @rdname fm_dof
#' @export
fm_dof.fm_lattice_Nd <- function(x) {
  prod(x$dims)
}
