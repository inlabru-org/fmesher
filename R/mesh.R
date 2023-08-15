#' @include deprecated.R

#' @title Generate lattice points covering a mesh
#'
#' @description Generate `terra`, `sf`, or `sp` lattice locations
#'
#' @export
#'
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
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
#' @param nx `r lifecycle::badge("deprecated")` Number of pixels in x direction,
#' or a numeric vector of x-values
#' @param ny `r lifecycle::badge("deprecated")` Number of pixels in y direction,
#' or a numeric vector of y-values
#' @return `sf`, `SpatRaster`, or `SpatialPixelsDataFrame` covering the mesh or
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
                      minimal = TRUE,
                      nx = deprecated(),
                      ny = deprecated()) {
  format <- match.arg(format, c("sf", "terra", "sp"))
  if (!fm_manifold(mesh, "R2")) {
    stop("fmesher::fm_pixels() currently works for R2 meshes only.")
  }
  if (is.null(mask)) {
    mask <- FALSE
  }

  x <- NULL
  if (lifecycle::is_present(nx)) {
    lifecycle::deprecate_soft(
      "0.0.1",
      "fm_pixels(nx)",
      "fm_pixels(dim)"
    )
    if (length(nx) == 1) {
      dims[1] <- nx
    } else {
      x <- nx
    }
  }
  y <- NULL
  if (lifecycle::is_present(ny)) {
    lifecycle::deprecate_soft(
      "0.0.1",
      "fm_pixels(ny)",
      "fm_pixels(dim)"
    )
    if (length(ny) == 1) {
      dims[2] <- ny
    } else {
      y <- ny
    }
  }

  if (!is.logical(mask)) {
    if (inherits(mask, "SpatialPolygonsDataFrame")) {
      mask <- as(mask, "SpatialPolygons")
    }
    mask <- sf::st_as_sf(mask)
    mask_bbox <- sf::st_bbox(mask)
  }

  if (is.null(x)) {
    if (is.null(xlim)) {
      xlim <- range(mesh$loc[, 1])
      if (!is.logical(mask) && minimal) {
        xlim <- c(max(xlim[1], mask_bbox[1]), min(xlim[2], mask_bbox[3]))
      }
    }
    x <- seq(xlim[1], xlim[2], length.out = dims[1])
  }
  if (is.null(y)) {
    if (is.null(ylim)) {
      ylim <- range(mesh$loc[, 2])
      if (!is.logical(mask) && minimal) {
        ylim <- c(max(ylim[1], mask_bbox[2]), min(ylim[2], mask_bbox[4]))
      }
    }
    y <- seq(ylim[1], ylim[2], length.out = dims[2])
  }

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
    pixels_within <- sf::st_within(pixels, mask)
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
#' @return A refined `fm_mesh_2d` object
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
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
  return(rmesh)
}



# Split triangles of a mesh into subtriangles
#
# @param mesh an inla.mesh object
# @param n number of added points along each edge
# @return A refined inla.mesh object
# @author Finn Lindgren \email{finn.lindgren@@gmail.com}
# @export

fm_subdivide <- function(mesh, n = 1) {
  if (n < 1) {
    return(mesh)
  }

  split.edges <- function(segm, n) {
    if (is.null(segm) || (nrow(segm$idx) == 0)) {
      return(segm)
    }
    n.loc <- nrow(segm$loc)
    n.idx <- nrow(segm$idx)
    loc <- do.call(
      rbind,
      c(
        list(segm$loc),
        lapply(
          seq_len(n),
          function(k) {
            (segm$loc[segm$idx[, 1], ] * k / (n + 1) +
              segm$loc[segm$idx[, 2], ] * (n - k + 1) / (n +
                1))
          }
        )
      )
    )
    idx <- do.call(
      rbind,
      c(
        list(cbind(
          segm$idx[, 1], n.loc + seq_len(n.idx)
        )),
        lapply(
          seq_len(n - 1),
          function(k) {
            cbind(
              n.loc * k + seq_len(n.idx),
              n.loc * (k + 1) + seq_len(n.idx)
            )
          }
        ),
        list(cbind(
          n.loc * n + seq_len(n.idx), segm$idx[, 2]
        ))
      )
    )

    segm2 <-
      fm_segm(
        loc = loc,
        idx = idx,
        grp = rep(segm$grp, n + 1),
        is.bnd = segm$is.bnd
      )

    segm2
  }

  p1 <- mesh$loc[mesh$graph$tv[, 1], ]
  p2 <- mesh$loc[mesh$graph$tv[, 2], ]
  p3 <- mesh$loc[mesh$graph$tv[, 3], ]

  tri.inner.loc <-
    do.call(
      rbind,
      lapply(
        seq_len(n + 2) - 1,
        function(k2) {
          do.call(
            rbind,
            lapply(
              seq_len(n + 2 - k2) - 1,
              function(k3) {
                w1 <- (n + 1 - k2 - k3) / (n + 1)
                w2 <- k2 / (n + 1)
                w3 <- k3 / (n + 1)
                p1 * w1 + p2 * w2 + p3 * w3
              }
            )
          )
        }
      )
    )

  n.tri <- nrow(p1)
  tri.edges <- fm_segm(
    loc = rbind(p1, p2, p3),
    idx = rbind(
      cbind(seq_len(n.tri), seq_len(n.tri) + n.tri),
      cbind(seq_len(n.tri) + n.tri, seq_len(n.tri) + 2 * n.tri),
      cbind(seq_len(n.tri) + 2 * n.tri, seq_len(n.tri))
    ),
    is.bnd = FALSE
  )
  # I think this code assumes fm_segm filters out duplicated points?
  new.loc <- rbind(tri.edges$loc, tri.inner.loc)

  boundary2 <- split.edges(fm_segm(mesh, boundary = TRUE), n = n)
  interior2 <- split.edges(fm_segm(mesh, boundary = FALSE), n = n)

  if (fm_manifold(mesh, "S2")) {
    radius <- mean(rowSums(mesh$loc^2)^0.5)
    renorm <- function(loc) {
      loc * (radius / rowSums(loc^2)^0.5)
    }
    new.loc <- renorm(new.loc)
    interior2$loc <- renorm(interior2$loc)
    boundary2$loc <- renorm(boundary2$loc)
  }

  mesh2 <- fm_rcdt_2d_inla(
    loc = new.loc,
    interior = interior2,
    boundary = boundary2,
    refine = list(
      min.angle = 0,
      max.edge = Inf
    ),
    crs = fm_crs(mesh)
  )

  mesh2
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
#' @param mesh `fm_mesh_2d` object to be intersected
#' @param poly `fm_segm` object with a closed polygon
#'   to intersect with the mesh
#' @returns An [fm_mesh_2d] object
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @keywords internal
#' @export
#' @examples
#' segm <- fm_segm(rbind(c(-4, -4), c(4, -4), c(0, 4)),
#'   is.bnd = TRUE
#' )
#' str(m <- fm_mesh_intersection(fmexample$mesh, segm))
#' plot(fmexample$mesh)
#' lines(segm, col = 4)
#' plot(m, edge.color = 2, add = TRUE)
fm_mesh_intersection <- function(mesh, poly) {
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

  mesh_poly <- fm_rcdt_2d_inla(boundary = poly)

  loc_tri <-
    (mesh_joint_cover$loc[mesh_joint_cover$graph$tv[, 1], , drop = FALSE] +
      mesh_joint_cover$loc[mesh_joint_cover$graph$tv[, 2], , drop = FALSE] +
      mesh_joint_cover$loc[mesh_joint_cover$graph$tv[, 3], , drop = FALSE]) / 3
  ok_tri <-
    fm_is_within(loc_tri, mesh) &
      fm_is_within(loc_tri, mesh_poly)
  if (any(ok_tri)) {
    loc_subset <- unique(sort(as.vector(mesh_joint_cover$graph$tv[ok_tri, , drop = FALSE])))
    new_idx <- integer(mesh$n)
    new_idx[loc_subset] <- seq_along(loc_subset)
    tv_subset <- matrix(new_idx[mesh_joint_cover$graph$tv[ok_tri, , drop = FALSE]],
      ncol = 3
    )
    loc_subset <- mesh_joint_cover$loc[loc_subset, , drop = FALSE]
    mesh_subset <- fm_rcdt_2d_inla(
      loc = loc_subset,
      tv = tv_subset,
      extend = FALSE
    )
  } else {
    mesh_subset <- NULL
  }

  mesh_subset
}



#' @title store points in different formats
#'
#' @description Convert a matrix of points into different formats.
#'
#' @param format character; `"sf"`, `"df"`, `"sp"`
#' @return
#' An `sf`, `data.frame`, or `SpatialPointsDataFrame` object, with
#' optional added information.
#' @export
#' @keywords internal
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
    points <- cbind(points, info)
  } else if (identical(format, "sp")) {
    points <- sp::SpatialPointsDataFrame(
      points,
      data = info,
      proj4string = fm_CRS(crs)
    )
  } else if (identical(format, "sf")) {
    points <- sf::st_as_sf(
      cbind(points, info),
      coords = seq_len(ncol(points)),
      crs = crs
    )
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
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
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
#' @param x An `fm_mesh_2d` or `inla.mesh` object.
#' @param format character; `"sf"`, `"df"`, `"sp"`
#' @return
#' An `sf`, `data.frame`, or `SpatialPointsDataFrame` object, with the vertex
#' coordinates, and a `.triangle` column with the triangle indices.
#'
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
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



# Convert loc information to raw matrix coordinates for the mesh
fm_onto_mesh <- function(mesh, loc, crs = NULL) {
  if (!is.matrix(loc) && !fm_crs_is_null(crs)) {
    warning("loc is non-matrix but crs specified; will be ignored")
  }
  if (inherits(loc, c(
    "SpatialPoints", "SpatialPointsDataFrame",
    "sf", "sfc", "sfg"
  ))) {
    crs <- fm_crs(loc)
  }
  mesh_crs <- fm_crs(mesh)

  loc_needs_normalisation <- FALSE
  if (!fm_crs_is_null(crs) && !fm_crs_is_null(mesh_crs)) {
    if (fm_crs_is_geocent(mesh_crs)) {
      if (!is.matrix(loc)) {
        if (!fm_crs_is_identical(crs, mesh_crs)) {
          loc <- fm_transform(loc, crs = mesh_crs, crs0 = crs)
        }
      }
    } else if (!fm_crs_is_identical(crs, mesh_crs)) {
      loc <- fm_transform(loc, crs = mesh_crs, crs0 = crs, passthrough = TRUE)
    }
  } else if (fm_manifold(mesh, "S2")) {
    loc_needs_normalisation <- TRUE
  }

  if (inherits(loc, c("SpatialPoints", "SpatialPointsDataFrame"))) {
    loc <- sp::coordinates(loc)
  } else if (inherits(loc, c("sf", "sfc", "sfg"))) {
    loc <- sf::st_coordinates(loc)
    c_names <- colnames(loc)
    c_names <- intersect(c_names, c("X", "Y", "Z"))
    loc <- loc[, c_names, drop = FALSE]
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
  x[["m"]]
}

#' @rdname fm_dof
#' @export
fm_dof.fm_mesh_2d <- function(x) {
  x[["n"]]
}




# Deprecated ####

#' @describeIn fmesher-deprecated Conversion to inla.mesh.segment
#' `r lifecycle::badge("deprecated")` in favour of [fm_as_segm()].
#' @returns An `fm_segm` object
#' @export
fm_as_inla_mesh_segment <-
  function(...) {
    lifecycle::deprecate_soft(
      "0.0.1",
      "fm_as_inla_mesh_segment()",
      "fm_as_segm()"
    )
    fm_as_segm(...)
  }

#' @describeIn fmesher-deprecated Conversion to inla.mesh.
#' `r lifecycle::badge("deprecated")` in favour of [fm_as_mesh_2d()].
#' @returns An `fm_mesh_2d` object
#' @export
fm_as_inla_mesh <- function(...) {
  lifecycle::deprecate_soft(
    "0.0.1",
    "fm_as_inla_mesh()",
    "fm_as_mesh_2d()"
  )
  fm_as_mesh_2d(...)
}
