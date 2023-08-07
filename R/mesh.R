#' @include deprecated.R

#' @title Generate lattice points covering a mesh
#'
#' @description Generate `terra`, `sf`, or `sp` lattice locations
#'
#' @export
#'
#' @author Fabian E. Bachl \email{bachlfab@@gmail.com} and
#'  Finn Lindgren \email{finn.lindgren@@gmail.com}
#'
#' @param mesh An `inla.mesh` object
#' @param nx Number of pixels in x direction
#' @param ny Number of pixels in y direction
#' @param mask If logical and TRUE, remove pixels that are outside the mesh.
#' If `mask` is an `sf` or `Spatial` object, only return pixels covered by this object.
#' @param format character; "sf", "terra" or "sp"
#' @return `sf`, `SpatRaster`, or `SpatialPixelsDataFrame` covering the mesh
#'
#' @examples
#' if (require("ggplot2", quietly = TRUE) &&
#'   require("tidyterra", quietly = TRUE)) {
#'   pxl <- fm_pixels(fmexample$mesh,
#'     nx = 50, ny = 50, mask = fmexample$boundary_sf[[1]],
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
fm_pixels <- function(mesh, nx = 150, ny = 150, mask = TRUE,
                      format = "sf") {
  format <- match.arg(format, c("sf", "terra", "sp"))
  if (!fm_manifold(mesh, "R2")) {
    stop("fmesher::fm_pixels() currently works for R2 meshes only.")
  }

  if (length(nx) == 1) {
    x <- seq(min(mesh$loc[, 1]), max(mesh$loc[, 1]), length.out = nx)
  } else {
    x <- nx
  }
  if (length(ny) == 1) {
    y <- seq(min(mesh$loc[, 2]), max(mesh$loc[, 2]), length.out = ny)
  } else {
    y <- ny
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
    pixels <- as(pixels, "Spatial")
    pixels <- as(pixels, "SpatialPixelsDataFrame")
    pixels$.mask <- TRUE
    pixels <- terra::rast(pixels)
  }

  pixels
}



#' Refine a 2d mesh
#'
#' @keywords internal
#'
#' @param mesh an fm_mesh_2d object
#' @param refine A list of refinement options passed on to
#' [fm_rcdt_2d_inla]
#' @return mesh A refined fm_mesh_2d object
#' @author Fabian E. Bachl \email{bachlfab@@gmail.com}

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
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @keywords internal
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
    interior = list(all_edges)
  )

  split_segm <- fm_split_lines(mesh_cover, segm = poly)

  joint_segm <- join_segm(split_segm, all_edges)

  mesh_joint_cover <- fm_rcdt_2d_inla(
    interior = list(joint_segm),
    extend = TRUE
  )

  mesh_poly <- fm_rcdt_2d_inla(boundary = poly)

  loc_tri <-
    (mesh_joint_cover$loc[mesh_joint_cover$graph$tv[, 1], , drop = FALSE] +
      mesh_joint_cover$loc[mesh_joint_cover$graph$tv[, 2], , drop = FALSE] +
      mesh_joint_cover$loc[mesh_joint_cover$graph$tv[, 3], , drop = FALSE]) / 3
  ok_tri <-
    fm_is_within(loc = loc_tri, mesh) &
      fm_is_within(loc = loc_tri, mesh_poly)
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
#' @param x An `inla.mesh` object.
#' @param format character; `"sf"`, `"df"`, `"sp"`
#' @return
#' An `sf`, `data.frame`, or `SpatialPointsDataFrame` object, with the vertex
#' coordinates, and a `.vertex` column with the vertex indices.
#'
#' @author Fabian E. Bachl \email{bachlfab@@gmail.com},
#' Finn Lindgren \email{finn.lindgren@@gmail.com}
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



# fm_bary ####

#' @title Compute barycentric coordinates
#'
#' @description Identify knot intervals or triangles and compute barycentric coordinates
#'
#' @param mesh `fm_mesh_1d` or `fm_mesh_2d` object
#' @param loc Points for which to identify the containing triangle, and
#' corresponding barycentric coordinates. May be a vector (for 1d) or
#' raw matrix coordinates, `sf`, or `sp` point information (for 2d).
#' @param \dots Arguments forwarded to sub-methods.
#'
#' @export
fm_bary <- function(mesh, loc, ...) {
  UseMethod("fm_bary")
}


#' @describeIn fm_bary Return a list with elements
#' `t` (start and endpoint knot indices) and `bary` (barycentric coordinates), both
#' 2-column matrices. For backwards compatibility with old inla code, a copy `index=t`
#' is also included in the list.
#' @param method character; method for defining the barycentric coordinates
#' @export
fm_bary.fm_mesh_1d <- function(mesh, loc, method = c("linear", "nearest"), ...) {
  method <- match.arg(method)

  if (method == "linear") {
    if (mesh$cyclic) {
      mloc <- c(mesh$loc - mesh$loc[1], diff(mesh$interval))
      loc <- (loc - mesh$loc[1]) %% diff(mesh$interval)
    } else {
      mloc <- c(mesh$loc - mesh$loc[1], diff(mesh$interval))
      loc <- pmax(0, pmin(diff(mesh$interval), loc - mesh$loc[1]))
    }
  } else {
    if (mesh$cyclic) {
      mloc <-
        c(
          mesh$loc[mesh$n] - diff(mesh$interval),
          mesh$loc,
          diff(mesh$interval)
        )
      mloc <- (mloc[-(mesh$n + 2)] + mloc[-1]) / 2
      loc <- (loc - mloc[1]) %% diff(mesh$interval)
      mloc <- mloc - mloc[1]
    } else {
      mloc <-
        c(
          0,
          (mesh$loc[1:(mesh$n - 1L)] +
            mesh$loc[2:mesh$n]) / 2 - mesh$loc[1],
          diff(mesh$interval)
        )
      loc <- pmax(0, pmin(diff(mesh$interval), loc - mesh$loc[1]))
    }
  }

  ## Binary split method:
  do.the.split <- function(knots, loc) {
    n <- length(knots)
    if (n <= 2L) {
      return(rep(1L, length(loc)))
    }
    split <- 1L + (n - 1L) %/% 2L ## Split point
    upper <- (loc >= knots[split])
    idx <- rep(0, length(loc))
    idx[!upper] <- do.the.split(knots[1:split], loc[!upper])
    idx[upper] <- split - 1L + do.the.split(knots[split:n], loc[upper])
    return(idx)
  }

  idx <- do.the.split(mloc, loc)

  if (method == "nearest") {
    u <- rep(0, length(loc))
    if (mesh$cyclic) {
      found <- which(idx == (mesh$n + 1L))
      idx[found] <- 1L
    }
  } else { ## (method=="linear") {
    u <- pmax(0, pmin(1, (loc - mloc[idx]) / (mloc[idx + 1L] - mloc[idx])))
    if (!mesh$cyclic) {
      found <- which(idx == mesh$n)
      idx[found] <- mesh$n - 1L
      u[found] <- 1
    }
  }

  if (mesh$cyclic) {
    index <- matrix(c(idx, (idx %% mesh$n) + 1L), length(idx), 2)
    bary <- matrix(c(1 - u, u), length(idx), 2)
  } else {
    index <- matrix(c(idx, idx + 1L), length(idx), 2)
    bary <- matrix(c(1 - u, u), length(idx), 2)
  }

  return(list(t = index, bary = bary, index = index))
}


#' @param crs Optional crs information for `loc`
#'
#' @describeIn fm_bary A list with elements `t` (vector of triangle indices) and `bary`
#' (3-column matrix of barycentric coordinates). Points that were not found
#' give `NA` entries in `t` and `bary`.
#'
#' @export
fm_bary.fm_mesh_2d <- function(mesh, loc, crs = NULL, ...) {
  loc <- fm_onto_mesh(mesh, loc, crs = crs)

  # Avoid sphere accuracy issues by scaling to unit sphere
  scale <- 1
  if (fm_manifold(mesh, "S2")) {
    scale <- 1 / mean(rowSums(mesh$loc^2)^0.5)
    loc <- loc / rowSums(loc^2)^0.5
  }

  pre_ok_idx <-
    which(rowSums(matrix(
      is.na(as.vector(loc)),
      nrow = nrow(loc),
      ncol = ncol(loc)
    )) == 0)
  result <- fmesher_bary(
    mesh_loc = mesh$loc * scale,
    mesh_tv = mesh$graph$tv - 1L,
    loc = loc[pre_ok_idx, , drop = FALSE],
    options = list()
  )
  tri <- rep(NA_integer_, nrow(loc))
  bary <- matrix(NA_real_, nrow(loc), 3)
  ok <- result$t >= 0
  tri[pre_ok_idx[ok]] <- result$t[ok] + 1L
  bary[pre_ok_idx[ok], ] <- result$bary[ok, ]
  list(t = tri, bary = bary)
}



#' @rdname fm_bary
#' @export
#' @method fm_bary inla.mesh
fm_bary.inla.mesh <- function(mesh, ...) {
  fm_bary.fm_mesh_2d(fm_as_mesh_2d(mesh), ...)
}

#' @rdname fm_bary
#' @export
#' @method fm_bary inla.mesh.1d
fm_bary.inla.mesh.1d <- function(mesh, ...) {
  fm_bary.fm_mesh_1d(fm_as_mesh_1d(mesh), ...)
}



# fm_fem ####

#' @title Compute finite element matrices
#'
#' @description (...)
#'
#' @param mesh `inla.mesh` object
#' @param order integer
#' @param ... Currently unused
#'
#' @return A list with elements `...`
#'
#' @export
fm_fem <- function(mesh, order = 2, ...) {
  UseMethod("fm_fem")
}

#' @rdname fm_fem
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
#' @export
fm_fem.fm_mesh_2d <- function(mesh, order = 2,
                              aniso = NULL,
                              ...) {
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



# fm_split_lines ####

#' @title Split lines at triangle edges
#'
#' @description Compute intersections between line segments and triangle edges,
#' and filter out segment of length zero.
#'
#' @param mesh An `fm_mesh_2d` or `inla.mesh` object
#' @param segm An [fm_segm()] object with segments to be split
#' @param ... Unused.
#' @return An [fm_segm()] object with the same crs as the mesh,
#' with an added field `origin`, that for each new segment gives the
#' originator index into to original `segm` object for each new line segment.
#'
#' @author Fabian E. Bachl \email{f.e.bachl@@bath.ac.uk}
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#'
#' @export
fm_split_lines <- function(mesh, ...) {
  UseMethod("fm_split_lines")
}


#' @rdname fm_split_lines
#' @export
fm_split_lines.fm_mesh_2d <- function(mesh, segm, ...) {
  segm <-
    fm_transform(
      fm_as_segm(segm),
      crs = fm_crs(mesh),
      crs0 = fm_crs(segm),
      passthrough = TRUE
    )
  origin <- seq_len(NROW(segm$idx))
  if (NROW(segm$loc) > 0) {
    # Filter out segments not on the mesh
    t1 <- fm_bary(mesh, loc = segm$loc, crs = fm_crs(segm))$t
    keep <- !(is.na(t1[segm$idx[, 1]]) | is.na(t1[segm$idx[, 2]]))
    # if (any(!keep)) { warning("points outside boundary! filtering...")}
    segm <- fm_segm(
      loc = segm$loc,
      idx = segm$idx[keep, , drop = FALSE],
      grp = segm$grp[keep],
      is.bnd = segm$is.bnd,
      crs = fm_crs(segm)
    )
    origin <- origin[keep]
  }

  if (NROW(segm$idx) == 0) {
    segm$origin <- origin
    return(segm)
  }

  # Split the segments into parts
  if (NCOL(segm$loc) == 2) {
    segm$loc <- cbind(segm$loc, rep(0, NROW(segm$loc)))
  }
  splt <- fmesher_split_lines(
    mesh_loc = mesh$loc,
    mesh_tv = mesh$graph$tv - 1L,
    loc = segm$loc,
    idx = segm$idx - 1L,
    options = list()
  )
  indexoutput <- list("split.idx", "split.t", "split.origin")
  for (name in intersect(names(splt), indexoutput)) {
    splt[[name]] <- splt[[name]] + 1L
  }

  segm.split <- fm_segm(
    loc = splt$split.loc,
    idx = splt$split.idx,
    grp = segm$grp[splt$split.origin],
    is.bnd = segm$is.bnd,
    crs = fm_crs(segm)
  )
  origin <- origin[splt$split.origin]

  #  plot(mesh)
  #  lines(segm)
  #  points(segm.split$loc[, 1:2], col="red", pch = 20)
  #  points(segm$loc[, 1:2], col="blue", pch = 20)

  # Filter out zero length segments
  keep <- rowSums((segm.split$loc[segm.split$idx[, 2], , drop = FALSE] -
    segm.split$loc[segm.split$idx[, 1], , drop = FALSE])^2) > 0
  segm.split <- fm_segm(
    loc = segm.split$loc,
    idx = segm.split$idx[keep, , drop = FALSE],
    grp = segm.split$grp[keep],
    is.bnd = segm.split$is.bnd,
    crs = fm_crs(segm)
  )
  segm.split$origin <- origin[keep]

  return(segm.split)
}

#' @rdname fm_split_lines
#' @export
fm_split_lines.inla.mesh <- function(mesh, ...) {
  fm_split_lines(fm_as_mesh_2d(mesh), ...)
}





# fm_nonconvex_hull_inla ####

#' @title Recursive curve simplification.
#'
#' @description
#' Helper from legacy `INLA::inla.simplify.curve()`
#'
#' Attempts to simplify a polygonal curve by joining nearly colinear segments.
#'
#' Uses a variation of the binary splitting Ramer-Douglas-Peucker algorithm,
#' with an ellipse of half-width `eps` ellipse instead of a rectangle, motivated by
#' prediction ellipse for Brownian bridge.
#'
#' @param loc Coordinate matrix.
#' @param idx Index vector into `loc` specifying a polygonal curve.
#' @param eps Absolute straightness tolerance. Default `NULL`, no constraint.
#' @param eps_rel Relative straightness tolerance. Default `NULL`, no constraint.
#' @return An index vector into `loc` specifying the simplified polygonal
#' curve.
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @details
#' Variation of Ramer-Douglas-Peucker.
#' Uses width epsilon ellipse instead of rectangle,
#' motivated by prediction ellipse for Brownian bridge.
#'
#' @examples
#'
#' theta <- seq(0, 2 * pi, length.out = 1000)
#' loc <- cbind(cos(theta), sin(theta))
#' idx <- fm_simplify_helper(loc = loc, idx = 1:nrow(loc), eps = 0.01)
#' print(c(nrow(loc), length(idx)))
#' plot(loc, type = "l")
#' lines(loc[idx, ], col = "red")
#' @export
#' @keywords internal
#' @family nonconvex inla legacy support
fm_simplify_helper <- function(loc, idx, eps = NULL, eps_rel = NULL) {
  n <- length(idx)
  if ((n <= 2) ||
    (is.null(eps) && is.null(eps_rel)) ||
    (min(eps, eps_rel) == 0)) {
    return(idx)
  }
  segm <- loc[idx[n], ] - loc[idx[1], ]
  segm.len <- sum(segm^2)^0.5
  if (segm.len <= 1e-12) {
    ## End point same as start; closed curve.  Split.
    len2 <- ((loc[idx[2:(n - 1)], 1] - loc[idx[1], 1])^2 +
      (loc[idx[2:(n - 1)], 2] - loc[idx[1], 2])^2)
    split <- which.max(len2) + 1L
  } else {
    segm.mid <- (loc[idx[n], ] + loc[idx[1], ]) / 2
    segm <- segm / segm.len
    segm.perp <- c(-segm[2], segm[1])
    vec <- (cbind(
      loc[idx[2:(n - 1)], 1] - segm.mid[1],
      loc[idx[2:(n - 1)], 2] - segm.mid[2]
    ))
    ## Always split if any point is outside the circle
    epsi <- min(c(eps, eps_rel * segm.len / 2, segm.len / 2))
    dist1 <- abs(vec[, 1] * segm[1] + vec[, 2] * segm[2]) / (segm.len / 2) * epsi
    dist2 <- abs(vec[, 1] * segm.perp[1] + vec[, 2] * segm.perp[2])
    dist <- (dist1^2 + dist2^2)^0.5

    ## Find the furthest point, in the ellipse metric, and
    ## check if it inside the ellipse (radius=segm.len/2)
    split <- which.max(dist) + 1L
    if (dist[split - 1L] < epsi) {
      ## Flat segment, eliminate.
      return(idx[c(1, n)])
    }
    ## Split at the furthest point.
    split <- which.max(dist) + 1L
  }

  ## Do the split recursively:
  return(c(
    fm_simplify_helper(loc, idx[1L:split], eps = eps, eps_rel = eps_rel),
    fm_simplify_helper(loc, idx[split:n], eps = eps, eps_rel = eps_rel)[-1L]
  ))
}

#' @title Recursive curve simplification.
#'
#' @description `r lifecycle::badge("experimental")`
#' Simplifies polygonal curve segments by joining nearly
#' co-linear segments.
#'
#' Uses a variation of the binary splitting Ramer-Douglas-Peucker algorithm,
#' with an ellipse of half-width `eps` ellipse instead of a rectangle, motivated by
#' prediction ellipse for Brownian bridge.
#'
#' @param x An [fm_segm()] object.
#' @param eps Absolute straightness tolerance. Default `NULL`, no constraint.
#' @param eps_rel Relative straightness tolerance. Default `NULL`, no constraint.
#' @param ... Currently unused.
#' @return The simplified [fm_segm()] object.
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @details
#' Variation of Ramer-Douglas-Peucker.
#' Uses width epsilon ellipse instead of rectangle,
#' motivated by prediction ellipse for Brownian bridge.
#' @references
#' Ramer, Urs (1972). "An iterative procedure for the polygonal approximation
#' of plane curves". *Computer Graphics and Image Processing*. **1** (3): 244–256.
#' [doi:10.1016/S0146-664X(72)80017-0](https://doi.org/10.1016/S0146-664X(72)80017-0)
#'
#' Douglas, David; Peucker, Thomas (1973). "Algorithms for the reduction of
#' the number of points required to represent a digitized line or its
#' caricature". *The Canadian Cartographer*. **10** (2): 112–122.
#' [doi:10.3138/FM57-6770-U75U-7727](https://doi.org/10.3138/FM57-6770-U75U-7727)
#'
#' @examples
#' theta <- seq(0, 2 * pi, length.out = 1000)
#' (segm <- fm_segm(cbind(cos(theta), sin(theta)),
#'   idx = seq_along(theta)
#' ))
#' (segm1 <- fm_simplify(segm, eps_rel = 0.1))
#' (segm2 <- fm_simplify(segm, eps_rel = 0.2))
#' plot(segm)
#' lines(segm1, col = 2)
#' lines(segm2, col = 3)
#'
#' (segm <- fm_segm(cbind(theta, sin(theta * 4)),
#'   idx = seq_along(theta)
#' ))
#' (segm1 <- fm_simplify(segm, eps_rel = 0.1))
#' (segm2 <- fm_simplify(segm, eps_rel = 0.2))
#' plot(segm)
#' lines(segm1, col = 2)
#' lines(segm2, col = 3)
#' @export
#' @family object creation and conversion
fm_simplify <- function(x, eps = NULL, eps_rel = NULL, ...) {
  if ((nrow(x$idx) <= 1) ||
    (is.null(eps) && is.null(eps_rel)) ||
    (min(eps, eps_rel) == 0)) {
    return(x)
  }

  segm_split <- list()
  k <- 0L

  not_handled_seg <- seq_len(nrow(x$idx))
  while (length(not_handled_seg) > 0) {
    next_seg <- not_handled_seg[1]
    seq_seg <- integer(0)
    seq_vtx <- x$idx[next_seg, 1]
    next_vtx <- x$idx[next_seg, 2]
    while (TRUE) {
      final <- next_vtx %in% seq_vtx
      seq_vtx <- c(seq_vtx, next_vtx)
      seq_seg <- c(seq_seg, next_seg)
      not_handled_seg <- setdiff(not_handled_seg, next_seg)
      if (final) {
        break
      }
      if (x$is.bnd) {
        next_seg <- not_handled_seg[which(x$idx[not_handled_seg, 1] == next_vtx)]
        if (length(next_seg) == 0) {
          break
        }
        next_seg <- next_seg[1]
        next_vtx <- x$idx[next_seg, 2]
      } else {
        next_seg1 <- not_handled_seg[which(x$idx[not_handled_seg, 1] == next_vtx)]
        next_seg2 <- not_handled_seg[which(x$idx[not_handled_seg, 2] == next_vtx)]
        if ((length(next_seg1) == 0) && (length(next_seg2) == 0)) {
          break
        }
        if (length(next_seg1) > 0) {
          next_seg <- next_seg1[1]
          next_vtx <- x$idx[next_seg, 2]
        } else {
          next_seg <- next_seg2[1]
          next_vtx <- x$idx[next_seg, 1]
        }
      }
    }

    # seq_vtx and seq_seg
    k <- k + 1
    # TODO: handle geocent data
    idx <- fm_simplify_helper(
      loc = x$loc,
      idx = seq_vtx,
      eps = eps,
      eps_rel = eps_rel
    )
    # TODO: improve granularity of group information.
    segm_split[[k]] <- fm_segm(
      loc = x$loc,
      idx = idx,
      grp = x$grp[seq_seg[1]],
      is.bnd = x$is.bnd,
      crs = fm_crs(x)
    )
  }
  segm_split <- fm_segm_join(fm_as_segm_list(segm_split))

  segm_split
}


#' Contour segment
#'
#' Helper from legacy `INLA::inla.contour.segment()`
#'
#' @export
#' @keywords internal
#' @family nonconvex inla legacy support
fm_segm_contour_helper <- function(x = seq(0, 1, length.out = nrow(z)),
                                   y = seq(0, 1, length.out = ncol(z)),
                                   z, nlevels = 10,
                                   levels = pretty(range(z, na.rm = TRUE), nlevels),
                                   groups = seq_len(length(levels)),
                                   positive = TRUE,
                                   eps = NULL,
                                   eps_rel = NULL,
                                   crs = NULL) {
  ## Input checking from contourLines:
  if (missing(z)) {
    if (!missing(x)) {
      if (is.list(x)) {
        z <- x$z
        y <- x$y
        x <- x$x
      } else {
        z <- x
        x <- seq.int(0, 1, length.out = nrow(z))
      }
    } else {
      stop("no 'z' matrix specified")
    }
  } else if (is.list(x)) {
    y <- x$y
    x <- x$x
  }

  if (is.null(eps)) {
    eps <- min(c(min(diff(x)), min(diff(y)))) / 8
  }
  ## End of input checking.

  ## Get contour pieces
  curves <- grDevices::contourLines(x, y, z, levels = levels)

  ## Make a mesh for easy gradient interpolation:
  latt <- fm_lattice_2d(x, y)
  mesh <-
    fm_rcdt_2d(
      lattice = latt,
      boundary = latt$segm,
      extend = (list(
        n = 3,
        offset = (max(
          diff(range(x)),
          diff(range(y))
        ) * 0.1)
      ))
    )
  ## Map function values to mesh indexing:
  zz <- rep(0, mesh$n)
  ## 2020-07-15:
  ## Bug in fmesher may lead to some lattice points missing. Ignore those.
  lattice_points_ok <- !is.na(mesh$idx$lattice)
  zz[mesh$idx$lattice[lattice_points_ok]] <- as.vector(z)[lattice_points_ok]

  ## Mapping from level to group value:
  level2grp <- function(level) {
    if (length(groups) == 1) {
      return(groups)
    }
    for (k in seq_along(groups)) {
      if (levels[k] == level) {
        return(groups[k])
      }
    }
    return(0)
  }

  ## Join all contour pieces into a single mesh.segment
  ## Different levels can later be identified via the grp indices.
  loc <- matrix(0, 0, 2)
  idx <- matrix(0, 0, 2)
  grp <- c()
  for (k in seq_len(length(curves))) {
    curve.loc <- cbind(curves[[k]]$x, curves[[k]]$y)
    curve.n <- nrow(curve.loc)

    ## Extract the rotated gradients along the curve
    curve.mid <- (curve.loc[1:(curve.n - 1), ] + curve.loc[2:curve.n, ]) / 2
    A <- fm_basis(mesh, loc = curve.mid, derivatives = TRUE)
    ## Gradients rotated 90 degrees CW, i.e. to the direction
    ## of CCW curves around positive excursions:
    grid.diff <- cbind(A$dy %*% zz, -A$dx %*% zz)

    ## Determine the CCW/CW orientation
    curve.diff <- diff(curve.loc)
    ccw <- (sum(curve.diff * grid.diff) >= 0) ## True if in CCW direction
    if ((ccw && positive) || (!ccw && !positive)) {
      curve.idx <- seq_len(curve.n)
    } else {
      curve.idx <- rev(seq_len(curve.n))
    }

    ## Filter short line segments:
    curve.idx <-
      fm_simplify_helper(curve.loc,
        curve.idx,
        eps = eps,
        eps_rel = eps_rel
      )

    ## Reorder, making sure any unused points are removed:
    curve.loc <- curve.loc[curve.idx, , drop = FALSE]
    curve.n <- nrow(curve.loc)
    curve.idx <- cbind(seq_len(curve.n - 1L), seq_len(curve.n - 1L) + 1L)

    ## Check if the curve is closed, and adjust if it is:
    if (max(abs(curve.loc[1, ] - curve.loc[curve.n, ])) < 1e-12) {
      curve.loc <- curve.loc[-curve.n, , drop = FALSE]
      curve.n <- nrow(curve.loc)
      curve.idx <-
        cbind(seq_len(curve.n), c(seq_len(curve.n - 1L) + 1L, 1L))
    }

    ## Add the curve:
    offset <- nrow(loc)
    loc <- rbind(loc, curve.loc)
    idx <- rbind(idx, curve.idx + offset)
    grp <- c(grp, rep(level2grp(curves[[k]]$level), curve.n - 1L))
  }

  fm_segm(loc = loc, idx = idx, grp = grp, is.bnd = FALSE, crs = crs)
}

#' @title Non-convex hull computation
#' @description Legacy method for `INLA::inla.nonconvex.hull()`
#' @seealso [fm_nonconvex_hull()]
#' @param resolution The internal computation resolution.  A warning will be
#' issued when this needs to be increased for higher accuracy, with the
#' required resolution stated.
#' @param eps,eps_rel The polygonal curve simplification tolerances used for
#' simplifying the resulting boundary curve.  See [fm_simplify_helper()] for
#' details.
#' @param \dots Unused.
#' @inheritParams fm_nonconvex_hull
#' @returns `fm_nonconvex_hull_inla()` returns an `fm_segm`/`inla.mesh.segment`
#' object, for compatibility with `inla.nonconvex.hull()`.
#' @export
#' @family nonconvex inla legacy support
fm_nonconvex_hull_inla <- function(x,
                                   convex = -0.15,
                                   concave = convex,
                                   resolution = 40,
                                   eps = NULL,
                                   eps_rel = NULL,
                                   crs = NULL,
                                   ...) {
  stopifnot(!is.null(x))
  if (inherits(x, c("SpatialPoints", "SpatialPointsDataFrame"))) {
    x <- fm_transform(
      sp::coordinates(x),
      crs0 = fm_crs(x),
      crs = fm_crs(crs),
      passthrough = TRUE
    )
  } else if (inherits(x, c("sf", "sfc"))) {
    x <- fm_transform(
      sf::st_coordinates(x),
      crs0 = fm_crs(x),
      crs = fm_crs(x),
      passthrough = TRUE
    )
  }

  if (length(resolution) == 1) {
    resolution <- rep(resolution, 2)
  }
  lim <- rbind(range(x[, 1]), range(x[, 2]))

  approx.diam <- max(diff(lim[1, ]), diff(lim[2, ]))
  if (convex < 0) {
    convex <- -convex * approx.diam
  }
  if (concave < 0) {
    concave <- -concave * approx.diam
  }
  if (concave == 0) {
    return(fm_nonconvex_hull_inla_basic(x, convex, resolution, eps,
      crs = crs
    ))
  }

  ex <- convex + concave
  domain <- c(diff(lim[1, ]), diff(lim[2, ])) + 2 * ex
  dif <- domain / (resolution - 1)
  if (max(dif) > min(convex, concave)) {
    req.res <- ceiling(domain / min(convex, concave) + 1)
    warning(paste("Resolution (",
      paste(resolution, collapse = ","),
      ") too small for convex/concave radius (",
      convex, ",", concave,
      ").\n",
      "Resolution >=(",
      paste(req.res, collapse = ","),
      ") required for more accurate results.",
      sep = ""
    ))
  }
  ax <-
    list(
      seq(lim[1, 1] - ex, lim[1, 2] + ex, length.out = resolution[1]),
      seq(lim[2, 1] - ex, lim[2, 2] + ex, length.out = resolution[2])
    )
  xy <- as.matrix(expand.grid(ax[[1]], ax[[2]]))

  requireNamespace("splancs")
  z <- (matrix(
    splancs::nndistF(x, xy),
    resolution[1], resolution[2]
  ))
  segm.dilation <-
    fm_segm_contour_helper(
      ax[[1]], ax[[2]], z,
      levels = c(convex + concave),
      positive = TRUE,
      eps = 0
    ) ## Don't simplify curve at this stage
  mesh.dilation <-
    fm_rcdt_2d(
      loc = xy,
      boundary = segm.dilation,
      extend = (list(
        n = 3,
        offset = (max(
          diff(ax[[1]]),
          diff(ax[[2]])
        ) * 0.1)
      ))
    )

  z <- (matrix(
    splancs::nndistF(mesh.dilation$loc, xy),
    resolution[1], resolution[2]
  ))
  segm.closing <-
    fm_segm_contour_helper(
      ax[[1]], ax[[2]], z,
      levels = c(concave),
      positive = TRUE,
      eps = eps,
      eps_rel = eps_rel
    )

  segm.closing$crs <- crs

  fm_as_segm(segm.closing)
}



#' @details Requires `splancs::nndistF()`
#' @export
#' @describeIn fm_nonconvex_hull_inla Special method for `convex = 0`.
## Based on an idea from Elias Teixeira Krainski
fm_nonconvex_hull_inla_basic <- function(x, convex = -0.15, resolution = 40,
                                         eps = NULL, crs = NULL) {
  stopifnot(!is.null(x))
  if (inherits(x, c("SpatialPoints", "SpatialPointsDataFrame"))) {
    x <- fm_transform(
      sp::coordinates(x),
      crs0 = fm_crs(x),
      crs = fm_crs(crs),
      passthrough = TRUE
    )
  } else if (inherits(x, c("sf", "sfc"))) {
    x <- fm_transform(
      sf::st_coordinates(x),
      crs0 = fm_crs(x),
      crs = fm_crs(x),
      passthrough = TRUE
    )
  }

  if (length(convex) == 1) {
    convex <- rep(convex, 2)
  }
  if (length(resolution) == 1) {
    resolution <- rep(resolution, 2)
  }

  lim <- rbind(range(x[, 1]), range(x[, 2]))
  ex <- convex
  if (convex[1] < 0) {
    ex[1] <- -convex[1] * diff(lim[1, ])
  }
  if (convex[2] < 0) {
    ex[2] <- -convex[2] * diff(lim[2, ])
  }

  domain <- c(diff(lim[1, ]), diff(lim[2, ])) + 2 * ex
  dif <- domain / (resolution - 1)
  if (any(dif > min(convex))) {
    req.res <- ceiling(domain / convex + 1)
    warning(paste("Resolution (",
      paste(resolution, collapse = ","),
      ") too small for convex (",
      paste(convex, collapse = ","),
      ").\n",
      "Resolution >=(",
      paste(req.res, collapse = ","),
      ") required for more accurate results.",
      sep = ""
    ))
  }

  ax <- list(
    seq(lim[1, 1] - ex[1], lim[1, 2] + ex[1], length.out = resolution[1]),
    seq(lim[2, 1] - ex[2], lim[2, 2] + ex[2], length.out = resolution[2])
  )
  xy <- as.matrix(expand.grid(ax[[1]], ax[[2]]))
  tr <- diag(c(1 / ex[1], 1 / ex[2]), nrow = 2, ncol = 2)

  requireNamespace("splancs")
  z <- matrix(
    splancs::nndistF(x %*% tr, xy %*% tr),
    resolution[1], resolution[2]
  )
  segm <- fm_segm_contour_helper(
    ax[[1]],
    ax[[2]],
    z,
    levels = c(1),
    positive = FALSE,
    eps = eps,
    crs = crs
  )
  return(segm)
}




# fm_nonconvex_hull ####

#' @title Compute an extension of a spatial object
#'
#' @description
#' Constructs a potentially nonconvex extension of a spatial object by
#' performing dilation by `convex + concave` followed by
#' erosion by `concave`. This is equivalent to dilation by `convex` followed
#' by closing (dilation + erosion) by `concave`.
#'
#' @describeIn fm_nonconvex_hull Basic nonconvex hull method.
#'
#' @details
#' Morphological dilation by `convex`, followed by closing by
#' `concave`, with minimum concave curvature radius `concave`.  If
#' the dilated set has no gaps of width between \deqn{2 \text{convex} (\sqrt{1+2
#' \text{concave}/\text{convex}} - 1)}{2*convex*(sqrt(1+2*concave/convex) - 1)}
#' and \eqn{2\text{concave}}{2*concave}, then the minimum convex curvature radius is
#' `convex`.
#'
#' The implementation is based on the identity \deqn{\text{dilation}(a) \&
#' \text{closing}(b) = \text{dilation}(a+b) \& \text{erosion}(b)}{
#' dilation(a) & closing(b) = dilation(a+b) & erosion(b)} where all operations
#' are with respect to disks with the specified radii.
#'
#' @param x A spatial object
#' @param ... Arguments passed on to the [fm_nonconvex_hull()] sub-methods
#' @param convex How much to extend
#' @param concave The minimum allowed reentrant curvature. Default equal to `convex`
#' @param preserveTopology logical; argument to `sf::st_simplify()`
#' @param dTolerance If not null, the `dTolerance` argument to `sf::st_simplify()`
#' @param crs Options crs object for the resulting polygon
#' @returns `fm_nonconvex_hull()` returns an extended object as an `sfc`
#' polygon object (regardless of the `x` class).
#' @references Gonzalez and Woods (1992), Digital Image Processing
#' @seealso [fm_nonconvex_hull_inla()]
#' @export
#' @examples
#' inp <- matrix(rnorm(20), 10, 2)
#' out <- fm_nonconvex_hull(inp, convex = 1)
#' plot(out)
#' points(inp, pch = 20)
fm_nonconvex_hull <- function(x, ...) {
  UseMethod("fm_nonconvex_hull")
}

#' @rdname fm_nonconvex_hull
#' @details Differs from `sf::st_buffer(x, convex)` followed by
#' `sf::st_concave_hull()` (available from GEOS 3.11)
#' in how the amount of allowed concavity is controlled.
#' @export
fm_nonconvex_hull.sfc <- function(x,
                                  convex = -0.15,
                                  concave = convex,
                                  preserveTopology = TRUE,
                                  dTolerance = NULL,
                                  crs = fm_crs(x),
                                  ...) {
  if ((convex < 0) || (concave < 0)) {
    diameter_bound <- fm_diameter(x)
    if (convex < 0) {
      convex <- diameter_bound * abs(convex)
    }
    if (concave < 0) {
      concave <- diameter_bound * abs(concave)
    }
  }

  nQuadSegs <- 64
  y <- sf::st_buffer(x, dist = convex + concave, nQuadSegs = nQuadSegs)
  y <- sf::st_union(y)
  if (concave > 0) {
    y <- sf::st_buffer(y, dist = -concave, nQuadSegs = nQuadSegs)
  }
  if (!is.null(dTolerance)) {
    y <- sf::st_simplify(y,
      preserveTopology = preserveTopology,
      dTolerance = dTolerance
    )
  }
  y <- sf::st_union(y)
  y <- sf::st_sfc(y, crs = fm_crs(x))
  if (!fm_crs_is_identical(fm_crs(y), crs)) {
    y <- fm_transform(y, crs = crs)
  }
  y
}


#' @describeIn fm_nonconvex_hull
#' Constructs a potentially nonconvex extension of a spatial object by
#' performing dilation by `convex + concave` followed by
#' erosion by `concave`. This is equivalent to dilation by `convex` followed
#' by closing (dilation + erosion) by `concave`.
#'
#' @param x A spatial object
#' @param convex numeric vector; How much to extend
#' @param concave numeric vector; The minimum allowed reentrant curvature. Default equal to `convex`
#' @param dTolerance If not null, the `dTolerance` argument to `sf::st_simplify()`,
#' passed on to [fm_nonconvex_hull()].
#' The default for `fm_extensions()` is `pmin(convex, concave) / 40`, chosen to
#' give approximately 4 or more subsegments per circular quadrant.
#' For `fm_nonconvex_hull()`, the default is `NULL`
#' (i.e. not to call `sf::st_simplify()`)
#' @returns `fm_extensions()` returns a list of `sfc` objects.
#' @export
#' @examples
#' if (TRUE) {
#'   inp <- sf::st_as_sf(as.data.frame(matrix(1:6, 3, 2)), coords = 1:2)
#'   bnd <- fm_extensions(inp, convex = c(0.75, 2))
#'   plot(fm_mesh_2d(boundary = bnd, max.edge = c(0.25, 1)), asp = 1)
#' }
fm_extensions <- function(x,
                          convex = -0.15,
                          concave = convex,
                          dTolerance = NULL,
                          ...) {
  if (any(convex < 0) || any(concave < 0) || any(dTolerance < 0)) {
    diameter_bound <- fm_diameter(x)
  }
  len <- max(length(convex), length(concave), length(dTolerance))
  scale_fun <- function(val) {
    if (any(val < 0)) {
      val[val < 0] <- diameter_bound * abs(val[val < 0])
    }
    if (length(val) < len) {
      val <- c(val, rep(val[length(val)], len - length(val)))
    }
    val
  }
  convex <- scale_fun(convex)
  concave <- scale_fun(concave)
  if (is.null(dTolerance)) {
    dTolerance <- pmin(convex, concave) / 40
  } else {
    dTolerance <- scale_fun(dTolerance)
  }

  y <- lapply(
    seq_along(convex),
    function(k) {
      fm_nonconvex_hull(
        x,
        convex = convex[k],
        concave = concave[k],
        dTolerance = dTolerance[k],
        ...
      )
    }
  )
  y
}



#' @rdname fm_nonconvex_hull
#' @export
fm_nonconvex_hull.matrix <- function(x, ...) {
  fm_nonconvex_hull.sfc(sf::st_multipoint(x), ...)
}

#' @rdname fm_nonconvex_hull
#' @export
fm_nonconvex_hull.sf <- function(x, ...) {
  fm_nonconvex_hull.sfc(sf::st_geometry(x), ...)
}

#' @rdname fm_nonconvex_hull
#' @export
fm_nonconvex_hull.Spatial <- function(x, ...) {
  fm_nonconvex_hull.sfc(sf::st_as_sfc(x), ...)
}

#' @rdname fm_nonconvex_hull
#' @export
fm_nonconvex_hull.sfg <- function(x, ...) {
  fm_nonconvex_hull.sfc(sf::st_sfc(x), ...)
}


# fm_diameter ####

#' Diameter bound for a geometric object
#'
#' Find an upper bound to the convex hull of a point set
#'
#' @param x A point set as an \eqn{n\times d}{n x d} matrix, or an
#' `fm_mesh_2d`/`1d`/`sf` related object.
#' @param manifold Character string specifying the manifold type. Default is to
#' treat the point set with Euclidean \eqn{R^d} metrics. Use
#' `manifold="S2"` for great circle distances on the unit sphere (this is
#' set automatically for `fm_fmesh_2d` objects).
#' @param \dots Additional parameters passed on to the submethods.
#' @return A scalar, upper bound for the diameter of the convex hull of the
#' point set.
#' @author Finn Lindgren <finn.lindgren@@gmail.com>
#' @examples
#'
#' fm_diameter(matrix(c(0, 1, 1, 0, 0, 0, 1, 1), 4, 2))
#' @export
fm_diameter <- function(x, ...) {
  if (is.null(x)) {
    return(0.0)
  }
  UseMethod("fm_diameter")
}

#' @rdname fm_diameter
#' @export
fm_diameter.matrix <- function(x, manifold = NULL, ...) {
  if (nrow(x) <= 1) {
    return(0)
  }
  if (ncol(x) == 1) {
    return(diff(range(x)))
  }

  if (identical(manifold, "S2")) {
    radius <- mean(rowSums(x^2)^0.5)
    x <- x / radius
    distance <- function(u, v) {
      2 * asin(pmin(
        1,
        ((u[1] - v[, 1])^2 + (u[2] - v[, 2])^2 + (u[3] - v[, 3])^2)^0.5 / 2
      ))
    }
    center <- colMeans(x)
    tmp <- sqrt(sum(center^2))
    if (tmp < 1e-6 * radius) {
      diam <- pi
    } else {
      center <- center / tmp
      diam <- min(pi, 2 * max(distance(center, x)))
    }
    diam <- diam * radius
    return(diam)
  }

  distance <- function(u, v) {
    d <- 0
    for (k in seq_len(ncol(v))) {
      d <- d + (u[k] - v[, k])^2
    }
    d^0.5
  }
  center <- rep(0, ncol(x))
  for (k in seq_len(ncol(x))) {
    center[k] <- mean(range(x[, k]))
  }
  diam <- 2 * max(distance(center, x))
  return(diam)
}

#' @rdname fm_diameter
#' @export
fm_diameter.sf <- function(x, ...) {
  fm_diameter.sfc(sf::st_geometry(x))
}

#' @rdname fm_diameter
#' @export
fm_diameter.sfg <- function(x, ...) {
  fm_diameter.sfc(sf::st_sfc(x))
}

#' @rdname fm_diameter
#' @export
fm_diameter.sfc <- function(x, ...) {
  fm_diameter.matrix(sf::st_coordinates(x))
}

#' @rdname fm_diameter
#' @export
fm_diameter.fm_lattice_2d <- function(x, ...) {
  fm_diameter.matrix(x$loc, manifold = fm_manifold(x), ...)
}

#' @rdname fm_diameter
#' @export
fm_diameter.fm_segm <- function(x, ...) {
  fm_diameter.matrix(x$loc, manifold = fm_manifold(x), ...)
}

#' @rdname fm_diameter
#' @export
fm_diameter.fm_mesh_2d <- function(x, ...) {
  fm_diameter.matrix(x$loc, manifold = fm_manifold(x), ...)
}

#' @rdname fm_diameter
#' @export
fm_diameter.fm_mesh_1d <- function(x, ...) {
  diff(x[["interval"]])
}

#' @rdname fm_diameter
#' @export
#' @method fm_diameter inla.mesh.1d
fm_diameter.inla.mesh.1d <- function(x, ...) {
  fm_diameter(fm_as_fm(x), ...)
}

#' @rdname fm_diameter
#' @export
#' @method fm_diameter inla.mesh.segment
fm_diameter.inla.mesh.segment <- function(x, ...) {
  fm_diameter(fm_as_fm(x), ...)
}

#' @rdname fm_diameter
#' @export
#' @method fm_diameter inla.mesh.lattice
fm_diameter.inla.mesh.lattice <- function(x, ...) {
  fm_diameter(fm_as_fm(x), ...)
}

#' @rdname fm_diameter
#' @export
#' @method fm_diameter inla.mesh
fm_diameter.inla.mesh <- function(x, ...) {
  fm_diameter(fm_as_fm(x), ...)
}








# fm_diameter ####

#' @title Query the mesh manifold type
#' @description
#' Extract a manifold definition string, or a logical for matching
#' manifold type
#' @param x A [fm_mesh_1d] or [fm_mesh_2d] object (or other object containing a
#' `manifold` element)
#' @param type `character`; if `NULL` (the default), returns the manifold definition string.
#' If `character`, returns `TRUE` if the manifold type of `x` matches at least
#' one of the character vector elements.
#' @export
fm_manifold <- function(x, type = NULL) {
  if (is.null(type)) {
    return(x[["manifold"]])
  }
  if (is.null(x) || is.null(x[["manifold"]])) {
    return(FALSE)
  }
  # Match exact manifold?
  if (x[["manifold"]] %in% type) {
    return(TRUE)
  }
  # Match space name or dimension?
  m <- intersect(c("M", "R", "S"), type)
  d <- intersect(as.character(seq_len(3)), type)
  if (length(c(m, d)) == 0) {
    return(FALSE)
  }
  grepl(paste0(c(m, d), collapse = "|"), x[["manifold"]])
}

#' @rdname fm_manifold
#' @export
fm_manifold_type <- function(x) {
  if (is.null(x) || is.null(x[["manifold"]])) {
    return(NULL)
  }
  # Match space name or dimension?
  m <- intersect(c("M", "R", "S"), strsplit(x[["manifold"]], "")[[1]])
  if (length(m) == 0) {
    return(NULL)
  }
  m
}

#' @rdname fm_manifold
#' @export
fm_manifold_dim <- function(x) {
  if (is.null(x) || is.null(x[["manifold"]])) {
    return(NULL)
  }
  # Match space name or dimension?
  d <- intersect(as.character(seq_len(3)), strsplit(x[["manifold"]], "")[[1]])
  if (length(d) == 0) {
    return(NULL)
  }
  as.integer(d)
}


# fm_segm ####

#' @title Make a spatial segment object
#' @describeIn fm_segm Create a new `fm_segm` object.
#' @export
#' @param ... Passed on to submethods
#' @family object creation and conversion
fm_segm <- function(...) {
  UseMethod("fm_segm")
}

#' @rdname fm_segm
#' @param loc Matrix of point locations, or `SpatialPoints`, or `sf`/`sfc` point
#' object.
#' @param idx Segment index sequence vector or index pair matrix.  The indices
#' refer to the rows of `loc`.  If `loc==NULL`, the indices will be
#' interpreted as indices into the point specification supplied to
#' [fm_rcdt_2d()].  If `is.bnd==TRUE`, defaults to linking
#' all the points in `loc`, as `c(1:nrow(loc),1L)`, otherwise
#' `1:nrow(loc)`.
#' @param grp Vector of group labels for each segment.  Set to `NULL` to
#' let the labels be chosen automatically in a call to
#' [fm_rcdt_2d()].
#' @param is.bnd `TRUE` if the segments are boundary segments, otherwise
#' `FALSE`.
#' @param crs An optional `fm_crs()`, `sf::st_crs()` or `sp::CRS()` object
#' @export
fm_segm.default <- function(loc = NULL, idx = NULL, grp = NULL, is.bnd = TRUE,
                            crs = NULL, ...) {
  if (is.null(loc) && is.null(idx)) {
    idx <- matrix(0L, 0, 2)
  }
  if (!is.null(loc)) {
    loc <- unify_loc_coords(loc)
    if (is.null(idx)) {
      idx <- if (is.bnd) {
        c(seq_len(nrow(loc)), 1)
      } else {
        seq_len(nrow(loc))
      }
    }
  }

  if (!is.null(idx)) {
    if (!is.vector(idx) && !is.matrix(idx)) {
      stop("'idx' must be a vector or a matrix")
    }
    if (is.vector(idx)) {
      idx <- as.matrix(idx, nrow = length(idx), ncol = 1)
    }
    if (ncol(idx) == 1) {
      if (nrow(idx) < 2) {
        if (nrow(idx) == 1) {
          warning("Segment specification must have at least 2, or 0, indices.")
        }
        idx <- matrix(0L, 0, 2)
      } else {
        idx <- matrix(c(idx[-nrow(idx)], idx[-1]), nrow(idx) - 1, 2)
      }
    }
    storage.mode(idx) <- "integer"
    if (!is.null(loc) &&
      (nrow(idx) > 0) &&
      (max(idx, na.rm = TRUE) > nrow(loc))) {
      warning(
        "Segment indices (max=", max(idx, na.rm = TRUE),
        ") exceed specified location list length (",
        nrow(loc), ")."
      )
    }
  }

  if (!is.null(grp)) {
    if (!is.vector(grp) && !is.matrix(grp)) {
      stop("'grp' must be a vector or a matrix")
    }
    grp <- as.vector(grp)
    if (length(grp) < nrow(idx)) {
      grp <- c(
        grp,
        rep(grp[length(grp)], nrow(idx) - length(grp))
      )
    }
    storage.mode(grp) <- "integer"
  }

  ## Filter away NAs in loc and idx
  if (!is.null(loc)) {
    idx[is.na(idx)] <- 0L ## Avoid R annoyances with logical+NA indexing
    while (sum(is.na(loc)) > 0) {
      i <- min(which(rowSums(is.na(loc)) > 0))
      loc <- loc[-i, , drop = FALSE]
      idx[idx == i] <- 0L
      idx[idx > i] <- idx[idx > i] - 1L
    }
    idx[idx == 0L] <- NA
  }
  while (sum(is.na(idx)) > 0) {
    i <- min(which(rowSums(is.na(idx)) > 0))
    idx <- idx[-i, , drop = FALSE]
    if (!is.null(grp)) {
      grp <- grp[-i]
    }
  }

  if (!is.null(loc)) {
    ## Identify unused locations and remap indices accordingly.
    idx.new <- rep(0L, nrow(loc))
    idx.new[as.vector(idx)] <- 1L
    loc <- loc[idx.new == 1L, , drop = FALSE]
    idx.new[idx.new == 1L] <- seq_len(sum(idx.new))
    idx <- matrix(idx.new[as.vector(idx)],
      nrow = nrow(idx),
      ncol = ncol(idx)
    )
  }

  ret <- list(loc = loc, idx = idx, grp = grp, is.bnd = is.bnd, crs = crs)
  class(ret) <- c("fm_segm", "inla.mesh.segment")
  return(ret)
}

#' @describeIn fm_segm Join multiple `fm_segm` objects into a single `fm_segm`
#' object.
#' @param grp When joining segments, use these group labels for segments
#' instead of the original group labels.
#' @param grp.default If `grp.default` is `NULL`, use these group labels for segments
#' with NULL group.
#' @export
fm_segm.fm_segm <- function(..., grp = NULL, grp.default = 0L) {
  segm <- fm_as_segm_list(list(...))
  fm_segm_join(segm, grp = grp, grp.default = grp.default)
}
#' @describeIn fm_segm Join `fm_segm` objects from a `fm_segm_list` into
#' a single `fm_segm` object. Equivalent to `fm_segm_join(x)`
#' @param grp When joining segments, use these group labels for segments
#' instead of the original group labels.
#' @param grp.default If `grp.default` is `NULL`, use these group labels for segments
#' with NULL group.
#' @param x A `fm_segm_list` object
#' @export
fm_segm.fm_segm_list <- function(x, grp = NULL, grp.default = 0L, ...) {
  fm_segm_join(x, grp = grp, grp.default = grp.default)
}
#' @describeIn fm_segm Join multiple `fm_segm` objects into a single `fm_segm`
#' object.
#' @export
fm_segm_join <- function(x, grp = NULL, grp.default = 0L) {
  segm <- fm_as_segm_list(x)
  segm <- lapply(seq_along(segm), function(k) {
    seg <- segm[[k]]
    if (!is.null(seg)) {
      if (is.null(grp)) {
        if (is.null(seg[["grp"]])) {
          seg[["grp"]] <-
            rep(
              grp.default[min(length(grp.default), k)],
              nrow(seg[["idx"]])
            )
        }
      } else {
        seg[["grp"]] <-
          rep(grp[min(length(grp), k)], nrow(seg[["idx"]]))
      }
    }
    seg
  })

  keep <- vapply(segm, function(x) !is.null(x), TRUE)
  segm <- segm[keep]
  if (!all(vapply(
    segm,
    function(x) inherits(x, "fm_segm"),
    TRUE
  ))) {
    stop("All objects must be of class 'fm_segm'.")
  }

  Nloc <- vapply(segm, function(x) nrow(x$loc), 0L)
  cumNloc <- c(0, cumsum(Nloc))
  Nidx <- vapply(segm, function(x) nrow(x$idx), 0L)

  loc <- do.call(rbind, lapply(segm, function(x) x$loc))
  idx <- do.call(rbind, lapply(
    seq_along(segm),
    function(k) segm[[k]]$idx + cumNloc[k]
  ))
  grp <- unlist(lapply(
    seq_along(segm),
    function(k) {
      if (is.null(segm[[k]]$grp)) {
        rep(grp.default[min(length(grp.default), k)], Nidx[k])
      } else {
        segm[[k]]$grp
      }
    }
  ))
  is.bnd <- vapply(segm, function(x) fm_is_bnd(x), TRUE)
  if (all(is.bnd) || all(!is.bnd)) {
    is.bnd <- all(is.bnd)
  } else {
    warning("Inconsistent 'is.bnd' attributes.  Setting 'is.bnd=FALSE'.")
    is.bnd <- FALSE
  }

  crs <- lapply(segm, function(x) fm_crs(x))
  if (!is.null(crs)) {
    crs <- crs[vapply(
      crs,
      function(x) !fm_crs_is_null(x),
      TRUE
    )]
    if (length(crs) > 0) {
      if (!all(vapply(
        crs,
        function(x) fm_crs_is_identical(crs[[1]], x),
        TRUE
      ))) {
        lapply(crs, function(x) show(x))
        stop("Inconsistent 'crs' attributes.")
      } else {
        crs <- crs[[1]]
      }
    } else {
      crs <- NULL
    }
  }

  fm_segm(
    loc = loc,
    idx = idx,
    grp = grp,
    is.bnd = is.bnd,
    crs = crs
  )
}
#' @describeIn fm_segm Split an `fm_segm` object by `grp` into a list of `fm_segm`
#' objects, optionally keeping only some groups.
#' @export
fm_segm_split <- function(x, grp = NULL, grp.default = 0L) {
  if (is.null(x[["grp"]])) {
    x[["grp"]] <- rep(grp.default, nrow(x[["idx"]]))
  }
  if (is.null(grp)) {
    grp <- sort(unique(x[["grp"]]))
  }
  segm_list <- lapply(
    grp,
    function(g) {
      keep <- x[["grp"]] == g
      fm_segm(
        loc = x[["loc"]],
        idx = x[["idx"]][keep, , drop = FALSE],
        grp = rep(g, sum(keep)),
        is.bnd = x[["is.bnd"]],
        crs = fm_crs(x)
      )
    }
  )
  return(segm_list)
}
#' @rdname fm_segm
#' @export
#' @method fm_segm inla.mesh.segment
fm_segm.inla.mesh.segment <- function(..., grp.default = 0) {
  do.call(
    fm_segm,
    c(
      fm_as_segm_list(list(...)),
      list(grp.default = grp.default)
    )
  )
}
#' @rdname fm_segm
#' @export
#' @method fm_segm inla.mesh
fm_segm.inla.mesh <- function(x, ...) {
  fm_segm.fm_mesh_2d(fm_as_mesh_2d(x), ...)
}
#' @describeIn fm_segm Extract the boundary or interior segments of a 2d mesh.
#' If `grp` is non-NULL, extracts only segments matching the matching the set
#' of groups given by `grp`.
#' @param x Mesh to extract segments from
#' @param boundary logical; if `TRUE`, extract the boundary segments,
#' otherwise interior constrain segments.
#' @export
fm_segm.fm_mesh_2d <- function(x, boundary = TRUE, grp = NULL, ...) {
  extract_segments <- function(mesh.loc,
                               segm,
                               grp = NULL,
                               is.bnd,
                               crs = NULL) {
    segments <- NULL
    if (nrow(segm[["idx"]]) == 0) {
      return(NULL)
    }
    if (is.null(grp)) {
      grp <- unique(sort(segm[["grp"]]))
    }
    extract <- (segm[["grp"]] %in% grp)
    if (!any(extract)) {
      return(NULL)
    }
    segments <-
      fm_segm(
        mesh.loc,
        idx = segm[["idx"]][extract, , drop = FALSE],
        grp = segm[["grp"]][extract, drop = FALSE],
        is.bnd = is.bnd,
        crs = crs
      )
    return(segments)
  }

  if (boundary) {
    segm <- x[["segm"]][["bnd"]]
  } else {
    segm <- x[["segm"]][["int"]]
  }

  extract_segments(
    mesh.loc = x[["loc"]],
    segm = segm,
    grp = NULL,
    is.bnd = boundary,
    crs = fm_crs(x)
  )
}


#' @title Convert objects to `fm_segm`
#' @describeIn fm_as_segm Convert an object to `fm_segm`.
#' @param x Object to be converted.
#' @param ... Arguments passed on to submethods
#' @export
#' @family object creation and conversion
fm_as_segm <- function(x, ...) {
  if (is.null(x)) {
    return(NULL)
  }
  UseMethod("fm_as_segm")
}

#' @describeIn fm_as_segm Convert each element, making a `fm_segm_list` object
#' @seealso [c.fm_segm()], [c.fm_segm_list()],
#' \code{\link[=[.fm_segm_list]{[.fm_segm_list()}}
#' @export
fm_as_segm_list <- function(x, ...) {
  if (inherits(x, c("fm_segm_list", "NULL"))) {
    return(x)
  }
  structure(
    fm_as_list(x, ..., .method = "fm_as_segm"),
    class = "fm_segm_list"
  )
}

#' @rdname fm_as_segm
#' @export
fm_as_segm.fm_segm <- function(x, ...) {
  x
}

#' @rdname fm_as_segm
#' @export
#' @method fm_as_segm inla.mesh.segment
fm_as_segm.inla.mesh.segment <- function(x, ...) {
  class(x) <- c("fm_segm", class(x))
  x
}

#' @rdname fm_segm
#' @export
`fm_is_bnd` <- function(x) {
  x[["is.bnd"]]
}

#' @rdname fm_segm
#' @param value logical
#' @export
`fm_is_bnd<-` <- function(x, value) {
  if (is.null(value)) {
    value <- FALSE
  }
  value <- as.logical(value)
  x[["is.bnd"]] <- rep(value, length(x[["is.bnd"]]))
  invisible(x)
}




#' Methods for fm_segm lists
#'
#' `fm_segm` lists can be combined into `fm_segm_list` list objects.
#' @name fm_segm_list
#' @param \dots Objects to be combined.
#' @seealso [fm_as_segm_list()]
NULL

#' @export
#' @describeIn fm_segm_list The `...` arguments should be `fm_segm`
#' objects, or coercible with `fm_as_segm_list(list(...))`.
#' @returns A `fm_segm_list` object
#' @examples
#' m <- c(A = fm_segm(1:2), B = fm_segm(3:4))
#' str(m)
#' str(m[2])
`c.fm_segm` <- function(...) {
  fm_as_segm_list(list(...))
}

#' @export
#' @describeIn fm_segm_list The `...` arguments should be `fm_segm_list`
#' objects.
`c.fm_segm_list` <- function(...) {
  y <- lapply(list(...), fm_as_segm_list)
  object <- do.call(NextMethod, list("c", y))
  class(object) <- "fm_segm_list"
  object
}

#' @export
#' @param x `fm_segm_list` object from which to extract element(s)
#' @param i indices specifying elements to extract
#' @describeIn fm_segm_list Extract sub-list
`[.fm_segm_list` <- function(x, i) {
  object <- NextMethod()
  class(object) <- "fm_segm_list"
  object
}





# fm_mesh ####

#' @title Convert objects to fmesher objects
#' @description
#' Used for conversion from general objects
#' (usually `inla.mesh` and other INLA specific classes)
#' to `fmesher` classes.
#'
#' @param x Object to be converted
#' @param ... Arguments forwarded to submethods
#' @rdname fm_as_fm
#' @export
#' @family object creation and conversion
fm_as_fm <- function(x, ...) {
  UseMethod("fm_as_fm")
}

# @description fm_as_list Convert each element of a list, or convert a single
# object and return in a list
# @param x list of objects to be converted.
# @param ... Arguments passed to each individual conversion call.
# @param .method character; name of a conversion generic to apply to each list
# element.
fm_as_list <- function(x, ..., .method) {
  if (is.null(x)) {
    return(list())
  }
  m_c <- method_classes(.method)
  if (inherits(x, setdiff(m_c, "list"))) {
    return(list(do.call(.method, list(x, ...))))
  }
  if (!inherits(x, "list")) {
    stop(paste0(
      "'list' object expected. Received '",
      paste0(class(x), collapse = ", "),
      "'."
    ))
  }
  if ("list" %in% m_c) {
    return(do.call(.method, list(x, ...)))
  }
  lapply(x, function(xx) do.call(.method, list(xx, ...)))
}
#' @describeIn fm_as_fm Convert each element of a list, or convert a single
#' object and return in a list
#' @export
fm_as_fm_list <- function(x, ...) {
  y <- fm_as_list(x, ..., .method = "fm_as_fm")
  if (all(vapply(y, function(xx) inherits(xx, "fm_segm"), TRUE))) {
    y <- fm_as_segm_list(y)
  }
  y
}
#' @rdname fm_as_fm
#' @usage
#' ## S3 method for class 'NULL'
#' fm_as_fm(x, ...)
#' @export
fm_as_fm.NULL <- function(x, ...) {
  NULL
}
#' @rdname fm_as_fm
#' @export
fm_as_fm.fm_mesh_1d <- function(x, ...) {
  #  class(x) <- c("fm_mesh_1d", setdiff(class(x), "fm_mesh_1d"))
  x
}
#' @rdname fm_as_fm
#' @export
fm_as_fm.fm_mesh_2d <- function(x, ...) {
  #  class(x) <- c("fm_mesh_2d", setdiff(class(x), "fm_mesh_2d"))
  x
}
#' @rdname fm_as_fm
#' @export
fm_as_fm.fm_segm <- function(x, ...) {
  #  class(x) <- c("fm_segm", setdiff(class(x), "fm_segm"))
  x
}
#' @rdname fm_as_fm
#' @export
fm_as_fm.fm_lattice_2d <- function(x, ...) {
  #  class(x) <- c("fm_lattice_2d", setdiff(class(x), "fm_lattice_2d"))
  x
}
#' @rdname fm_as_fm
#' @export
fm_as_fm.crs <- function(x, ...) {
  fm_crs(x)
}
#' @rdname fm_as_fm
#' @export
fm_as_fm.CRS <- function(x, ...) {
  fm_crs(x)
}
#' @rdname fm_as_fm
#' @export
fm_as_fm.fm_crs <- function(x, ...) {
  fm_crs(x)
}
#' @rdname fm_as_fm
#' @export
#' @method fm_as_fm inla.CRS
fm_as_fm.inla.CRS <- function(x, ...) {
  fm_crs(x)
}
#' @rdname fm_as_fm
#' @export
#' @method fm_as_fm inla.mesh.1d
fm_as_fm.inla.mesh.1d <- function(x, ...) {
  fm_as_mesh_1d(x, ...)
}
#' @rdname fm_as_fm
#' @export
#' @method fm_as_fm inla.mesh
fm_as_fm.inla.mesh <- function(x, ...) {
  fm_as_mesh_2d(x, ...)
}

#' @rdname fm_as_fm
#' @export
#' @method fm_as_fm inla.mesh.segment
fm_as_fm.inla.mesh.segment <- function(x, ...) {
  fm_as_segm(x, ...)
}

#' @rdname fm_as_fm
#' @export
#' @method fm_as_fm inla.mesh.lattice
fm_as_fm.inla.mesh.lattice <- function(x, ...) {
  fm_as_lattice_2d(x, ...)
}

# fm_mesh_1d ####

`match.arg.vector` <- function(arg = NULL,
                               choices,
                               length = NULL) {
  ## Like match.arg, but for a vector of options 'arg'
  if (is.null(length)) {
    length <- ifelse(is.null(arg), 1L, length(arg))
  }
  if (is.null(arg)) {
    arg <- match.arg(arg, choices)
  } else {
    for (k in seq_along(arg)) {
      arg[k] <- match.arg(arg[k], choices)
    }
  }
  if (length(arg) < length) {
    arg <- c(arg, rep(arg, length - length(arg)))
  } else if (length(arg) > length) {
    stop("Option list too long.")
  }
  return(arg)
}

#' @title Make a 1D mesh object
#' @description
#' Create a `fm_mesh_1d` object.
#'
#' @param loc B-spline knot locations.
#' @param interval Interval domain endpoints.
#' @param boundary Boundary condition specification.  Valid conditions are
#' `c('neumann', 'dirichlet', 'free', 'cyclic')`.  Two separate values can
#' be specified, one applied to each endpoint.
#' @param degree The B-spline basis degree.  Supported values are 0, 1, and 2.
#' @param free.clamped If `TRUE`, for `'free'` boundaries, clamp the
#' basis functions to the interval endpoints.
#' @param \dots Additional options, currently unused.
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @export
#' @family object creation and conversion
fm_mesh_1d <- function(loc,
                       interval = range(loc),
                       boundary = NULL,
                       degree = 1,
                       free.clamped = FALSE,
                       ...) {
  ## Note: do not change the order of these options without also
  ## changing 'basis.reduction' below.
  boundary.options <- c("neumann", "dirichlet", "free", "cyclic")

  boundary <- match.arg.vector(boundary, boundary.options, length = 2)
  cyclic <- !is.na(pmatch(boundary[1], "cyclic"))
  if (cyclic && is.na(pmatch(boundary[2], "cyclic"))) {
    stop("Inconsistent boundary specification 'boundary=c(",
      paste(boundary, collapse = ","), ")'.",
      sep = ""
    )
  }

  loc.orig <- loc
  if (cyclic) {
    loc <-
      (sort(unique(c(0, loc - interval[1]) %% diff(interval))) +
        interval[1])
  } else {
    loc <-
      (sort(unique(c(
        interval,
        pmax(interval[1], pmin(interval[2], loc))
      ))))
  }

  n <- length(loc)

  if (loc[1] < interval[1]) {
    stop("All 'loc' must be >= interval[1].")
  }
  if (loc[n] > interval[2]) {
    stop("All 'loc' must be <= interval[2].")
  }

  if ((degree < 0) || (degree > 2)) {
    stop(paste("'degree' must be 0, 1, or 2.  'degree=",
      degree,
      "' is not supported.",
      sep = ""
    ))
  }

  if (length(free.clamped) == 1L) {
    free.clamped <- rep(free.clamped, 2)
  }


  ## Number of basis functions
  if (degree == 0) {
    basis.reduction <- c(0, 1, 0, 1 / 2) ## neu, dir, free, cyclic
  } else if (degree == 1) {
    basis.reduction <- c(0, 1, 0, 1 / 2) ## neu, dir, free, cyclic
  } else {
    basis.reduction <- c(1, 1, 0, 1)
  }
  m <- (n + cyclic + (degree == 2) * 1
    - basis.reduction[pmatch(boundary[1], boundary.options)]
    - basis.reduction[pmatch(boundary[2], boundary.options)])
  ## if (m < 1+max(1,degree)) {
  if (m < 1L) {
    stop("Degree ", degree,
      " meshes must have at least ", 1L,
      " basis functions, not 'm=", m, "'.",
      sep = ""
    )
  }

  ## Compute representative basis midpoints.
  if ((degree == 0) || (degree == 1)) {
    mid <- loc
    if (boundary[1] == "dirichlet") {
      mid <- mid[-1]
    }
    if (boundary[2] == "dirichlet") {
      mid <- mid[-(m + 1)]
    }
  } else { ## degree==2
    if (cyclic) {
      mid <- (loc + c(loc[-1], interval[2])) / 2
    } else {
      mid <- c(loc[1], (loc[-n] + loc[-1]) / 2, loc[n])
      mid <-
        switch(boundary[1],
          neumann = mid[-1],
          dirichlet = mid[-1],
          free = mid
        )
      mid <-
        switch(boundary[2],
          neumann = mid[-(m + 1)],
          dirichlet = mid[-(m + 1)],
          free = mid
        )
    }
  }

  mesh <-
    structure(
      list(
        n = n,
        m = m,
        loc = loc,
        mid = mid,
        interval = interval,
        boundary = boundary,
        cyclic = cyclic,
        manifold = ifelse(cyclic, "S1", "R1"),
        degree = degree,
        free.clamped = free.clamped,
        idx = list(loc = NULL)
      ),
      class = c("fm_mesh_1d", "inla.mesh.1d")
    )

  if (degree < 2) {
    mesh$idx$loc <-
      fm_bary(mesh, loc.orig, method = "nearest")$index[, 1]
  } else {
    if (length(mid) >= 2) {
      mesh$idx$loc <-
        fm_bary(fm_mesh_1d(mid, degree = 0),
          loc.orig,
          method = "nearest"
        )$index[, 1]
    } else {
      mesh$idx$loc <- rep(1, length(loc.orig))
    }
  }

  return(mesh)
}



#' @title Convert objects to `fm_segm`
#' @describeIn fm_as_mesh_1d Convert an object to `fm_mesh_1d`.
#' @param x Object to be converted.
#' @param ... Arguments passed on to submethods
#' @export
#' @family object creation and conversion
#' @export
fm_as_mesh_1d <- function(x, ...) {
  if (is.null(x)) {
    return(NULL)
  }
  UseMethod("fm_as_mesh_1d")
}
#' @describeIn fm_as_mesh_1d Convert each element of a list
#' @export
fm_as_mesh_1d_list <- function(x, ...) {
  fm_as_list(x, ..., .method = "fm_as_mesh_1d")
}
#' @rdname fm_as_mesh_1d
#' @param x Object to be converted
#' @export
fm_as_mesh_1d.fm_mesh_1d <- function(x, ...) {
  #  class(x) <- c("fm_mesh_1d", setdiff(class(x), "fm_mesh_1d"))
  x
}
#' @rdname fm_as_mesh_1d
#' @param x Object to be converted
#' @export
#' @method fm_as_mesh_1d inla.mesh.1d
fm_as_mesh_1d.inla.mesh.1d <- function(x, ...) {
  class(x) <- c("fm_mesh_1d", class(x))
  x
}


# fm_mesh_2d ####

unify_loc_coords <- function(loc, crs = crs) {
  if (is.null(loc) || (NROW(loc) == 0)) {
    return(matrix(c(0.0), 0, 3))
  }
  if (inherits(loc, c(
    "SpatialPoints",
    "SpatialPointsDataFrame"
  ))) {
    loc <- fm_transform(
      sp::coordinates(loc),
      crs0 = fm_crs(loc),
      crs = crs,
      passthrough = TRUE
    )
  } else if (inherits(loc, c("sf", "sfg", "sfc"))) {
    loc <- fm_transform(
      sf::st_coordinates(loc),
      crs0 = fm_crs(loc),
      crs = crs,
      passthrough = TRUE
    )
  }
  if (!is.matrix(loc)) {
    if (is.vector(loc)) {
      loc <- matrix(loc, 1, length(loc))
    } else {
      loc <- as.matrix(loc)
    }
  }
  if (ncol(loc) < 3) {
    while (ncol(loc) < 3) {
      loc <- cbind(loc, 0.0)
    }
  } else if (ncol(loc) > 3) {
    stop("Coordinates can have at mots 3 columns.")
  }
  loc
}

unify_segm_coords <- function(segm, crs = NULL) {
  if (is.null(segm)) {
    return(NULL)
  }
  unify.one.segm <- function(segm, crs = NULL) {
    segm <- fm_transform(segm, crs, passthrough = TRUE)
    if (ncol(segm$loc) == 2) {
      segm$loc <- cbind(segm$loc, 0.0)
    }
    segm
  }
  if (inherits(segm, "fm_segm")) {
    segm <- unify.one.segm(segm, crs = crs)
  } else {
    for (j in seq_along(segm)) {
      if (!is.null(segm[[j]])) {
        segm[[j]] <- unify.one.segm(segm[[j]], crs = crs)
      }
    }
  }
  segm
}




handle_rcdt_options_inla <- function(
    ...,
    quality.spec = NULL,
    cutoff = 1e-12,
    extend = NULL,
    refine = NULL,
    .n,
    .loc) {
  options <- list(cutoff = cutoff)
  if (is.null(quality.spec)) {
    quality <- NULL
  } else {
    quality <- rep(NA, .n$segm + .n$lattice + .n$loc)
    ## Order must be same as in .loc
    if (!is.null(quality.spec$segm)) {
      quality[seq_len(.n$segm)] <- quality.spec$segm
    }
    if (!is.null(quality.spec$lattice)) {
      quality[.n$segm + seq_len(.n$lattice)] <- quality.spec$lattice
    }
    if (!is.null(quality.spec$loc)) {
      quality[.n$segm + .n$lattice + seq_len(.n$loc)] <- quality.spec$loc
    }
    ## NA:s will be replaced with max.edge settings below.

    options <- c(options, list(quality = quality))
  }

  cet_sides <- NULL
  cet_margin <- NULL
  if (inherits(extend, "list")) {
    cet_sides <- ifelse(is.null(extend$n), 16L, as.integer(extend$n))
    cet_margin <- ifelse(is.null(extend$offset), -0.1, extend$offset)
  }
  options <- c(options, list(cet_sides = cet_sides, cet_margin = cet_margin))

  if (inherits(refine, "list")) {
    rcdt_min_angle <- 0
    rcdt_max_edge <- 0
    # Multiply by 2 to cover S2; could remove if supplied manifold info
    max.edge.default <- fm_diameter(.loc) * 2
    if ((inherits(extend, "list")) && (!is.null(extend$offset))) {
      max.edge.default <- (max.edge.default +
        max(0, 2 * extend$offset))
      max.edge.default <- (max.edge.default *
        (1 + max(0, -2 * extend$offset)))
    }
    rcdt_min_angle <- ifelse(is.null(refine$min.angle), 21, refine$min.angle)
    rcdt_max_edge <- ifelse(
      is.null(refine$max.edge) ||
        is.na(refine$max.edge),
      max.edge.default,
      refine$max.edge
    )
    max_edge_extra <- ifelse(
      is.null(refine$max.edge.extra) ||
        is.na(refine$max.edge.extra),
      rcdt_max_edge,
      refine$max.edge.extra
    )

    if (!is.null(refine[["max.n.strict"]]) &&
      !is.na(refine$max.n.strict)) {
      rcdt_max_n0 <- as.integer(refine$max.n.strict)
    } else {
      rcdt_max_n0 <- -1L
    }
    if (!is.null(refine[["max.n"]]) &&
      !is.na(refine$max.n)) {
      rcdt_max_n1 <- as.integer(refine$max.n)
    } else {
      rcdt_max_n1 <- -1L
    }

    if (!is.null(options[["quality"]])) {
      options[["quality"]][is.na(options[["quality"]])] <- max_edge_extra
    }

    options <-
      c(
        options,
        list(
          rcdt_min_angle = rcdt_min_angle,
          rcdt_max_edge = rcdt_max_edge,
          rcdt_max_n0 = rcdt_max_n0,
          rcdt_max_n1 = rcdt_max_n1
        )
      )
  }

  options
}



#' @title Refined Constrained Delaunay Triangulation
#'
#' @description
#' Computes a refined constrained Delaunay triangulation on R2 or S2.
#'
#' @param loc Input coordinates that should be part of the mesh
#' @param tv Initial triangulation, as a N-by-3 indec vector into `loc`
#' @param boundary,interior Objects supported by [fm_as_segm()].
#' If `boundary` is `numeric`, `fm_nonconvex_hull(loc, convex = boundary)` is
#' used.
#' @param extend `logical` or `list` specifying whether to extend the
#' data region, with parameters \describe{ \item{list("n")}{the number of edges
#' in the extended boundary (default=8)} \item{list("offset")}{the extension
#' distance.  If negative, interpreted as a factor relative to the approximate
#' data diameter (default=-0.10)} } Setting to `FALSE` is only useful in
#' combination `lattice` or `boundary`.
#' @param refine `logical` or `list` specifying whether to refine the
#' triangulation, with parameters \describe{ \item{list("min.angle")}{the
#' minimum allowed interior angle in any triangle.  The algorithm is guaranteed
#' to converge for `min.angle` at most 21 (default=`21`)}
#' \item{list("max.edge")}{the maximum allowed edge length in any triangle.  If
#' negative, interpreted as a relative factor in an ad hoc formula depending on
#' the data density (default=`Inf`)} \item{list("max.n.strict")}{the
#' maximum number of vertices allowed, overriding `min.angle` and
#' `max.edge` (default=-1, meaning no limit)} \item{list("max.n")}{the
#' maximum number of vertices allowed, overriding `max.edge` only
#' (default=-1, meaning no limit)} }
#' @param lattice An `fm_lattice_2d` object, generated by
#' [fm_lattice_2d()], specifying points on a regular lattice.
#' @param cutoff The minimum allowed distance between points.  Point at most as
#' far apart as this are replaced by a single vertex prior to the mesh
#' refinement step.
#' @param globe If non-NULL, an integer specifying the level of subdivision
#' for global mesh points, used with [fmesher_globe_points()]
#' @param quality.spec List of vectors of per vertex `max.edge` target
#' specification for each location in `loc`, `boundary/interior`
#' (`segm`), and `lattice`.  Only used if refining the mesh.
#' @param crs Optional crs object
#' @param ... Currently passed on to `fm_mesh_2d_inla` or converted to
#' [fmesher_rcdt()] options.
#' @examples
#' m <- fm_rcdt_2d(
#'   boundary = fm_nonconvex_hull(cbind(0, 0), convex = 5),
#'   rcdt_max_edge = 1
#' )
#'
#' @export
fm_rcdt_2d <-
  function(...) {
    fm_rcdt_2d_inla(...)
  }

#' @describeIn fm_rcdt_2d Legacy method for the `INLA::inla.mesh.create()`
#' interface
#' @export
fm_rcdt_2d_inla <-
  function(loc = NULL,
           tv = NULL,
           boundary = NULL,
           interior = NULL,
           extend = (missing(tv) || is.null(tv)),
           refine = FALSE,
           lattice = NULL,
           globe = NULL,
           cutoff = 1e-12,
           quality.spec = NULL,
           crs = NULL,
           ...) {
    crs.target <- crs
    if (!fm_crs_is_null(crs) &&
      fm_crs_is_geocent(crs)) {
      ## Build all geocentric meshes on a sphere, and transform afterwards,
      ## to allow general geoids.
      crs <- fm_crs("sphere")
    }

    if (!is.null(loc) && !is.matrix(loc)) {
      crs.loc <- fm_crs(loc)
    } else {
      crs.loc <- NULL
    }
    loc <- unify_loc_coords(loc)
    if (!fm_crs_is_null(crs.loc) && !fm_crs_is_null(crs)) {
      loc <- fm_transform(loc, crs = crs, passthrough = TRUE, crs0 = crs.loc)
      loc <- unify_loc_coords(loc)
    }

    if (!is.null(globe)) {
      loc.globe <- fmesher_globe_points(globe = globe)
      crs.globe <- fm_crs("sphere")
      if (!fm_crs_is_null(crs.globe) && !fm_crs_is_null(crs)) {
        loc.globe <- fm_transform(loc.globe, crs = crs, passthrough = TRUE, crs0 = crs.globe)
        loc.globe <- unify_loc_coords(loc.globe)
      }
      loc <- rbind(loc, loc.globe)
    }
    loc.n <- max(0L, nrow(loc))

    lattice.boundary <- NULL
    if (is.null(lattice) || !is.null(tv)) {
      if (!is.null(lattice)) {
        warning("Both 'lattice' and 'tv' specified.  Ignoring 'lattice'.")
      }
      lattice <- list(loc = NULL, segm = NULL)
      lattice.n <- 0L
    } else {
      lattice <- fm_as_lattice_2d(lattice)

      if (!fm_crs_is_null(fm_crs(lattice))) {
        lattice <- fm_transform(
          lattice,
          crs = crs,
          passthrough = TRUE
        )
      }
      if (NCOL(lattice$loc) == 2) {
        lattice$loc <- cbind(lattice$loc, 0.0)
      }
      if (NCOL(lattice$segm$loc) == 2) {
        lattice$segm$loc <- cbind(lattice$segm$loc, 0.0)
      }

      if (is.logical(extend) && !extend) {
        lattice.boundary <- lattice$segm
      }
    }
    lattice.n <- max(0L, nrow(lattice$loc))

    segm.n <- 0L
    if (!is.null(lattice.boundary) && is.null(boundary)) {
      boundary <- lattice.boundary
      lattice.boundary <- NULL
    }
    if (is.null(boundary)) {
      bnd <- NULL
      bnd_grp <- NULL
      loc.bnd <- matrix(0.0, 0, 3)
    } else {
      if (is.numeric(boundary)) {
        boundary <- fm_nonconvex_hull(loc, convex = boundary)
      } else {
        boundary <- fm_as_segm(boundary)
      }
      if (is.null(boundary$loc)) {
        boundary$loc <- loc
        boundary$crs <- fm_crs(crs)
      } else if (!fm_crs_is_null(crs)) {
        boundary <- fm_transform(boundary, crs = crs, passthrough = TRUE)
      }
      if (!is.null(lattice.boundary)) {
        boundary <-
          fm_segm_join(fm_as_segm_list(list(boundary, lattice.boundary)))
      }

      bnd <- segm.n + boundary$idx
      bnd_grp <- boundary$grp
      if (ncol(boundary$loc) == 2) {
        boundary$loc <- cbind(boundary$loc, 0.0)
      }
      segm.n <- segm.n + max(0L, nrow(boundary$loc))
      loc.bnd <- boundary$loc
    }

    if (is.null(interior)) {
      int <- NULL
      int_grp <- NULL
      loc.int <- matrix(0.0, 0, 3)
    } else {
      interior <- fm_as_segm(interior)
      if (is.null(interior$loc)) {
        interior$loc <- loc
        interior$crs <- fm_crs(crs)
      } else if (!fm_crs_is_null(crs)) {
        interior <- fm_transform(interior, crs = crs, passthrough = TRUE)
      }
      int <- segm.n + interior$idx
      int_grp <- interior$grp
      if (ncol(interior$loc) == 2) {
        interior$loc <- cbind(interior$loc, 0.0)
      }
      segm.n <- segm.n + max(0L, nrow(interior$loc))
      loc.int <- interior$loc
    }

    loc <- rbind(loc.bnd, loc.int, lattice$loc, loc)

    options <- handle_rcdt_options_inla(
      extend = extend,
      refine = refine,
      cutoff = cutoff,
      qulity.spec = quality.spec,
      ...,
      .n = list(
        segm = segm.n,
        lattice = lattice.n,
        loc = loc.n
      ),
      .loc = loc
    )

    if (!is.null(tv)) {
      tv <- tv + segm.n + lattice.n - 1L
    }
    if (!is.null(bnd)) {
      bnd <- bnd - 1L
    }
    if (!is.null(int)) {
      int <- int - 1L
    }
    result <- fmesher_rcdt(
      options = options,
      loc = loc, tv = tv,
      boundary = bnd, interior = int,
      boundary_grp = bnd_grp, interior_grp = int_grp
    )

    idx_C2R <- function(x) {
      x <- x + 1L
      x[x == 0] <- NA
      x
    }

    if (!fm_crs_is_null(crs) &&
      !fm_crs_is_identical(crs, crs.target)) {
      ## Target is a non-spherical geoid
      result[["s"]] <- fm_transform(result[["s"]], crs0 = crs, crs = crs.target)
      crs <- crs.target
    }

    split_idx <- function(idx, splits) {
      cumulative_splits <- c(0L, cumsum(splits))
      idx <- lapply(
        seq_along(splits),
        function(k) {
          if (splits[k] > 0) {
            idx[cumulative_splits[k] + seq_len(splits[k])]
          } else {
            NULL
          }
        }
      )
      names(idx) <- names(splits)
      idx
    }
    idx.all <- idx_C2R(result[["idx"]])
    idx <- split_idx(idx.all, c(segm = segm.n, lattice = lattice.n, loc = loc.n))

    m <- structure(
      list(
        meta = list(
          is.refined = !is.null(options[["rcdt_max_edge"]])
        ),
        manifold = result[["manifold"]],
        n = nrow(result[["s"]]),
        loc = result[["s"]],
        graph = list(
          tv = idx_C2R(result[["tv"]]),
          vt = idx_C2R(result[["vt"]]),
          tt = idx_C2R(result[["tt"]]),
          tti = idx_C2R(result[["tti"]]),
          vv = fm_as_dgCMatrix(result[["vv"]])
        ),
        segm = list(
          int = fm_segm(
            idx = idx_C2R(result[["segm.int.idx"]]),
            grp = result[["segm.int.grp"]],
            is.bnd = FALSE
          ),
          bnd = fm_segm(
            idx = idx_C2R(result[["segm.bnd.idx"]]),
            grp = result[["segm.bnd.grp"]],
            is.bnd = TRUE
          )
        ),
        idx = idx,
        crs = fm_crs(crs),
        n = nrow(result[["s"]])
      ),
      class = c("fm_mesh_2d", "inla.mesh")
    )

    remap_unused <- function(mesh) {
      ## Remap indices to remove unused vertices
      used <- !is.na(mesh$graph$vt)
      if (!all(used)) {
        used <- which(used)
        idx.map <- rep(NA, nrow(mesh$loc))
        idx.map[used] <- seq_len(length(used))
        mesh$loc <- mesh$loc[used, , drop = FALSE]
        mesh$graph$tv <-
          matrix(idx.map[as.vector(mesh$graph$tv)], nrow(mesh$graph$tv), 3)
        mesh$graph$vt <- mesh$graph$vt[used, , drop = FALSE]
        ## graph$tt  ## No change needed
        ## graph$tti ## No change needed
        mesh$graph$vv <- mesh$graph$vv[used, used, drop = FALSE]
        if (!is.null(mesh$idx$loc)) {
          mesh$idx$loc <- idx.map[mesh$idx$loc]
        }
        if (!is.null(mesh$idx$lattice)) {
          mesh$idx$lattice <- idx.map[mesh$idx$lattice]
        }
        if (!is.null(mesh$idx$segm)) {
          mesh$idx$segm <- idx.map[mesh$idx$segm]
        }
        mesh$segm$bnd$idx <-
          matrix(idx.map[mesh$segm$bnd$idx], nrow(mesh$segm$bnd$idx), 2)
        mesh$segm$int$idx <-
          matrix(idx.map[mesh$segm$int$idx], nrow(mesh$segm$int$idx), 2)
        mesh$segm$bnd$idx[mesh$segm$bnd$idx == 0L] <- NA
        mesh$segm$int$idx[mesh$segm$int$idx == 0L] <- NA
      }
      mesh
    }

    m <- remap_unused(m)

    m
  }

#' @describeIn fm_rcdt_2d Construct a plain Delaunay triangulation.
#' @export
fm_delaunay_2d <- function(loc, crs = NULL, ...) {
  if (is.null(crs) && !is.matrix(loc)) {
    crs <- fm_crs(loc)
  }
  loc <- unify_loc_coords(loc, crs = crs)

  hull <- grDevices::chull(loc[, 1], loc[, 2])
  bnd <- fm_segm(
    loc = loc[hull[rev(seq_along(hull))], , drop = FALSE],
    is.bnd = TRUE
  )
  mesh <- fm_rcdt_2d_inla(
    loc = loc,
    boundary = bnd,
    extend = list(n = 3),
    refine = FALSE,
    crs = crs,
    ...
  )
  return(mesh)
}


#' @title Make a 2D mesh object
#' @export
#' @param ... Currently passed on to `fm_mesh_2d_inla`
#' @family object creation and conversion
fm_mesh_2d <- function(...) {
  fm_mesh_2d_inla(...)
}

#' @describeIn fm_mesh_2d Legacy method for `INLA::inla.mesh.2d()`
#' Create a triangle mesh based on initial point locations, specified or
#' automatic boundaries, and mesh quality parameters.
#' @export
#'
#' @param loc Matrix of point locations to be used as initial triangulation
#' nodes.  Can alternatively be a `SpatialPoints` or
#' `SpatialPointsDataFrame` object.
#' @param loc.domain Matrix of point locations used to determine the domain
#' extent.  Can alternatively be a `SpatialPoints` or
#' `SpatialPointsDataFrame` object.
#' @param offset The automatic extension distance.  One or two values, for an
#' inner and an optional outer extension.  If negative, interpreted as a factor
#' relative to the approximate data diameter (default=-0.10???)
#' @param n The number of initial nodes in the automatic extensions
#' (default=16)
#' @param boundary one or more (as list) of [fm_segm()] objects, or objects
#' supported by [fm_as_segm()]
#' @param interior one object supported by [fm_as_segm()]
#' @param max.edge The largest allowed triangle edge length.  One or two
#' values.
#' @param min.angle The smallest allowed triangle angle.  One or two values.
#' (Default=21)
#' @param cutoff The minimum allowed distance between points.  Point at most as
#' far apart as this are replaced by a single vertex prior to the mesh
#' refinement step.
#' @param max.n.strict The maximum number of vertices allowed, overriding
#' `min.angle` and `max.edge` (default=-1, meaning no limit).  One or
#' two values, where the second value gives the number of additional vertices
#' allowed for the extension.
#' @param max.n The maximum number of vertices allowed, overriding
#' `max.edge` only (default=-1, meaning no limit).  One or two values,
#' where the second value gives the number of additional vertices allowed for
#' the extension.
#' @param plot.delay On Linux (and Mac if appropriate X11 libraries are
#' installed), specifying a nonnegative numeric value activates a rudimentary
#' plotting system in the underlying `fmesher` program, showing the
#' triangulation algorithm at work, with waiting time factor `plot.delay`
#' between each step.
#'
#' On all systems, specifying any negative value activates displaying the
#' result after each step of the multi-step domain extension algorithm.
#' @param crs An optional [fm_crs()], `sf::crs` or `sp::CRS` object
#' @return An `inla.mesh` object.
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @seealso [fm_rcdt_2d()], [fm_mesh_2d()], [fm_delaunay_2d()],
#' [fm_nonconvex_hull()], [fm_extensions()], [fm_refine()]
fm_mesh_2d_inla <- function(loc = NULL, ## Points to include in final triangulation
                            loc.domain = NULL, ## Points that determine the automatic domain
                            offset = NULL, ## Size of automatic extensions
                            n = NULL, ## Sides of automatic extension polygons
                            boundary = NULL, ## User-specified domains (list of length 2)
                            interior = NULL, ## User-specified constraints for the inner domain
                            max.edge = NULL,
                            min.angle = NULL, ## Angle constraint for the entire domain
                            cutoff = 1e-12, ## Only add input points further apart than this
                            max.n.strict = NULL,
                            max.n = NULL,
                            plot.delay = NULL,
                            crs = NULL,
                            ...) {
  ## plot.delay: Do plotting.
  ## NULL --> No plotting
  ## <0  --> Intermediate meshes displayed at the end
  ## >0   --> Dynamical fmesher plotting


  if ((missing(max.edge) || is.null(max.edge)) &&
    (missing(max.n.strict) || is.null(max.n.strict)) &&
    (missing(max.n) || is.null(max.n))) {
    max.edge <- NA
    #    stop("At least one of max.edge, max.n.strict, and max.n must be specified")
  }

  if (!is.null(crs)) {
    issphere <- fm_crs_is_identical(crs, fm_crs("sphere"))
    isgeocentric <- fm_crs_is_geocent(crs)
    if (isgeocentric) {
      crs.target <- crs
      crs <- fm_CRS("sphere")
    }
  }

  loc <- unify_loc_coords(loc, crs = crs)
  loc.domain <- unify_loc_coords(loc.domain, crs = crs)

  boundary <- fm_as_segm_list(boundary)
  interior <- fm_as_segm(interior)

  if (length(boundary) == 0) {
    list(NULL)
  }

  if (missing(offset) || is.null(offset)) {
    if (length(boundary) < 2) {
      offset <- -0.05
    } else {
      offset <- c(-0.05, -0.15)
    }
  }
  if (any(offset < 0) &&
    (fm_diameter(loc) +
      fm_diameter(loc.domain) +
      fm_diameter(interior) == 0.0)) {
    offset[offset < 0] <- 1
  }
  if (missing(n) || is.null(n)) {
    n <- c(8)
  }
  if (missing(max.edge) || is.null(max.edge)) {
    max.edge <- c(NA)
  }
  if (missing(min.angle) || is.null(min.angle)) {
    min.angle <- c(21)
  }
  if (missing(max.n.strict) || is.null(max.n.strict)) {
    max.n.strict <- c(NA)
  }
  if (missing(max.n) || is.null(max.n)) {
    max.n <- c(NA)
  }
  if (missing(cutoff) || is.null(cutoff)) {
    cutoff <- 1e-12
  }
  if (missing(plot.delay) || is.null(plot.delay)) {
    plot.delay <- NULL
  }

  num.layers <-
    max(c(
      length(boundary), length(offset), length(n),
      length(min.angle), length(max.edge),
      length(max.n.strict), length(max.n)
    ))
  if (num.layers > 2) {
    warning(paste("num.layers=", num.layers, " > 2 detected.  ",
      "Excess information ignored.",
      sep = ""
    ))
    num.layers <- 2
  }

  if (length(boundary) < num.layers) {
    boundary <- c(boundary, rep(list(NULL), num.layers - length(boundary)))
  }
  if (length(min.angle) < num.layers) {
    min.angle <- c(min.angle, min.angle)
  }
  if (length(max.n.strict) < num.layers) {
    max.n0 <- c(max.n.strict, max.n.strict)
  }
  if (length(max.n) < num.layers) {
    max.n <- c(max.n, max.n)
  }
  if (length(max.edge) < num.layers) {
    max.edge <- c(max.edge, max.edge)
  }
  if (length(offset) < num.layers) {
    offset <- c(offset, -0.15)
  }
  if (length(n) < num.layers) {
    n <- c(n, 16)
  }
  if (length(n) < num.layers) {
    n <- c(n, 16)
  }

  ## Unify the dimensionality of the boundary&interior segments input
  ## and optionally transform coordinates.
  boundary <- unify_segm_coords(boundary, crs = crs)
  interior <- unify_segm_coords(interior, crs = crs)

  ## Triangulate to get inner domain boundary
  ## Constraints included only to get proper domain extent
  ## First, attach the loc points to the domain definition set
  if (!is.null(loc) && !is.null(loc.domain)) {
    loc.domain <- rbind(loc.domain, loc)
  }
  mesh1 <-
    fm_rcdt_2d(
      loc = loc.domain,
      boundary = boundary[[1]],
      interior = interior,
      cutoff = cutoff,
      extend = list(n = n[1], offset = offset[1]),
      refine = FALSE,
      plot.delay = plot.delay,
      crs = crs
    )

  ## Save the resulting boundary
  boundary1 <- fm_segm(mesh1, boundary = TRUE)
  interior1 <- fm_segm(mesh1, boundary = FALSE)

  if (!is.null(plot.delay) && (plot.delay < 0)) {
    plot(mesh1)
  }

  ## Triangulate inner domain
  mesh2 <-
    fm_rcdt_2d(
      loc = loc,
      boundary = boundary1,
      interior = interior1,
      cutoff = cutoff,
      extend = FALSE, ## Should have no effect
      refine =
        list(
          min.angle = min.angle[1],
          max.edge = max.edge[1],
          max.edge.extra = max.edge[1],
          max.n.strict = max.n.strict[1],
          max.n = max.n[1]
        ),
      plot.delay = plot.delay,
      crs = crs
    )

  boundary2 <- fm_segm(mesh2, boundary = TRUE)
  interior2 <- fm_segm(mesh2, boundary = FALSE)

  if (!is.null(plot.delay) && (plot.delay < 0)) {
    plot(mesh2)
  }

  if (num.layers == 1) {
    if (!is.null(crs) && isgeocentric && !issphere) {
      mesh2$loc <- fm_transform(mesh2$loc, crs0 = mesh2$crs, crs = crs.target)
      mesh2$crs <- crs.target
    }
    return(mesh2)
  }

  ## Triangulate inner+outer domain
  mesh3 <-
    fm_rcdt_2d(
      loc = rbind(loc, mesh2$loc),
      boundary = boundary[[2]],
      interior = fm_segm(boundary2, interior2),
      cutoff = cutoff,
      extend = list(n = n[2], offset = offset[2]),
      refine =
        list(
          min.angle = min.angle[2],
          max.edge = max.edge[2],
          max.edge.extra = max.edge[2],
          max.n.strict = mesh2$n + max.n.strict[2],
          max.n = mesh2$n + max.n[2]
        ),
      plot.delay = plot.delay,
      crs = crs
    )

  ## Hide generated points, to match regular fm_rcdt_2d_inla output
  mesh3$idx$loc <- mesh3$idx$loc[seq_len(nrow(loc))]

  ## Obtain the corresponding segm indices.
  segm.loc <- matrix(0.0, 0, 3)
  for (k in seq_along(boundary)) {
    if (!is.null(boundary[[k]])) {
      segm.loc <- rbind(segm.loc, boundary[[k]]$loc)
    }
  }
  for (k in seq_along(interior)) {
    if (!is.null(interior[[k]])) {
      segm.loc <- rbind(segm.loc, interior[[k]]$loc)
    }
  }
  if (nrow(segm.loc) > 0) {
    proj <- fm_evaluator(mesh3, loc = segm.loc)$proj
    mesh3$idx$segm <- rep(NA, nrow(segm.loc))
    if (any(proj$ok)) {
      t.idx <- proj$t[proj$ok]
      tv.idx <- max.col(proj$bary[proj$ok, , drop = FALSE],
        ties.method = "first"
      )
      mesh3$idx$segm[proj$ok] <-
        mesh3$graph$tv[t.idx + nrow(mesh3$graph$tv) * (tv.idx - 1)]
    }
  } else {
    mesh3$idx$segm <- NULL
  }

  if (!is.null(crs) && isgeocentric && !issphere) {
    mesh3$loc <-
      fm_transform(mesh3$loc, crs0 = mesh3$crs, crs = crs.target)
    mesh3$crs <- crs.target
  }

  if (!is.null(plot.delay) && (plot.delay < 0)) {
    plot(mesh3)
  }

  return(mesh3)
}

#' @title Convert objects to `fm_mesh_2d`
#' @describeIn fm_as_mesh_2d Convert an object to `fm_mesh_2d`.
#' @param x Object to be converted.
#' @param ... Arguments passed on to submethods
#' @export
#' @family object creation and conversion
#' @export
fm_as_mesh_2d <- function(x, ...) {
  if (is.null(x)) {
    return(NULL)
  }
  UseMethod("fm_as_mesh_2d")
}
#' @describeIn fm_as_mesh_2d Convert each element of a list
#' @export
fm_as_mesh_2d_list <- function(x, ...) {
  fm_as_list(x, ..., .method = "fm_as_mesh_2d")
}
#' @rdname fm_as_mesh_2d
#' @param x Object to be converted
#' @export
fm_as_mesh_2d.fm_mesh_2d <- function(x, ...) {
  #  class(x) <- c("fm_mesh_2d", setdiff(class(x), "fm_mesh_2d"))
  x
}
#' @rdname fm_as_mesh_2d
#' @export
#' @method fm_as_mesh_2d inla.mesh
fm_as_mesh_2d.inla.mesh <- function(x, ...) {
  x[["crs"]] <- fm_crs(x[["crs"]])
  if (!is.null(x$segm$bnd)) {
    x$segm$bnd <- fm_as_fm(x$segm$bnd)
  }
  if (!is.null(x$segm$int)) {
    x$segm$int <- fm_as_fm(x$segm$int)
  }
  class(x) <- c("fm_mesh_2d", class(x))
  x
}


# fm_tensor ####

#' @title Make a tensor product function space
#' @export
#' @param x list of function space objects, such as [fm_mesh_2d()].
#' @param ... Currently unused
#' @family object creation and conversion
fm_tensor <- function(x, ...) {
  nn <- names(x)
  if (is.null(nn)) {
    nn <- as.character(seq_along(x))
  } else if (any(nn == "")) {
    stop("all or no elements of the list of function space objects need to be named.")
  }
  structure(
    list(fun_spaces = lapply(x, fm_as_fm)),
    class = "fm_tensor"
  )
}

#' @title Convert objects to `fm_tensor`
#' @describeIn fm_as_tensor Convert an object to `fm_tensor`.
#' @param x Object to be converted.
#' @param ... Arguments passed on to submethods
#' @export
#' @family object creation and conversion
#' @export
fm_as_tensor <- function(x, ...) {
  if (is.null(x)) {
    return(NULL)
  }
  UseMethod("fm_as_tensor")
}
#' @describeIn fm_as_tensor Convert each element of a list
#' @export
fm_as_tensor_list <- function(x, ...) {
  fm_as_list(x, ..., .method = "fm_as_tensor")
}
#' @rdname fm_as_tensor
#' @param x Object to be converted
#' @export
fm_as_tensor.fm_tensor <- function(x, ...) {
  #  class(x) <- c("fm_tensor", setdiff(class(x), "fm_tensor"))
  x
}





# fm_lattice_2d ####

#' Special coordinate mappings for `fm_mesh_2d` projections.
#'
#' Calculates coordinate mappings for `fm_mesh_2d` projections.
#' This is an internal function not intended for general use.
#'
#' @keywords internal
#' @param loc Coordinates to be mapped.
#' @param projection The projection type.
#' @param inverse If `TRUE`, `loc` are map coordinates and
#' coordinates in the mesh domain are calculated.  If `FALSE`, `loc`
#' are coordinates in the mesh domain and the forward map projection is
#' calculated.
#' @return For `fm_mesh_2d_map_lim`, a list:
#' \item{xlim }{X axis limits in the map domain}
#' \item{ylim }{Y axis limits in the map domain}
#' No attempt is
#' made to find minimal limits for partial spherical domains.
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @seealso [fm_evaluator()]
#' @export
fm_mesh_2d_map <- function(loc,
                           projection =
                             c("default", "longlat", "longsinlat", "mollweide"),
                           inverse = TRUE) {
  projection <- match.arg(projection)
  if (identical(projection, "default")) {
    return(loc)
  } else if (identical(projection, "longlat")) {
    if (inverse) {
      proj <-
        cbind(
          cos(loc[, 1] * pi / 180) * cos(loc[, 2] * pi / 180),
          sin(loc[, 1] * pi / 180) * cos(loc[, 2] * pi / 180),
          sin(loc[, 2] * pi / 180)
        )
    } else {
      proj <-
        cbind(
          atan2(loc[, 2], loc[, 1]) * 180 / pi,
          asin(pmax(-1, pmin(+1, loc[, 3]))) * 180 / pi
        )
    }
  } else if (identical(projection, "longsinlat")) {
    if (inverse) {
      coslat <- sqrt(pmax(0, 1 - loc[, 2]^2))
      proj <-
        cbind(
          cos(loc[, 1] * pi / 180) * coslat,
          sin(loc[, 1] * pi / 180) * coslat,
          loc[, 2]
        )
    } else {
      proj <-
        cbind(
          atan2(loc[, 2], loc[, 1]) * 180 / pi,
          loc[, 3]
        )
    }
  } else if (identical(projection, "mollweide")) {
    if (inverse) {
      ok <- ((loc[, 1]^2 + 4 * loc[, 2]^2) <= 4)
      cos.theta <- sqrt(pmax(0, 1 - loc[ok, 2]^2))
      theta <- atan2(loc[ok, 2], cos.theta)
      sin.lat <- (2 * theta + sin(2 * theta)) / pi
      cos.lat <- sqrt(pmax(0, 1 - sin.lat^2))
      lon <- loc[ok, 1] * pi / 2 / (cos.theta + (cos.theta == 0))
      lon[cos.theta == 0] <- pi / 2 * sign(theta[cos.theta == 0])
      proj <- matrix(NA, nrow(loc), 3)
      proj[ok, ] <- cbind(cos(lon) * cos.lat, sin(lon) * cos.lat, sin.lat)
    } else {
      lon <- atan2(loc[, 2], loc[, 1])
      z <- pmin(1, pmax(-1, loc[, 3]))
      sin.theta <- z
      cos.theta <- sqrt(pmax(0, 1 - sin.theta^2))
      ## NR-solver for sin.theta.
      ## Typically finishes after at most 7 iterations.
      ## When cos.theta=0, sin.theta is already correct, +/- 1.
      nook <- (cos.theta > 0)
      for (k in 1:20) {
        if (any(nook)) {
          delta <-
            (atan2(sin.theta[nook], cos.theta[nook]) +
              sin.theta[nook] * cos.theta[nook] - pi / 2 * z[nook]) /
              (2 * cos.theta[nook])
          sin.theta[nook] <- sin.theta[nook] - delta
          cos.theta[nook] <- sqrt(1 - sin.theta[nook]^2)
          nook[nook] <- (abs(delta) > 1e-14)
        }
      }
      proj <- cbind(2 * lon / pi * cos.theta, sin.theta)
    }
  } else {
    stop(paste("Unknown projection '", projection, "'.", sep = ""))
  }
  return(proj)
}


#' @export
#' @describeIn fm_mesh_2d_map Projection extent limit calculations
fm_mesh_2d_map_lim <- function(loc = NULL,
                               projection =
                                 c("default", "longlat", "longsinlat", "mollweide")) {
  projection <- match.arg(projection)
  if (identical(projection, "default")) {
    if (is.null(loc)) {
      lim <- list(xlim = c(0, 1), ylim = c(0, 1))
    } else {
      lim <-
        list(
          xlim = range(loc[, 1], na.rm = TRUE),
          ylim = range(loc[, 2], na.rm = TRUE)
        )
    }
  } else if (identical(projection, "longlat")) {
    lim <- list(xlim = c(-180, 180), ylim = c(-90, 90))
  } else if (identical(projection, "longsinlat")) {
    lim <- list(xlim = c(-180, 180), ylim = c(-1, 1))
  } else if (identical(projection, "mollweide")) {
    lim <- list(xlim = c(-2, 2), ylim = c(-1, 1))
  } else {
    stop(paste("Unknown projection '", projection, "'.", sep = ""))
  }
  return(lim)
}




#' @title Make a lattice object
#' @export
#' @param ... Currently passed on to `inla.mesh.lattice`
#' @family object creation and conversion
fm_lattice_2d <- function(...) {
  UseMethod("fm_lattice_2d")
}

#' Lattice grids for inla.mesh
#'
#' Construct a lattice grid for [fm_mesh_2d()]
#'
#' @param x vector or grid matrix of x-values
#' @param y vector of grid matrix of y-values
#' @param z if x is a matrix, a grid matrix of z-values
#' @param dims the size of the grid, length 2 vector
#' @param units One of `c("default", "longlat", "longsinlat", "mollweide")`
#' or NULL (equivalent to `"default"`).
#' @param crs An optional `fm_crs`, `sf::st_crs`, or `sp::CRS` object
#' @return An `fm_lattice_2d` object.
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @seealso [fm_mesh_2d()]
#' @examples
#'
#' lattice <- fm_lattice_2d(seq(0, 1, length.out = 17), seq(0, 1, length.out = 10))
#'
#' ## Use the lattice "as-is", without refinement:
#' mesh <- fm_rcdt_2d_inla(lattice = lattice, boundary = lattice$segm)
#' mesh <- fm_rcdt_2d_inla(lattice = lattice, extend = FALSE)
# plot(mesh)
#'
#' ## Refine the triangulation, with limits on triangle angles and edges:
#' mesh <- fm_rcdt_2d(
#'   lattice = lattice,
#'   refine = list(max.edge = 0.08),
#'   extend = FALSE
#' )
# plot(mesh)
#'
#' ## Add an extension around the lattice, but maintain the lattice edges:
#' mesh <- fm_rcdt_2d(
#'   lattice = lattice,
#'   refine = list(max.edge = 0.08),
#'   interior = lattice$segm
#' )
# plot(mesh)
#'
#' ## Only add extension:
#' mesh <- fm_rcdt_2d(lattice = lattice, refine = list(max.edge = 0.08))
# plot(mesh)
#' @rdname fm_lattice_2d
#' @export
fm_lattice_2d.default <- function(
    x = seq(0, 1, length.out = 2),
    y = seq(0, 1, length.out = 2),
    z = NULL,
    dims =
      if (is.matrix(x)) {
        dim(x)
      } else {
        c(length(x), length(y))
      },
    units = NULL,
    crs = NULL,
    ...) {
  if (is.null(crs)) {
    units <- match.arg(units, c("default", "longlat", "longsinlat", "mollweide"))

    lim <- fm_mesh_2d_map_lim(projection = units)
    xlim <- lim$xlim
    ylim <- lim$ylim
  } else { ## !is.null(crs)
    if (!is.null(units)) {
      stop("Only one of 'units' and 'crs' can be non-null.")
    }
    bounds <- fm_crs_bounds(crs, warn.unknown = TRUE)
    xlim <- bounds$xlim
    ylim <- bounds$ylim
  }

  if (missing(x) && !missing(dims)) {
    x <- seq(xlim[1], xlim[2], length.out = dims[1])
  }
  if (missing(y) && !missing(dims)) {
    y <- seq(ylim[1], ylim[2], length.out = dims[2])
  }
  dims <- as.integer(dims)

  if (is.matrix(x)) {
    if (!identical(dims, dim(x)) ||
      !identical(dims, dim(y)) ||
      (is.matrix(z) && !identical(dims, dim(z)))) {
      stop("The size of matrices 'x', 'y', and 'z' must match 'dims'.")
    }
    loc <- cbind(as.vector(x), as.vector(y), as.vector(z))
    x <- NULL
    y <- NULL
  } else {
    if (!identical(dims[1], length(x)) ||
      !identical(dims[2], length(y))) {
      stop(paste("The lengths of vectors 'x' and 'y' (",
        length(x), ",", length(y),
        ") must match 'dims' (", dims[1], ",", dims[2], ").",
        sep = ""
      ))
    }
    loc <- (cbind(
      rep(x, times = dims[2]),
      rep(y, each = dims[1])
    ))
  }
  if (!is.double(loc)) {
    storage.mode(loc) <- "double"
  }

  if (is.null(crs)) {
    loc <- fm_mesh_2d_map(loc = loc, projection = units, inverse = TRUE)
  }

  ## Construct lattice boundary
  segm.idx <- (c(
    1:(dims[1] - 1),
    dims[1] * (1:(dims[2] - 1)),
    dims[1] * dims[2] - (0:(dims[1] - 2)),
    dims[1] * ((dims[2] - 1):1) + 1
  ))
  segm.grp <- (c(
    rep(1L, dims[1] - 1),
    rep(2L, dims[2] - 1),
    rep(3L, dims[1] - 1),
    rep(4L, dims[2] - 1)
  ))

  segm <- fm_segm(
    loc = loc[segm.idx, , drop = FALSE],
    grp = segm.grp,
    is.bnd = TRUE,
    crs = crs
  )

  lattice <- list(dims = dims, x = x, y = y, loc = loc, segm = segm, crs = crs)
  class(lattice) <- c("fm_lattice_2d", "inla.mesh.lattice")
  return(lattice)
}

#' @title Convert objects to `fm_lattice_2d`
#' @describeIn fm_as_lattice_2d Convert an object to `fm_lattice_2d`.
#' @param x Object to be converted.
#' @param ... Arguments passed on to submethods
#' @export
#' @family object creation and conversion
#' @export
fm_as_lattice_2d <- function(...) {
  UseMethod("fm_as_lattice_2d")
}
#' @describeIn fm_as_lattice_2d Convert each element of a list
#' @export
fm_as_lattice_2d_list <- function(x, ...) {
  fm_as_list(x, ..., .method = "fm_as_lattice_2d")
}
#' @rdname fm_as_lattice_2d
#' @param x Object to be converted
#' @export
fm_as_lattice_2d.fm_lattice_2d <- function(x, ...) {
  #  class(x) <- c("fm_lattice_2d", setdiff(class(x), "fm_lattice_2d"))
  x
}
#' @rdname fm_as_lattice_2d
#' @param x Object to be converted
#' @export
#' @method fm_as_lattice_2d inla.mesh.lattice
fm_as_lattice_2d.inla.mesh.lattice <- function(x, ...) {
  x[["crs"]] <- fm_crs(x[["crs"]])
  x[["segm"]] <- fm_as_fm(x[["segm"]])
  class(x) <- c("fm_lattice_2d", class(x))
  x
}



# fm_bbox


#' @title Bounding box class
#'
#' @description
#' Simple class for handling bounding box information
#' @param x Object from which to extract bounding box information
#' @param ... Passed on to sub-methods
#' @export
#' @examples
#' fm_bbox(matrix(1:6, 3, 2))
fm_bbox <- function(...) {
  UseMethod("fm_bbox")
}

#' @describeIn fm_bbox Construct a bounding box from
#' precomputed interval information, stored as a list,
#' `list(lim = list(xlim, ylim))`.
#' @export
fm_bbox.list <- function(x, ...) {
  structure(
    list(
      lim = x
    ),
    class = "fm_bbox"
  )
}

#' @rdname fm_bbox
#' @export
fm_bbox.matrix <- function(x, ...) {
  fm_bbox(lapply(
    seq_len(ncol(x)),
    function(k) {
      if (all(is.na(x[, k]))) {
        c(NA_real_, NA_real_)
      } else {
        range(x[, k], na.rm = TRUE)
      }
    }
  ))
}

#' @rdname fm_bbox
#' @export
fm_bbox.fm_bbox <- function(x, ...) {
  x
}

#' @rdname fm_bbox
#' @export
fm_bbox.fm_mesh_2d <- function(x, ...) {
  fm_bbox(x[["loc"]])
}

#' @rdname fm_bbox
#' @export
fm_bbox.fm_segm <- function(x, ...) {
  if (is.null(x[["loc"]])) {
    return(fm_bbox(list()))
  }
  fm_bbox(x[["loc"]])
}

#' @rdname fm_bbox
#' @export
fm_bbox.fm_lattice_2d <- function(x, ...) {
  fm_bbox(x[["loc"]])
}

#' @rdname fm_bbox
#' @export
fm_bbox.sf <- function(x, ...) {
  fm_bbox(sf::st_geometry(x))
}

#' @rdname fm_bbox
#' @export
fm_bbox.sfg <- function(x, ...) {
  fm_bbox(sf::st_sfc(x))
}

#' @rdname fm_bbox
#' @export
fm_bbox.sfc <- function(x, ...) {
  loc <- sf::st_coordinates(x)
  loc <- loc[, intersect(colnames(loc), c("X", "Y", "Z", "M")), drop = FALSE]
  fm_bbox(loc)
}

#' @rdname fm_bbox
#' @export
fm_bbox.bbox <- function(x, ...) {
  # sf bbox objects are length 4 vectors
  fm_bbox(matrix(x, 2, 2, byrow = TRUE))
}

#' @rdname fm_bbox
#' @export
fm_bbox.inla.mesh <- function(x, ...) {
  fm_bbox(fm_as_fm(x))
}

#' @rdname fm_bbox
#' @export
fm_bbox.inla.mesh.segment <- function(x, ...) {
  fm_bbox(fm_as_fm(x))
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
