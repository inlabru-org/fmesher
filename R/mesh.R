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
#' \donttest{
#' if (require("ggplot2", quietly = TRUE) &&
#'   require("inlabru", quietly = TRUE)) {
#'   data("mrsea", package = "inlabru")
#'   pxl <- fm_pixels(mrsea$mesh,
#'     nx = 50, ny = 50, mask = mrsea$boundary,
#'     format = "terra"
#'   )
#'   ggplot() +
#'     gg(pxl, fill = "grey", alpha = 0.5) +
#'     gg(mrsea$mesh)
#' }
#' }
#'
fm_pixels <- function(mesh, nx = 150, ny = 150, mask = TRUE,
                      format = "sf") {
  format <- match.arg(format, c("sf", "terra", "sp"))
  if (!identical(mesh$manifold, "R2")) {
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
    pixels <- pixels[lengths(pixels_within) > 0, , drop = FALSE]
  }

  if (identical(format, "sp")) {
    pixels <- as(pixels, "Spatial")
    pixels <- as(pixels, "SpatialPixelsDataFrame")
  } else if (identical(format, "terra")) {
    pixels <- as(pixels, "Spatial")
    pixels <- as(pixels, "SpatialPixels")
    pixels <- terra::rast(pixels)
    pixels$.mask <- as.vector(as.matrix(pixels_within))
  }

  pixels
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

  boundary2 <- split.edges(fm_segm(mesh), n = n)

  if (identical(mesh$manifold, "S2")) {
    radius <- mean(rowSums(mesh$loc^2)^0.5)
    renorm <- function(loc) {
      loc * (radius / rowSums(loc^2)^0.5)
    }
    tri.inner.loc <- renorm(tri.inner.loc)
    tri.edges2$loc <- renorm(tri.edges2$loc)
    boundary2$loc <- renorm(boundary2$loc)
  }

  mesh2 <- INLA::inla.mesh.create(
    loc = tri.inner.loc,
    interior = tri.edges,
    boundary = boundary2,
    refine = list(
      min.angle = 0,
      max.edge = Inf
    ),
    crs = fm_CRS(mesh)
  )

  mesh2
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


#' @title Extract vertex locations from an `inla.mesh`
#'
#' @description Extracts the vertices of an `inla.mesh` object.
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
#' \donttest{
#' if (require("ggplot2", quietly = TRUE) &&
#'   require("inlabru", quietly = TRUE)) {
#'   data("mrsea", package = "inlabru")
#'   vrt <- fm_vertices(mrsea$mesh, format = "sp")
#'   ggplot() +
#'     gg(mrsea$mesh) +
#'     gg(vrt, color = "red")
#' }
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
#' \donttest{
#' if (require("ggplot2", quietly = TRUE) &&
#'   require("inlabru", quietly = TRUE)) {
#'   data("mrsea", package = "inlabru")
#'   vrt <- fm_centroids(mrsea$mesh, format = "sp")
#'   ggplot() +
#'     gg(mrsea$mesh) +
#'     gg(vrt, color = "red")
#' }
#' }
#'
fm_centroids <- function(x, format = NULL) {
  ## Extract triangle centroids
  loc <- (x$loc[x$graph$tv[, 1], , drop = FALSE] +
    x$loc[x$graph$tv[, 2], , drop = FALSE] +
    x$loc[x$graph$tv[, 3], , drop = FALSE]) / 3

  if (identical(x$manifold, "S2")) {
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
  } else if (identical(mesh$manifold, "S2")) {
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





#' @title Compute barycentric coordinates
#'
#' @description Identify triangles and compute barycentric coordinates
#'
#' @param mesh `fm_mesh_2d` or `inla.mesh` object
#' @param loc Points for which to identify the containing triangle, and
#' corresponding barycentric coordinates. May be raw matrix coordinates, `sf`, or `sp`
#' point information.
#' @param crs Optional crs information for `loc`
#'
#' @return A list with elements `t` (vector of triangle indices) and `bary`
#' (3-column matrix of barycentric coordinates). Points that were not found
#' give `NA` entries in `t` and `bary`.
#'
#' @export
fm_bary <- function(mesh, loc, crs = NULL) {
  loc <- fm_onto_mesh(mesh, loc, crs = crs)

  # Avoid sphere accuracy issues by scaling to unit sphere
  scale <- 1
  if (identical(mesh$manifold, "S2")) {
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
#' @method fm_fem inla.mesh
fm_fem.inla.mesh <- function(mesh, order = 2, ...) {
  fm_fem(fm_as_fm(mesh), order = order, ...)
}

#' @rdname fm_fem
#' @export
fm_fem.fm_mesh_2d <- function(mesh, order = 2, ...) {
  result <- fmesher_fem(
    mesh_loc = mesh$loc,
    mesh_tv = mesh$graph$tv - 1L,
    fem_order_max = order,
    options = list()
  )
  result
}



#' @title Split lines at triangle edges
#'
#' @description Compute intersections between line segments and triangle edges,
#' and filter out segment of length zero.
#'
#' @param mesh An `fm_mesh_2d` or `inla.mesh` object
#' @param sp Start points of lines
#' @param ep End points of lines
#' @param ... Passed on to the method for `fm_mesh_2d`.
#' @return List of start and end points resulting from splitting the given lines:
#' `list(sp, ep, split.origin, idx, split.loc)`
#'
#' @author Fabian E. Bachl \email{f.e.bachl@@bath.ac.uk}
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @keywords internal
#'
#' @export
fm_split_lines <- function(mesh, ...) {
  UseMethod("fm_split_lines")
}


#' @rdname fm_split_lines
#' @export
fm_split_lines.inla.mesh <- function(mesh, ...) {
  fm_split_lines(fm_as_mesh_2d(mesh), ...)
}

#' @rdname fm_split_lines
#' @export
fm_split_lines.fm_mesh_2d <- function(mesh, sp, ep, ...) {
  idx <- seq_len(NROW(sp))
  if (NROW(sp) > 0) {
    # Filter out segments not on the mesh
    t1 <- fm_bary(mesh, loc = as.matrix(sp))$t
    t2 <- fm_bary(mesh, loc = as.matrix(ep))$t
    keep <- !(is.na(t1) | is.na(t2))
    # if (any(!keep)) { warning("points outside boundary! filtering...")}
    sp <- sp[keep, , drop = FALSE]
    ep <- ep[keep, , drop = FALSE]
    idx <- idx[keep]
  }

  if (NROW(sp) == 0) {
    return(list(
      sp = sp, ep = ep,
      split.origin = NULL,
      idx = idx,
      split.loc = NULL
    ))
  }

  loc <- as.matrix(rbind(sp, ep))

  # Split the segments into parts
  if (NCOL(loc) == 2) {
    loc <- cbind(loc, rep(0, NROW(loc)))
  }
  np <- dim(sp)[1]
  sp.idx <- cbind(seq_len(np), np + seq_len(np))
  splt <- fmesher_split_lines(
    mesh_loc = mesh$loc,
    mesh_tv = mesh$graph$tv - 1L,
    loc = loc,
    idx = sp.idx - 1L,
    options = list()
  )
  indexoutput <- list("split.idx", "split.t", "split.origin")
  for (name in intersect(names(splt), indexoutput)) {
    splt[[name]] <- splt[[name]] + 1L
  }
  # plot(data$mesh)
  # points(loc)
  # points(splt$split.loc,col="blue)

  # Start points of new segments
  sp <- splt$split.loc[splt$split.idx[, 1], seq_len(dim(sp)[2]), drop = FALSE]
  # End points of new segments
  ep <- splt$split.loc[splt$split.idx[, 2], seq_len(dim(ep)[2]), drop = FALSE]
  idx <- idx[splt$split.idx[, 1]]
  origin <- splt$split.origin

  # Filter out zero length segments
  sl <- apply((ep - sp)^2, MARGIN = 1, sum)
  sp <- sp[!(sl == 0), , drop = FALSE]
  ep <- ep[!(sl == 0), , drop = FALSE]
  origin <- origin[!(sl == 0)]
  idx <- idx[!(sl == 0)]

  return(list(
    sp = sp, ep = ep,
    split.origin = origin,
    idx = idx,
    split.loc = splt$split.loc
  ))
}







#' @title Compute an extension of a spatial object
#'
#' @description
#' Constructs a potentially nonconvex extension of a spatial object by
#' performing dilation by `convex + concave` followed by
#' erosion by `concave`. This is equivalent to dilation by `convex` followed
#' by closing (dilation + erosion) by `concave`.
#'
#' @param x A spatial object
#' @param ... Arguments passed on to the sub-methods
#' @param convex How much to extend
#' @param concave The minimum allowed reentrant curvature. Default equal to `convex`
#' @param preserveTopology logical; argument to `sf::st_simplify()`
#' @param dTolerance If not null, the `dTolerance` argument to `sf::st_simplify()`
#' @returns An extended object
#' @references Gonzalez and Woods (1992), Digital Image Processing
#' @export
#' @examples
#' inp <- sf::st_as_sf(as.data.frame(matrix(1:6, 3, 2)), coords = 1:2)
#' out <- fm_nonconvex_hull(inp, convex = 1)
#' plot(out)
fm_nonconvex_hull <- function(x, ...) {
  UseMethod("fm_nonconvex_hull")
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

#' @describeIn fm_nonconvex_hull Legacy method for `INLA::inla.nonconvex.hull()`
#' @export
fm_nonconvex_hull_inla <- function(...) {
  fm_as_segm(INLA::inla.nonconvex.hull(...))
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
                                  ...) {
  if ((convex < 0) || (concave < 0)) {
    approx_diameter <- fm_diameter(x)
    if (convex < 0) {
      convex <- approx_diameter * abs(convex)
    }
    if (concave < 0) {
      concave <- approx_diameter * abs(concave)
    }
  }

  nQuadSegs <- 64
  y <- sf::st_buffer(x, dist = convex + concave, nQuadSegs = nQuadSegs)
  y <- sf::st_union(y)
  y <- sf::st_buffer(y, dist = -concave, nQuadSegs = nQuadSegs)
  if (!is.null(dTolerance)) {
    y <- sf::st_simplify(y,
      preserveTopology = preserveTopology,
      dTolerance = dTolerance
    )
  }
  y <- sf::st_union(y)
  y
}

#' @title Compute approximate spatial object diameter
#' @param x A spatial object
#' @param ... Currently unused
#' @export
fm_diameter <- function(x, ...) {
  UseMethod("fm_diameter")
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
fm_diameter.matrix <- function(x, ...) {
  if (ncol(x) == 1) {
    lim <- range(x)
    approx_diameter <- diff(lim)
  } else if (ncol(x) == 2) {
    lim <- rbind(range(x[, 1]), range(x[, 2]))
    approx_diameter <- max(diff(lim[1, ]), diff(lim[2, ]))
  } else if (ncol(x) >= 3) {
    lim <- rbind(range(x[, 1]), range(x[, 2]), range(x[, 3]))
    approx_diameter <- max(diff(lim[1, ]), diff(lim[2, ]), diff(lim[3, ]))
  }
  approx_diameter
}


#' @title Compute extensions of a spatial object
#'
#' @description
#' Constructs a potentially nonconvex extension of a spatial object by
#' performing dilation by `convex + concave` followed by
#' erosion by `concave`. This is equivalent to dilation by `convex` followed
#' by closing (dilation + erosion) by `concave`.
#' @seealso [fm_nonconvex_hull()]
#'
#' @param x A spatial object
#' @param convex numeric vector; How much to extend
#' @param concave numeric vector; The minimum allowed reentrant curvature. Default equal to `convex`
#' @param dTolerance If not null, the `dTolerance` argument to `sf::st_simplify()`,
#' passed on to [fm_nonconvex_hull()].
#' The default is `pmin(convex, concave) / 40`, chosen to
#' give approximately 4 or more subsegments per circular quadrant.
#' @param ... Optional further arguments to pass on to [fm_nonconvex_hull()].
#' @returns A list of `sfg` objects.
#' @export
#' @examples
#' if (fm_safe_inla()) {
#'   inp <- sf::st_as_sf(as.data.frame(matrix(1:6, 3, 2)), coords = 1:2)
#'   bnd <- fm_extensions(inp, convex = c(0.75, 2))
#'   plot(fm_mesh_2d(boundary = bnd, max.edge = c(0.25, 1)), asp = 1)
#' }
#'
fm_extensions <- function(x,
                          convex = -0.15,
                          concave = convex,
                          dTolerance = NULL,
                          ...) {
  if (any(convex < 0) || any(concave < 0) || any(dTolerance < 0)) {
    approx_diameter <- fm_diameter(x)
  }
  len <- max(length(convex), length(concave), length(dTolerance))
  scale_fun <- function(val) {
    if (any(val < 0)) {
      val[val < 0] <- approx_diameter * abs(val[val < 0])
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
  # Match exact manifold?
  if (x[["manifold"]] %in% type) {
    return(TRUE)
  }
  # Match space name or dimension?
  m <- intersect(c("M", "R", "S"), type)
  d <- intersect(as.character(seq_len(3)), type)
  return(grepl(paste0(c(m, d), collapse = "|"), x[["manifold"]]))
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
#' [inla.mesh.create()].  If `is.bnd==TRUE`, defaults to linking
#' all the points in `loc`, as `c(1:nrow(loc),1L)`, otherwise
#' `1:nrow(loc)`.
#' @param grp Vector of group labels for each segment.  Set to `NULL` to
#' let the labels be chosen automatically in a call to
#' [inla.mesh.create()].
#' @param is.bnd `TRUE` if the segments are boundary segments, otherwise
#' `FALSE`.
#' @param crs An optional `fm_crs()`, `sf::st_crs()` or `sp::CRS()` object
#' @export
fm_segm.default <- function(loc = NULL, idx = NULL, grp = NULL, is.bnd = TRUE,
                            crs = NULL, ...) {
  if (is.null(loc) && is.null(idx)) {
    stop("At most one of 'loc' and 'idx' may be missing or null.")
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
    grp <- matrix(grp, min(length(grp), nrow(idx)), 1)
    if (nrow(grp) < nrow(idx)) {
      grp <- (matrix(
        c(
          as.vector(grp),
          rep(grp[nrow(grp)], nrow(idx) - length(grp))
        ),
        nrow(idx), 1
      ))
    }
    storage.mode(grp) <- "integer"
  } else {
    grp <- NULL
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
      grp <- grp[-i, , drop = FALSE]
    }
  }

  if (!is.null(loc)) {
    ## Identify unused locations and remap indices accordingly.
    idx.new <- rep(0L, nrow(loc))
    idx.new[as.vector(idx)] <- 1L
    loc <- loc[idx.new == 1L, , drop = FALSE]
    idx.new[idx.new == 1L] <- seq_len(sum(idx.new))
    idx <- (matrix(idx.new[as.vector(idx)],
      nrow = nrow(idx),
      ncol = ncol(idx)
    ))
  }

  ret <- list(loc = loc, idx = idx, grp = grp, is.bnd = is.bnd, crs = crs)
  class(ret) <- c("fm_segm", "inla.mesh.segment")
  return(ret)
}

#' @describeIn fm_segm Join multiple `fm_segm` objects into a single `fm_segm`
#' object.
#' @param grp.default When joining segments, use this group label for segments
#' that have `grp == NULL`.
#' @export
fm_segm.fm_segm <- function(..., grp.default = 0) {
  segm <- fm_as_segm_list(list(...))

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
    function(x) segm[[x]]$idx + cumNloc[x]
  ))
  grp <- unlist(lapply(
    seq_along(segm),
    function(x) {
      if (is.null(segm[[x]]$grp)) {
        rep(grp.default, Nidx[x])
      } else {
        segm[[x]]$grp
      }
    }
  ))
  is.bnd <- vapply(segm, function(x) x$is.bnd, TRUE)
  if (!all(is.bnd) || (any(!is.bnd) && !all(!is.bnd))) {
    warning("Inconsistent 'is.bnd' attributes.  Setting 'is.bnd=FALSE'.")
    is.bnd <- FALSE
  } else {
    is.bnd <- all(is.bnd)
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

#' @describeIn fm_as_segm Convert each element of a list
#' @export
fm_as_segm_list <- function(x, ...) {
  fm_as_list(x, ..., .method = "fm_as_segm")
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

#' @rdname fm_as_segm
#' @export
fm_as_segm.matrix <- function(x, crs = NULL, ...) {
  fm_segm(loc = x, crs = fm_CRS(crs))
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
  if (is.null(x)) {
    return(NULL)
  }
  UseMethod("fm_as_fm")
}
# @description fm_as_list Convert each element of a list, or convert a single
# object and return in a list
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
  fm_as_list(x, ..., .method = "fm_as_fm")
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

#' @title Make a 1D mesh object
#' @export
#' @param ... Currently passed on to `inla.mesh.1d`
#' @family object creation and conversion
fm_mesh_1d <- function(...) {
  fm_as_mesh_1d(INLA::inla.mesh.1d(...))
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
  if (is.null(loc) || (nrow(loc) == 0)) {
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
    loc <- as.matrix(loc)
  }
  if (ncol(loc) == 2) {
    loc <- cbind(loc, 0.0)
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
    cet_sides <- ifelse(is.null(extend$n), 16, extend$n)
    cet_margin <- ifelse(is.null(extend$offset), -0.1, extend$offset)
  }
  options <- c(options, list(cet_sides = cet_sides, cet_margin = cet_margin))

  if (inherits(refine, "list")) {
    rcdt_min_angle <- 0
    rcdt_max_edge <- 0
    ## Multiply by to ensure cover S2 domains
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
      rcdt_max_n0 <- refine$max.n.strict
    } else {
      rcdt_max_n0 <- -1
    }
    if (!is.null(refine[["max.n"]]) &&
      !is.na(refine$max.n)) {
      rcdt_max_n1 <- refine$max.n
    } else {
      rcdt_max_n1 <- -1
    }
    is.refined <- TRUE

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

#' @describeIn fm_rcdt_2d Legacy method for `INLA::inla.mesh.create()`
#' @export
fm_rcdt_2d_inla <-
  function(loc = NULL,
           tv = NULL,
           boundary = NULL,
           interior = NULL,
           crs = NULL,
           globe = NULL,
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

    segm.n <- 0L
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
      if (!fm_crs_is_null(crs)) {
        boundary <- fm_transform(boundary, crs = crs, passthrough = TRUE)
      }
      bnd <- segm.n + boundary$idx - 1L
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
      if (!fm_crs_is_null(crs)) {
        interior <- fm_transform(interior, crs = crs, passthrough = TRUE)
      }
      int <- segm.n + interior$idx - 1L
      int_grp <- interior$grp
      if (ncol(interior$loc) == 2) {
        interior$loc <- cbind(interior$loc, 0.0)
      }
      segm.n <- segm.n + max(0L, nrow(interior$loc))
      loc.int <- interior$loc
    }

    lattice.n <- 0L
    loc.lattice <- matrix(0, 0, 3)

    loc <- rbind(loc.bnd, loc.int, loc.lattice, loc)

    options <- handle_rcdt_options_inla(...,
                                        .n = list(
                                          segm = segm.n,
                                          lattice = lattice.n,
                                          loc = loc.n
                                        ),
                                        .loc = loc)

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
#' @seealso [fm_rcdt_2d()], [fm_mesh_2d()], [inla.delaunay()],
#' [fm_nonconvex_hull()], [fm_extensions()]
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
  ##  fm_as_mesh_2d(INLA::inla.mesh.2d(..., boundary = boundary, interior = interior))

  ## plot.delay: Do plotting.
  ## NULL --> No plotting
  ## <0  --> Intermediate meshes displayed at the end
  ## >0   --> Dynamical fmesher plotting


  if ((missing(max.edge) || is.null(max.edge)) &&
    (missing(max.n.strict) || is.null(max.n.strict)) &&
    (missing(max.n) || is.null(max.n))) {
    stop("At least one of max.edge, max.n.strict, and max.n must be specified")
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
    boundary <- c(boundary, list(NULL))
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
    return(invisible(mesh2))
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

  ## Hide generated points, to match regular inla.mesh.create output
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

  return(invisible(mesh3))
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

#' @title Make a lattice object
#' @export
#' @param ... Currently passed on to `inla.mesh.lattice`
#' @family object creation and conversion
fm_lattice_2d <- function(...) {
  UseMethod("fm_lattice_2d")
}

#' @rdname fm_lattice_2d
#' @export
fm_lattice_2d.default <- function(...) {
  fm_as_lattice_2d(INLA::inla.mesh.lattice(...))
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
  class(x) <- c("fm_lattice_2d", class(x))
  x
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
