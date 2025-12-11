#' @include mesh.R
#' @include deprecated.R


# fm_as_sfc ####

#' @title Conversion methods from mesh related objects to sfc
#' @rdname fm_as_sfc
#' @family fm_as
#' @param x An object to be coerced/transformed/converted into another class
#' @param ... Arguments passed on to other methods
#' @export
#' @family object creation and conversion
#' @examples
#' fm_as_sfc(fmexample$mesh)
#' fm_as_sfc(fmexample$mesh, multi = TRUE)
#' fm_as_sfc(fmexample$mesh, format = "loc")
#'
#' # Boundary edge conversion to polygons is supported from version 0.4.0.9002:
#' fm_as_sfc(fmexample$mesh, format = "bnd")
#'
fm_as_sfc <- function(x, ...) {
  UseMethod("fm_as_sfc")
}

#' @describeIn fm_as_sfc `r lifecycle::badge("experimental")`
#'
#' @param format One of "mesh", "int", "bnd", or "loc". Default
#'   "mesh".
#' @param multi logical; if `TRUE`, attempt to a
#'   `sfc_MULTIPOLYGON/LINESTRING/POINT/GEOMETRYCOLLECTION`, otherwise a set of
#'   `sfc_POLYGON/LINESTRING/POINT`. Default `FALSE`
#' @returns * `fm_as_sfc`: An
#' `sfc_MULTIPOLYGON/LINESTRING/POINT/GEOMETRYCOLLECTION` or
#' `sfc_POLYGON/LINESTRING/POINT` object
#' @export
fm_as_sfc.fm_mesh_2d <- function(x,
                                 ...,
                                 format = NULL,
                                 multi = FALSE) {
  stopifnot(inherits(x, "fm_mesh_2d"))
  format <- match.arg(format, c("mesh", "int", "bnd", "loc"))
  if (identical(format, "mesh")) {
    if (multi) {
      geom <- sf::st_sfc(
        sf::st_multipolygon(
          lapply(
            seq_len(nrow(x$graph$tv)),
            function(k) {
              list(x$loc[x$graph$tv[k, c(1, 2, 3, 1)], , drop = FALSE])
            }
          ),
          dim = "XYZ"
        ),
        check_ring_dir = TRUE,
        crs = fm_crs(x$crs)
      )
    } else {
      geom <- sf::st_sfc(
        lapply(
          seq_len(nrow(x$graph$tv)),
          function(k) {
            sf::st_polygon(
              list(x$loc[x$graph$tv[k, c(1, 2, 3, 1)], , drop = FALSE]),
              dim = "XYZ"
            )
          }
        ),
        crs = fm_crs(x$crs)
      )
    }
  } else if (format %in% c("int", "bnd")) {
    geom <- fm_as_sfc(
      fm_segm(x, boundary = identical(format, "bnd")),
      multi = multi
    )
  } else if (identical(format, "loc")) {
    geom <- sf::st_sfc(
      sf::st_multipoint(x$loc, dim = "XYZ"),
      crs = fm_crs(x$crs)
    )
    if (!multi) {
      geom <- sf::st_geometry(sf::st_cast(sf::st_sf(geometry = geom), "POINT"))
    }
  } else {
    stop("Unsupported mesh conversion format '", format, "'.")
  }
  geom
}


#' @describeIn fm_as_sfc `r lifecycle::badge("experimental")`
#'
#' @export
fm_as_sfc.fm_segm <- function(x, ..., multi = FALSE) {
  stopifnot(inherits(x, "fm_segm"))

  segm_comp <- fm_components(x)

  segm_bnd <- vapply(
    seq_along(segm_comp),
    function(k) {
      all(fm_is_bnd(segm_comp[[k]]))
    },
    logical(1L)
  )
  segm_bnd <- segm_comp[segm_bnd]
  segm_int <- vapply(
    seq_along(segm_comp),
    function(k) {
      all(!fm_is_bnd(segm_comp[[k]]))
    },
    logical(1L)
  )
  segm_int <- segm_comp[segm_int]

  geom_bnd <- NULL
  if (length(segm_bnd) > 0) {
    geom_bnd <-
      lapply(
        seq_along(segm_bnd),
        function(k) {
          # If all segments are boundary segments, we can use a polygon
          # with the last point repeated to close the ring.
          # For each positive area polygon, subtract negative area polygons
          # that are fully contained in the positive area polygon.
          sf::st_polygon(
            list(
              segm_bnd[[k]]$loc[
                c(
                  segm_bnd[[k]]$idx[, 1],
                  segm_bnd[[k]]$idx[nrow(segm_bnd[[k]]$idx), 2]
                ), ,
                drop = FALSE
              ]
            ),
            dim = "XYZ"
          )
        }
      )
    areas <- fm_area(segm_bnd) # fm_area.fm_segm_list
    geom_bnd_pos <- geom_bnd[areas > 0]
    geom_bnd_neg <- geom_bnd[areas < 0]
    contains <- sf::st_contains(
      sf::st_sfc(geom_bnd_pos, crs = fm_crs(x)),
      sf::st_sfc(geom_bnd_neg, crs = fm_crs(x))
    )
    geom_list <- lapply(
      seq_along(geom_bnd_pos),
      function(k) {
        if (length(contains[[k]]) > 0) {
          poly <- sf::st_difference(
            geom_bnd_pos[[k]],
            sf::st_sfc(geom_bnd_neg[contains[[k]]])
          )
        } else {
          poly <- geom_bnd_pos[[k]]
        }
        poly
      }
    )
    geom_bnd <- sf::st_sfc(geom_list, crs = fm_crs(x))
  } else {
    geom_bnd <- sf::st_sfc(crs = fm_crs(x))
  }

  geom_int <- NULL
  if (length(segm_int) > 0) {
    geom_int <- sf::st_sfc(
      lapply(
        seq_along(segm_int),
        function(k) {
          sf::st_linestring(
            segm_int[[k]]$loc[
              c(
                segm_int[[k]]$idx[, 1],
                segm_int[[k]]$idx[nrow(segm_int[[k]]$idx), 2]
              ), ,
              drop = FALSE
            ],
            dim = "XYZ"
          )
        }
      ),
      crs = fm_crs(x)
    )
  }

  if (length(geom_bnd) > 0 && length(geom_int) > 0) {
    geom <- sf::st_sfc(
      c(geom_bnd, geom_int),
      crs = fm_crs(x)
    )
  } else if (length(geom_bnd) > 0) {
    geom <- geom_bnd
  } else if (length(geom_int) > 0) {
    geom <- geom_int
  } else {
    geom <- sf::st_sfc(crs = fm_crs(x))
  }

  if (multi) {
    geom <- sf::st_union(geom)
  }

  geom
}

#' @rdname fm_as_sfc
#'
#' @export
fm_as_sfc.fm_segm_list <- function(x, ...) {
  do.call(c, lapply(x, fm_as_sfc, ...))
}

#' @rdname fm_as_sfc
#'
#' @export
fm_as_sfc.sfc <- function(x, ...) {
  x
}

#' @rdname fm_as_sfc
#' @export
fm_as_sfc.sf <- function(x, ...) {
  sf::st_geometry(x)
}


# fm_as_mesh_2d ####

#' @rdname fm_as_mesh_2d
#' @export
fm_as_mesh_2d.sfg <-
  function(x, ...) {
    fm_as_mesh_2d(sf::st_sfc(x), ...)
  }

#' @rdname fm_as_mesh_2d
#'
#' @export
fm_as_mesh_2d.sfc_MULTIPOLYGON <- function(x, ...) {
  if (length(x) > 1) {
    warning(
      "More than one MULTIPOLYGON detected,",
      " but conversion method only uses one.",
      immediate. = TRUE
    )
  }
  # Ensure correct CCW ring orientation; sf doesn't take into account
  # that geos has CW as canonical orientation
  x <- sf::st_sfc(x, check_ring_dir = TRUE)
  tv <- matrix(seq_len(3 * length(x[[1]])), length(x[[1]]), 3, byrow = TRUE)
  loc <- do.call(
    rbind,
    lapply(
      x[[1]],
      function(xx) {
        if ((length(xx) > 1) ||
          (nrow(xx[[1]]) > 4)) {
          stop("Invalid geometry; non-triangle detected.")
        }
        xx[[1]][1:3, , drop = FALSE]
      }
    )
  )
  crs <- fm_CRS(sf::st_crs(x))
  mesh <- fm_rcdt_2d_inla(
    loc = loc,
    tv = tv,
    ...,
    crs = crs
  )
  mesh
}

#' @rdname fm_as_mesh_2d
#'
#' @export
fm_as_mesh_2d.sfc_POLYGON <- function(x, ...) {
  # Ensure correct CCW ring orientation; sf doesn't take into account
  # that geos has CW as canonical orientation
  sfc <- sf::st_sfc(x, check_ring_dir = TRUE)
  tv <- matrix(seq_len(3 * NROW(x)), NROW(x), 3, byrow = TRUE)
  loc <- do.call(
    rbind,
    lapply(
      x,
      function(xx) {
        if ((length(xx) > 1) ||
          (nrow(xx[[1]]) > 4)) {
          stop("Invalid geometry; non-triangle detected.")
        }
        xx[[1]][1:3, , drop = FALSE]
      }
    )
  )
  crs <- fm_CRS(sf::st_crs(x))
  mesh <- fm_rcdt_2d_inla(
    loc = loc,
    tv = tv,
    ...,
    crs = crs
  )
  mesh
}

#' @rdname fm_as_mesh_2d
#' @export
fm_as_mesh_2d.sf <-
  function(x, ...) {
    fm_as_mesh_2d(sf::st_geometry(x), ...)
  }


# fm_as_segm ####

#' @rdname fm_as_segm
#' @export
fm_as_segm.sfg <-
  function(x, ...) {
    fm_as_segm(sf::st_sfc(x), ...)
  }

#' @rdname fm_as_segm
#' @param reverse logical; When TRUE, reverse the order of the input points.
#'   Default `FALSE`
#' @param grp if non-null, should be an integer vector of grouping labels for
#'   one for each segment.
#'    Default `NULL`
#' @param is.bnd logical; if `TRUE`, set the boundary flag for the segments.
#'   Default `TRUE`
#' @export
fm_as_segm.sfc_POINT <-
  function(x, reverse = FALSE, grp = NULL, is.bnd = TRUE, ...) {
    sfc <- x
    crs <- sf::st_crs(sfc)

    loc <- sf::st_coordinates(sfc)
    coord_names <- intersect(c("X", "Y", "Z"), colnames(loc))
    loc <- unname(loc[, coord_names, drop = FALSE])

    n <- dim(loc)[1L]
    if (all(is.bnd)) {
      idx <- c(seq_len(n), 1L)
    } else {
      idx <- seq_len(n)
    }
    if (reverse) {
      idx <- rev(idx)
      if (!is.null(grp)) {
        grp <- rev(grp)
      }
    }
    fm_segm(
      loc = loc, idx = idx, grp = grp, is.bnd = all(is.bnd),
      crs = crs
    )
  }

#' @rdname fm_as_segm
#' @param join logical; if `TRUE`, join input segments with common vertices.
#'    Default `TRUE`
#' @export
#' @examples
#' (segm <- fm_segm(fmexample$mesh, boundary = FALSE))
#' (segm_sfc <- fm_as_sfc(segm))
#' (fm_as_segm(segm_sfc))
#'
fm_as_segm.sfc_LINESTRING <-
  function(x, join = TRUE, grp = NULL, reverse = FALSE, ...) {
    sfc <- x

    crs <- sf::st_crs(sfc)

    segm <- list()
    if (is.null(grp)) {
      grp <- seq_along(sfc)
    } else {
      grp <- c(grp, rep(grp[length(grp)], length(sfc) - length(grp)))
    }
    for (k in seq_along(sfc)) {
      loc <- sf::st_coordinates(sfc[k])
      coord_names <- intersect(c("X", "Y", "Z"), colnames(loc))
      if (nrow(loc) == 0) {
        segm[[k]] <- fm_segm(
          loc = matrix(0, 0, length(coord_names)),
          idx = matrix(0L, 0, 2),
          grp = integer(0),
          is.bnd = FALSE,
          crs = crs
        )
        next
      }
      loc <- unname(loc[, coord_names, drop = FALSE])

      n <- dim(loc)[1L]
      if (reverse) {
        idx <- seq(n, 1L, length.out = n)
      } else {
        idx <- seq_len(n)
      }
      segm[[k]] <- fm_segm(
        loc = loc,
        idx = idx,
        grp = grp[k],
        is.bnd = FALSE,
        crs = crs
      )
    }

    if (join) {
      segm <- fm_segm_join(segm)
    }
    segm
  }

#' @rdname fm_as_segm
#' @param join logical; if `TRUE`, join input segments with common vertices.
#'    Default `TRUE`
#' @export
fm_as_segm.sfc_MULTILINESTRING <-
  function(x, join = TRUE, grp = NULL, reverse = FALSE, ...) {
    sfc <- x

    crs <- sf::st_crs(sfc)

    segm <- list()
    if (is.null(grp)) {
      grp <- seq_along(sfc)
    } else {
      grp <- c(grp, rep(grp[length(grp)], length(sfc) - length(grp)))
    }
    for (k in seq_along(sfc)) {
      loc <- sf::st_coordinates(sfc[k])
      coord_names <- intersect(c("X", "Y", "Z"), colnames(loc))
      if (nrow(loc) == 0) {
        segm[[k]] <- fm_segm(
          loc = matrix(0, 0, length(coord_names)),
          idx = matrix(0L, 0, 2),
          grp = integer(0),
          is.bnd = FALSE,
          crs = crs
        )
        next
      }
      coord_names <- intersect(c("X", "Y", "Z"), colnames(loc))
      Linfo <- loc[, c("L1", "L2"), drop = FALSE]
      uniqueLinfo <- unique(Linfo)
      loc <- unname(loc[, coord_names, drop = FALSE])

      segm_k <-
        lapply(
          seq_len(nrow(uniqueLinfo)),
          function(i) {
            subset <- which((Linfo[, 1] == uniqueLinfo[i, 1]) &
              (Linfo[, 2] == uniqueLinfo[i, 2]))
            idx <- seq_along(subset)
            if (reverse) {
              idx <- rev(idx)
            }
            fm_segm(
              loc = loc[subset, , drop = FALSE],
              idx = idx,
              grp = grp[k],
              is.bnd = FALSE,
              crs = crs
            )
          }
        )
      segm[[k]] <- fm_segm_join(segm_k)
    }

    if (join) {
      segm <- fm_segm_join(segm)
    }
    segm
  }

#' @rdname fm_as_segm
#' @export
fm_as_segm.sfc_POLYGON <-
  function(x, join = TRUE, grp = NULL, ...) {
    # Ensure correct CCW ring orientation; sf doesn't take into account
    # that geos has CW as canonical orientation
    sfc <- sf::st_sfc(x, check_ring_dir = TRUE)
    crs <- sf::st_crs(sfc)

    segm <- list()
    if (is.null(grp)) {
      grp <- seq_along(sfc)
    } else {
      grp <- c(grp, rep(grp[length(grp)], length(sfc) - length(grp)))
    }
    for (k in seq_along(sfc)) {
      loc <- sf::st_coordinates(sfc[k])
      coord_names <- intersect(c("X", "Y", "Z"), colnames(loc))
      if (nrow(loc) == 0) {
        segm[[k]] <- fm_segm(
          loc = matrix(0, 0, length(coord_names)),
          idx = matrix(0L, 0, 2),
          grp = integer(0),
          is.bnd = TRUE,
          crs = crs
        )
        next
      }
      coord_names <- intersect(c("X", "Y", "Z"), colnames(loc))
      L1info <- loc[, "L1", drop = TRUE]
      L2info <- loc[, "L2", drop = TRUE]
      loc <- unname(loc[, coord_names, drop = FALSE])
      # If winding directions are correct, all info is already available
      # For 3D, cannot check winding, so must assume correct.
      segm_k <-
        lapply(
          unique(L1info),
          function(i) {
            subset <- which(L1info == i)
            # sfc_POLYGON repeats the initial point within each L1
            n <- length(subset) - 1
            subset <- subset[-(n + 1)]
            idx <- c(seq_len(n), 1L)
            fm_segm(
              loc = loc[subset, , drop = FALSE],
              idx = idx,
              grp = grp[k],
              is.bnd = TRUE,
              crs = crs
            )
          }
        )
      segm[[k]] <- fm_segm_join(segm_k)
    }

    if (join) {
      segm <- fm_segm_join(segm)
    }
    segm
  }

#' @rdname fm_as_segm
#' @export
fm_as_segm.sfc_MULTIPOLYGON <-
  function(x, join = TRUE, grp = NULL, ...) {
    # Ensure correct CCW ring orientation; sf doesn't take into account
    # that geos has CW as canonical orientation
    sfc <- sf::st_sfc(x, check_ring_dir = TRUE)
    crs <- sf::st_crs(sfc)

    segm <- list()
    if (is.null(grp)) {
      grp <- seq_along(sfc)
    } else {
      grp <- c(grp, rep(grp[length(grp)], length(sfc) - length(grp)))
    }
    for (k in seq_along(sfc)) {
      loc <- sf::st_coordinates(sfc[k])
      coord_names <- intersect(c("X", "Y", "Z"), colnames(loc))
      if (nrow(loc) == 0) {
        segm[[k]] <- fm_segm(
          loc = matrix(0, 0, length(coord_names)),
          idx = matrix(0L, 0, 2),
          grp = integer(0),
          is.bnd = TRUE,
          crs = crs
        )
        next
      }
      coord_names <- intersect(c("X", "Y", "Z"), colnames(loc))
      Linfo <- loc[, c("L1", "L2"), drop = FALSE]
      loc <- unname(loc[, coord_names, drop = FALSE])
      # If winding directions are correct, all info is already available
      # For 3D, cannot check winding, so must assume correct.
      uniqueLinfo <- unique(Linfo)
      segm_k <-
        lapply(
          seq_len(nrow(uniqueLinfo)),
          function(i) {
            subset <- which((Linfo[, 1] == uniqueLinfo[i, 1]) &
              (Linfo[, 2] == uniqueLinfo[i, 2]))
            # sfc_POLYGON repeats the initial point
            n <- length(subset) - 1
            subset <- subset[-(n + 1)]
            idx <- c(seq_len(n), 1L)
            fm_segm(
              loc = loc[subset, , drop = FALSE],
              idx = idx,
              grp = grp[k],
              is.bnd = TRUE,
              crs = crs
            )
          }
        )
      segm[[k]] <- fm_segm_join(segm_k)
    }

    if (join) {
      segm <- fm_segm_join(segm)
    }
    segm
  }

#' @rdname fm_as_segm
#' @export
fm_as_segm.sfc_GEOMETRY <-
  function(x, grp = NULL, join = TRUE, ...) {
    if (is.null(grp)) {
      grp <- seq_along(x)
    } else {
      grp <- c(grp, rep(grp[length(grp)], length(x) - length(grp)))
    }
    segm <-
      lapply(
        seq_along(x),
        function(k) {
          fm_as_segm(x[k], grp = grp[k], join = join, ...)
        }
      )
    if (join) {
      segm <- fm_segm_join(segm)
    }
    segm
  }

#' @rdname fm_as_segm
#' @export
fm_as_segm.sf <-
  function(x, ...) {
    sfc <- sf::st_geometry(x)
    fm_as_segm(sfc, ...)
  }
