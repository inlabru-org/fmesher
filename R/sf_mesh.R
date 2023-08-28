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
#'
fm_as_sfc <- function(x, ...) {
  UseMethod("fm_as_sfc")
}

#' @describeIn fm_as_sfc `r lifecycle::badge("experimental")`
#'
#' @param multi logical; if `TRUE`, attempt to a `sfc_MULTIPOLYGON`, otherwise
#' a set of `sfc_POLYGON`. Default `FALSE`
#' @returns * `fm_as_sfc`: An `sfc_MULTIPOLYGON` or `sfc_POLYGON` object
#' @exportS3Method fm_as_sfc inla.mesh
#' @export
fm_as_sfc.inla.mesh <- function(x, ..., multi = FALSE) {
  fm_as_sfc.fm_mesh_2d(fm_as_mesh_2d(x), ..., multi = multi)
}

#' @describeIn fm_as_sfc `r lifecycle::badge("experimental")`
#'
#' @param multi logical; if `TRUE`, attempt to a `sfc_MULTIPOLYGON`, otherwise
#' a set of `sfc_POLYGON`. Default `FALSE`
#' @returns * `fm_as_sfc`: An `sfc_MULTIPOLYGON` or `sfc_POLYGON` object
#' @exportS3Method fm_as_sfc inla.mesh
#' @export
fm_as_sfc.fm_mesh_2d <- function(x, ..., multi = FALSE) {
  stopifnot(inherits(x, "fm_mesh_2d"))
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
      check_ring_dir = TRUE
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
      )
    )
  }
  sf::st_crs(geom) <- fm_crs(x$crs)
  geom
}


#' @describeIn fm_as_sfc `r lifecycle::badge("experimental")`
#'
#' @exportS3Method fm_as_sfc inla.mesh.segment
#' @export
fm_as_sfc.inla.mesh.segment <- function(x, ..., multi = FALSE) {
  fm_as_sfc.fm_segm(fm_as_segm(x), ..., multi = multi)
}
#' @describeIn fm_as_sfc `r lifecycle::badge("experimental")`
#'
#' @export
fm_as_sfc.fm_segm <- function(x, ..., multi = FALSE) {
  stopifnot(inherits(x, "fm_segm"))

  if (all(fm_is_bnd(x))) {
    warning("fm_as_sfc currently only supports (multi)linestring output")
  }

  group_segments <- list()
  used_seg <- c()
  active_group <- 0L
  group <- integer(nrow(x$idx))
  closed_loop <- logical(0)
  while (any(group == 0L)) {
    active_group <- active_group + 1L
    closed_loop <- c(closed_loop, FALSE)
    curr_seg <- which.min(group)
    group[curr_seg] <- active_group
    group_segments[[active_group]] <- curr_seg
    used_seg <- c(used_seg, curr_seg)
    repeat {
      next_seg <- which(x$idx[, 1] == x$idx[curr_seg, 2])
      if (length(next_seg) == 0) {
        break
      }
      if (any(next_seg %in% used_seg)) {
        closed_loop[active_group] <- TRUE
        break
      }
      curr_seg <- min(next_seg)
      group[curr_seg] <- active_group
      group_segments[[active_group]] <-
        c(group_segments[[active_group]], curr_seg)
      used_seg <- c(used_seg, curr_seg)
    }
  }

  if (multi) {
    geom <- sf::st_sfc(
      sf::st_multilinestring(
        lapply(
          seq_along(group_segments),
          function(k) {
            x$loc[
              c(
                x$idx[group_segments[[k]], 1],
                x$idx[group_segments[[k]][length(group_segments[[k]])], 2]
              ), ,
              drop = FALSE
            ]
          }
        ),
        dim = "XYZ"
      )
    )
  } else {
    geom <- sf::st_sfc(
      lapply(
        seq_along(group_segments),
        function(k) {
          sf::st_linestring(
            x$loc[
              c(
                x$idx[group_segments[[k]], 1],
                x$idx[group_segments[[k]][length(group_segments[[k]])], 2]
              ), ,
              drop = FALSE
            ],
            dim = "XYZ"
          )
        }
      )
    )
  }

  sf::st_crs(geom) <- fm_crs(x$crs)
  geom
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
    warning("More than one MULTIPOLYGON detected, but conversion method only uses one.",
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
      grp <- seq_len(length(sfc))
    } else {
      grp <- c(grp, rep(grp[length(grp)], length(sfc) - length(grp)))
    }
    for (k in seq_len(length(sfc))) {
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
      grp <- seq_len(length(sfc))
    } else {
      grp <- c(grp, rep(grp[length(grp)], length(sfc) - length(grp)))
    }
    for (k in seq_len(length(sfc))) {
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
            idx <- seq_len(length(subset))
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
      grp <- seq_len(length(sfc))
    } else {
      grp <- c(grp, rep(grp[length(grp)], length(sfc) - length(grp)))
    }
    for (k in seq_len(length(sfc))) {
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
      grp <- seq_len(length(sfc))
    } else {
      grp <- c(grp, rep(grp[length(grp)], length(sfc) - length(grp)))
    }
    for (k in seq_len(length(sfc))) {
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
      grp <- seq_len(length(x))
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
