#' @include deprecated.R

# fm_segm ####

#' @title Make a spatial segment object
#' @describeIn fm_segm Create a new `fm_segm` object.
#' @export
#' @param ... Passed on to submethods
#' @returns An `fm_segm` or `fm_segm_list` object
#' @family object creation and conversion
#' @examples
#' fm_segm(rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1)), is.bnd = FALSE)
#' fm_segm(rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1)), is.bnd = TRUE)
#'
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
    loc <- fm_unify_coords(loc)
    if (is.null(idx)) {
      idx <- if (all(is.bnd)) {
        c(seq_len(NROW(loc)), 1)
      } else {
        seq_len(NROW(loc))
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
      (max(idx, na.rm = TRUE) > NROW(loc))) {
      warning(
        "Segment indices (max=", max(idx, na.rm = TRUE),
        ") exceed specified location list length (",
        NROW(loc), ")."
      )
    }
  }

  if (!is.null(grp)) {
    if (!is.vector(grp) && !is.matrix(grp)) {
      stop("'grp' must be a vector or a matrix")
    }
    grp <- as.vector(grp)
    if (length(grp) < NROW(idx)) {
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
    idx.new <- rep(0L, NROW(loc))
    idx.new[as.vector(idx)] <- 1L
    loc <- loc[idx.new == 1L, , drop = FALSE]
    idx.new[idx.new == 1L] <- seq_len(sum(idx.new))
    idx <- matrix(idx.new[as.vector(idx)],
      nrow = nrow(idx),
      ncol = ncol(idx)
    )
  }

  if (length(is.bnd) == 1L) {
    is.bnd <- rep(is.bnd, nrow(idx))
  } else {
    if (any(is.bnd) && any(!is.bnd)) {
      warning("Mixed is.bnd status within a single fm_segm is not yet supported, replacing with FALSE.")
    }
    is.bnd <- rep(FALSE, nrow(idx))
  }

  ret <- structure(
    list(loc = loc, idx = idx, grp = grp, is.bnd = is.bnd, crs = crs),
    class = c("fm_segm", "inla.mesh.segment")
  )
  return(ret)
}

#' @describeIn fm_segm Join multiple `fm_segm` objects into a single `fm_segm`
#' object. If `is.bnd` is non-NULL, it overrides the input segment information.
#' Otherwise, it checks if the inputs are consistent.
#' @param grp When joining segments, use these group labels for segments
#' instead of the original group labels.
#' @param grp.default If `grp.default` is `NULL`, use these group labels for segments
#' with NULL group.
#' @export
fm_segm.fm_segm <- function(..., grp = NULL, grp.default = 0L, is.bnd = NULL) {
  segm <- fm_as_segm_list(list(...))
  fm_segm_join(segm, grp = grp, grp.default = grp.default, is.bnd = is.bnd)
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
#' object. If `is.bnd` is non-NULL, it overrides the segment information.
#' Otherwise it checks for consistency.
#' @export
#' @examples
#' fm_segm_join(fmexample$boundary_fm)
#'
fm_segm_join <- function(x, grp = NULL, grp.default = 0L, is.bnd = NULL) {
  segm <- fm_as_segm_list(x)
  segm <- lapply(seq_along(segm), function(k) {
    seg <- segm[[k]]
    if (!is.null(seg)) {
      if (is.null(grp)) {
        if (is.null(seg[["grp"]])) {
          seg[["grp"]] <-
            rep(
              grp.default[min(length(grp.default), k)],
              NROW(seg[["idx"]])
            )
        }
      } else {
        seg[["grp"]] <-
          rep(grp[min(length(grp), k)], NROW(seg[["idx"]]))
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

  Nloc <- vapply(segm, function(x) NROW(x$loc), 0L)
  cumNloc <- c(0, cumsum(Nloc))
  Nidx <- vapply(segm, function(x) NROW(x$idx), 0L)

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
      } else if (Nidx[k] > 0) {
        segm[[k]]$grp
      } else {
        integer(0)
      }
    }
  ))
  if (is.null(is.bnd)) {
    segm_ <- segm[vapply(segm, function(x) nrow(x[["idx"]]) > 0, TRUE)]
    if (length(segm_) == 0) {
      is.bnd <- TRUE
    } else {
      is.bnd <- vapply(segm_, function(x) all(fm_is_bnd(x)), TRUE)
      not.is.bnd <- vapply(segm_, function(x) all(!fm_is_bnd(x)), TRUE)
      if ((all(is.bnd) && all(!not.is.bnd)) ||
        (all(!is.bnd) && all(not.is.bnd))) {
        is.bnd <- all(is.bnd) && all(!not.is.bnd)
      } else {
        warning("Inconsistent 'is.bnd' attributes.  Setting 'is.bnd=FALSE'.")
        is.bnd <- FALSE
      }
    }
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
#' @describeIn fm_segm Split an `fm_segm` object by `grp` into an `fm_segm_list`
#' object, optionally keeping only some groups.
#' @export
fm_segm_split <- function(x, grp = NULL, grp.default = 0L) {
  if (is.null(x[["grp"]])) {
    x[["grp"]] <- rep(grp.default, NROW(x[["idx"]]))
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
  return(fm_as_segm_list(segm_list))
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
#' @examples
#' fm_segm(fmexample$mesh, boundary = TRUE)
#' fm_segm(fmexample$mesh, boundary = FALSE)
#'
fm_segm.fm_mesh_2d <- function(x, boundary = TRUE, grp = NULL, ...) {
  extract_segments <- function(mesh.loc,
                               segm,
                               grp = NULL,
                               is.bnd,
                               crs = NULL) {
    segments <- NULL
    if (NROW(segm[["idx"]]) == 0) {
      return(fm_segm(is.bnd = is.bnd, crs = crs))
    }
    if (is.null(grp)) {
      grp <- unique(sort(segm[["grp"]]))
    }
    extract <- (segm[["grp"]] %in% grp)
    if (!any(extract)) {
      return(fm_segm(is.bnd = is.bnd, crs = crs))
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


# fm_as_segm ####

#' @title Convert objects to `fm_segm`
#' @describeIn fm_as_segm Convert an object to `fm_segm`.
#' @param x Object to be converted.
#' @param ... Arguments passed on to submethods
#' @returns An `fm_segm` or `fm_segm_list` object
#' @export
#' @family object creation and conversion
#' @examples
#' fm_as_segm_list(list(
#'   fm_segm(fmexample$mesh),
#'   fm_segm(fmexample$mesh, boundary = FALSE)
#' ))
#'
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
  fm_as_list(x, ..., .class_stub = "segm")
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

# fm_is_bnd ####

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
  y <- lapply(list(...), fm_as_segm_list)
  return(do.call("c", y))
}

#' @export
#' @describeIn fm_segm_list The `...` arguments should be coercible to `fm_segm_list`
#' objects.
`c.fm_segm_list` <- function(...) {
  if (!all(vapply(
    list(...),
    function(xx) is.null(xx) || inherits(xx, "fm_segm_list"),
    TRUE
  ))) {
    y <- lapply(list(...), fm_as_segm_list)
    return(do.call("c", y))
  }
  object <- NextMethod()
  class(object) <- c("fm_segm_list", "fm_list", "list")
  object
}

#' @export
#' @param x `fm_segm_list` object from which to extract element(s)
#' @param i indices specifying elements to extract
#' @describeIn fm_segm_list Extract sub-list
`[.fm_segm_list` <- function(x, i) {
  object <- NextMethod()
  class(object) <- class(x)
  object
}
