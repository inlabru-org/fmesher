#' @include deprecated.R

# fm_bbox ####


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
#' precomputed interval information, stored as a list of 2-vector ranges,
#' `list(xlim, ylim, ...)`.
#' @export
fm_bbox.list <- function(x, ...) {
  if (!all(vapply(x, function(xx) {
    length(xx) == 2
  }, TRUE))) {
    stop("List not coercible to fm_bbox; some list element has length != 2.")
  }
  structure(
    x,
    class = "fm_bbox"
  )
}

#' @rdname fm_bbox
#' @usage
#' ## S3 method for class 'NULL'
#' fm_bbox(...)
#' @export
fm_bbox.NULL <- function(...) {
  fm_bbox(list())
}

#' @rdname fm_bbox
#' @export
fm_bbox.numeric <- function(x, ...) {
  x <- as.vector(x)
  fm_bbox(
    list(
      if (all(is.na(x))) {
        c(NA_real_, NA_real_)
      } else {
        range(x, na.rm = TRUE)
      }
    )
  )
}

#' @rdname fm_bbox
#' @export
fm_bbox.matrix <- function(x, ...) {
  do.call(
    c,
    c(
      lapply(
        seq_len(ncol(x)),
        function(k) {
          fm_bbox(x[, k, drop = TRUE])
        }
      ),
      list(.join = TRUE)
    )
  )
}

#' @rdname fm_bbox
#' @export
fm_bbox.Matrix <- function(x, ...) {
  do.call(
    c,
    c(
      lapply(
        seq_len(ncol(x)),
        function(k) {
          fm_bbox(x[, k, drop = TRUE])
        }
      ),
      list(.join = TRUE)
    )
  )
}

#' @rdname fm_bbox
#' @export
fm_bbox.fm_bbox <- function(x, ...) {
  x
}

#' @rdname fm_bbox
#' @export
fm_bbox.fm_mesh_1d <- function(x, ...) {
  fm_bbox(x[["interval"]])
}

#' @rdname fm_bbox
#' @export
fm_bbox.fm_mesh_2d <- function(x, ...) {
  if (fm_manifold(x, "R2")) {
    d <- 2L
  } else {
    d <- min(NCOL(x[["loc"]]), 3L)
  }
  box <- fm_bbox(x[["loc"]][, seq_len(d), drop = FALSE])
  ranges <- vapply(
    seq_along(box),
    function(k) {
      diff(box[[k]])
    },
    numeric(1)
  )
  if (all(ranges <= 0)) {
    return(box)
  }
  d <- max(1L, max(which(ranges > 0)))
  box <- box[seq_len(d)]
  box
}

#' @rdname fm_bbox
#' @export
fm_bbox.fm_mesh_3d <- function(x, ...) {
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
fm_bbox.fm_lattice_Nd <- function(x, ...) {
  fm_bbox(x[["loc"]])
}

#' @rdname fm_bbox
#' @export
fm_bbox.fm_tensor <- function(x, ...) {
  do.call(c, c(lapply(x[["fun_spaces"]], fm_bbox), list(.join = TRUE)))
}

#' @rdname fm_bbox
#' @export
fm_bbox.fm_collect <- function(x, ...) {
  do.call(c, c(lapply(x[["fun_spaces"]], fm_bbox), list(.join = FALSE)))
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
fm_as_bbox <- function(x, ...) {
  fm_bbox(x, ...)
}


#' @export
#' @param x `fm_bbox` object from which to extract element(s)
#' @param i indices specifying elements to extract
#' @describeIn fm_bbox Extract sub-list
`[.fm_bbox` <- function(x, i) {
  object <- NextMethod()
  class(object) <- class(x)
  object
}

#' @export
#' @describeIn fm_bbox The `...` arguments should be `fm_bbox` objects, or
#'   coercible with `fm_as_bbox(list(...))`.
#' @param .join logical; if `TRUE`, concatenate the bounding boxes into a single
#'   multi-dimensional bounding box. Default is `FALSE`.
#' @returns For [c.fm_bbox()], a `fm_bbox_list` object if `.join = FALSE` (the
#'   default) or an `fm_bbox` object if `.join = TRUE`.
#' @examples
#' m <- c(A = fm_bbox(cbind(1, 2)), B = fm_bbox(cbind(3, 4)))
#' str(m)
#' str(m[2])
`c.fm_bbox` <- function(..., .join = FALSE) {
  y <- lapply(list(...), function(xx) fm_as_list(xx, .class_stub = "bbox"))
  y <- do.call("c", y)
  if (.join) {
    y <- unclass(y)
    y <- lapply(y, function(xx) unclass(xx))
    y <- do.call("c", y)
    y <- fm_bbox(y)
  }
  y
}

#' @describeIn fm_bbox Convert a list to a `fm_bbox_list` object, with
#' each element converted to an `fm_bbox` object.
#' @export
#' @examples
#' m <- fm_as_bbox_list(list(
#'   A = fm_bbox(cbind(1, 2)),
#'   B = fm_bbox(cbind(3, 4))
#' ))
#' str(fm_as_bbox_list(m))
fm_as_bbox_list <- function(x, ...) {
  if (is.list(x) && !inherits(x, c("fm_bbox", "fm_bbox_list"))) {
    y <- lapply(x, function(xx) fm_as_bbox(xx))
    y <- do.call("c", y)
  } else {
    y <- fm_as_list(x, ..., .class_stub = "bbox")
  }
  y
}
