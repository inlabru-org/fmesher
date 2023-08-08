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
#' @describeIn fm_bbox The `...` arguments should be `fm_bbox`
#' objects, or coercible with `fm_as_bbox(list(...))`.
#' @returns A `fm_bbox_list` object
#' @examples
#' m <- c(A = fm_bbox(cbind(1, 2), B = fm_bbox(cbind(3, 4))))
#' str(m)
#' str(m[2])
`c.fm_bbox` <- function(...) {
  y <- lapply(list(...), function(xx) fm_as_list(xx, .class_stub = "bbox"))
  return(do.call("c", y))
}
