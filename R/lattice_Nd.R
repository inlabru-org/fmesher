#' @include deprecated.R

# fm_lattice_Nd ####

#' @title Lattice grids for N dimensions
#'
#' @description
#' Construct an N-dimensional lattice grid
#'
#' @export
#' @param ... Passed on to submethods
#' @family object creation and conversion
#' @rdname fm_lattice_Nd
fm_lattice_Nd <- function(x = NULL, ...) {
  # Need to specify the dispatch object explicitly to handle the NULL case:
  UseMethod("fm_lattice_Nd", x)
}

#' @param x `list`, `data.frame`, `matrix`, `fm_bbox` or `NULL`. If a list of
#'   vectors, `as.matrix(expand.grid(x))` is used to create a full grid
#'   coordinates. `data.frame` and `matrix` input is assumed to follow the same
#'   ordering convention as the output of `expand.grid()`. of length N of
#'   vectors or grid matrices of coordinate values. List vector values are
#'   sorted before use.
#' @param dims numeric; the size of the grid of dimension `length(dims)`
#' @param values list of grid axis values
#' @returns An `fm_lattice_Nd` object with elements
#' \describe{
#' \item{dims}{integer vector}
#' \item{values}{the grid coordinate axis values}
#' \item{loc}{matrix of constructed grid coordinates}
#' }
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @seealso [fm_mesh_3d()]
#' @examples
#' (lattice <- fm_lattice_Nd(
#'   list(
#'     seq(0, 1, length.out = 3),
#'     seq(0, 1, length.out = 4),
#'     seq(0, 1, length.out = 2)
#'   )
#' ))
#'
#' if (requireNamespace("geometry", quietly = TRUE)) {
#'   (mesh <- fm_delaunay_3d(lattice$loc))
#' }
#' @rdname fm_lattice_Nd
#' @export
fm_lattice_Nd.matrix <- function(
    x = NULL,
    dims = NULL,
    values = NULL,
    ...) {
  if (is.null(dims)) {
    stop("The 'dims' argument must be specified.")
  }
  if (!(nrow(x) == prod(dims))) {
    stop("The number of rows in 'x' must match the product of 'dims'.")
  }
  if (is.null(values)) {
    values <- lapply(seq_len(ncol(x)), function(k) {
      sort(unique(x[, k]))
    })
    names(values) <- colnames(x)
    the_dims <- lengths(values)
    if (!identical(dims, the_dims)) {
      warning(paste0(
        "The number of unique values in each dimension doesn't ",
        "match the 'dims' argument."
      ))
    }
  } else {
    the_dims <- lengths(values)
    if (!identical(dims, the_dims)) {
      stop(paste0(
        "The number of values in each dimension doesn't ",
        "match the 'dims' argument."
      ))
    }
  }
  fm_lattice_Nd_create(loc = x, dims = dims, values = values)
}


#' @rdname fm_lattice_Nd
#' @export
fm_lattice_Nd.data.frame <- function(
    x = NULL,
    ...) {
  fm_lattice_Nd(as.matrix(x), ...)
}

#' @rdname fm_lattice_Nd
#' @export
fm_lattice_Nd.list <- function(
    x = NULL,
    dims = NULL,
    ...) {
  mat <- vapply(x, function(xx) is.array(xx), logical(1))
  if (all(mat)) {
    if (is.null(dims)) {
      dims <- dim(x[[1]])
    }
    if (!all(vapply(x, function(xx) identical(dim(xx), dims), logical(1)))) {
      stop(
        "The x elements must have identical array dimensions and match dims."
      )
    }
    loc <- do.call(rbind, lapply(x, as.vector))
    values <- lapply(seq_len(ncol(loc)), function(k) sort(unique(loc[, k])))
    names(values) <- names(x)
  } else {
    if (!all(!mat)) {
      stop("The x elements must all be arrays or all be vectors.")
    }
    if (is.null(dims)) {
      dims <- vapply(x, function(xx) length(xx), integer(1))
    } else {
      the_dims <- vapply(x, function(xx) length(xx), integer(1))
      for (k in seq_along(x)) {
        if (dims[k] != the_dims[k]) {
          if (!(the_dims[k] %in% c(1L, 2L))) {
            stop("The x element lengths must match dims or be 1 or 2.")
          }
          x[[k]] <- seq(x[[k]][1], x[[k]][2], length.out = dims[k])
        }
      }
    }
    values <- lapply(x, function(xx) sort(unique(xx)))
    loc <- expand.grid(values, stringsAsFactors = FALSE)
  }

  fm_lattice_Nd_create(loc = loc, dims = dims, values = values)
}

#' @rdname fm_lattice_Nd
#' @export
fm_lattice_Nd.fm_bbox <- function(
    x = NULL,
    dims = NULL,
    ...) {
  if (is.null(dims)) {
    dims <- rep(2L, length(x))
  }
  values <- lapply(seq_len(length(x)), function(k) {
    seq(x[[k]][1], x[[k]][2], length.out = dims[k])
  })
  loc <- as.matrix(expand.grid(values, stringsAsFactors = FALSE))

  fm_lattice_Nd_create(loc = loc, dims = dims, values = values)
}

#' @describeIn fm_lattice_Nd Ignores the `NULL` `x` and creates a lattice
#' based on `values` (if non-NULL) and `dims` unit hypercube
#' lattice grid with `dims` dimensions.
#' @export
fm_lattice_Nd.NULL <- function(x = NULL, ..., dims = NULL) {
  if (is.null(dims)) {
    dims <- c(2L, 2L)
  }
  fm_lattice_Nd(lapply(
    dims,
    function(k) {
      seq(0.0, 1.0, length.out = k)
    }
  ))
}

fm_lattice_Nd_create <- function(loc, dims, values) {
  structure(
    list(
      loc = as.matrix(loc),
      dims = dims,
      values = values,
      manifold = paste0("R", ncol(loc))
    ),
    class = "fm_lattice_Nd"
  )
}


#' @title Convert objects to `fm_lattice_Nd`
#' @describeIn fm_as_lattice_Nd Convert an object to `fm_lattice_Nd`.
#' @param x Object to be converted.
#' @param ... Arguments passed on to submethods
#' @returns An `fm_lattice_Md` or `fm_lattice_Nd_list` object
#' @export
#' @family object creation and conversion
#' @export
#' @examples
#' (fm_as_lattice_Nd_list(list(
#'   fm_lattice_Nd(list(1:3, 1:2)),
#'   fm_lattice_Nd(list(1:4))
#' )))
#'
fm_as_lattice_Nd <- function(...) {
  UseMethod("fm_as_lattice_Nd")
}
#' @describeIn fm_as_lattice_Nd Convert each element of a list
#' @export
fm_as_lattice_Nd_list <- function(x, ...) {
  fm_as_list(x, ..., .class_stub = "lattice_Nd")
}
#' @rdname fm_as_lattice_Nd
#' @param x Object to be converted
#' @export
fm_as_lattice_Nd.fm_lattice_Nd <- function(x, ...) {
  #  class(x) <- c("fm_lattice_Nd", setdiff(class(x), "fm_lattice_Nd"))
  x
}
