#' @include deprecated.R

# fm ####

#' @title Convert objects to fmesher objects
#' @description
#' Used for conversion from general objects
#' (usually `inla.mesh` and other INLA specific classes)
#' to `fmesher` classes.
#'
#' @param x Object to be converted
#' @param ... Arguments forwarded to submethods
#' @rdname fm_as_fm
#' @returns An object of some `fm_*` class
#' @export
#' @family object creation and conversion
fm_as_fm <- function(x, ...) {
  UseMethod("fm_as_fm")
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
fm_as_fm.fm_bbox <- function(x, ...) {
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
