#' @include deprecated.R

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
  fm_as_list(x, ..., .class_stub = "tensor")
}
#' @rdname fm_as_tensor
#' @param x Object to be converted
#' @export
fm_as_tensor.fm_tensor <- function(x, ...) {
  #  class(x) <- c("fm_tensor", setdiff(class(x), "fm_tensor"))
  x
}
