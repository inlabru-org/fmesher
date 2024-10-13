#' @include deprecated.R

# fm_tensor ####

#' @title Make a tensor product function space
#' @description `r lifecycle::badge("experimental")`
#' Tensor product function spaces. The interface and object storage model
#' is experimental and may change.
#' @export
#' @param x list of function space objects, such as [fm_mesh_2d()].
#' @param ... Currently unused
#' @returns A `fm_tensor` or `fm_tensor_list` object. Elements of `fm_tensor`:
#' \describe{
#' \item{fun_spaces}{`fm_list` of function space objects}
#' \item{manifold}{character; manifold type summary. Regular subset of Rd "Rd",
#' if all function spaces have type "R",
#' torus connected "Td" if all function spaces have type "S", and otherwise "Md"
#' In all cases, `d` is the sum of the manifold dimensions of the function
#' spaces.}
#' }
#' @family object creation and conversion
#' @examples
#' m <- fm_tensor(list(
#'   space = fmexample$mesh,
#'   time = fm_mesh_1d(1:5)
#' ))
#' m2 <- fm_as_tensor(m)
#' m3 <- fm_as_tensor_list(list(m, m))
#' c(fm_dof(m$fun_spaces$space) * fm_dof(m$fun_spaces$time), fm_dof(m))
#' str(fm_evaluator(m, loc = list(space = cbind(0, 0), time = 2.5)))
#' str(fm_basis(m, loc = list(space = cbind(0, 0), time = 2.5)))
#' str(fm_fem(m))
fm_tensor <- function(x, ...) {
  nn <- names(x)
  if (is.null(nn)) {
    nn <- as.character(seq_along(x))
  } else if (any(nn == "")) {
    stop(
      "All or no elements of the list of",
      " function space objects need to be named."
    )
  }
  m <- structure(
    list(
      fun_spaces = fm_as_list(x),
      manifold = ""
    ),
    class = "fm_tensor"
  )
  type <- vapply(m$fun_spaces, fm_manifold_type, character(1))
  d <- vapply(m$fun_spaces, fm_manifold_dim, integer(1))
  type <- unique(type)
  if ((length(type) == 1L) && (type %in% c("R", "S", "T", "M"))) {
    m$manifold <- paste0(list(
      R = "R",
      S = "T",
      T = "T",
      M = "M"
    )[[type]], sum(d))
  } else {
    m$manifold <- paste0("M", sum(d))
  }
  m
}

#' @title Convert objects to `fm_tensor`
#' @describeIn fm_as_tensor Convert an object to `fm_tensor`.
#' @param x Object to be converted.
#' @param ... Arguments passed on to submethods
#' @returns An `fm_tensor` object
#' @export
#' @family object creation and conversion
#' @export
#' @examples
#' fm_as_tensor_list(list(fm_tensor(list())))
#'
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
