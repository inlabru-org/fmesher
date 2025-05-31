#' @include deprecated.R

# fm_collect ####

#' @title Make a collection function space
#' @description `r lifecycle::badge("experimental")`
#' Collection function spaces. The interface and object storage model
#' is experimental and may change.
#' @export
#' @param x list of function space objects, such as [fm_mesh_2d()], all of the
#' same type.
#' @param ... Currently unused
#' @returns A `fm_collect` or `fm_collect_list` object.
#'   Elements of `fm_collect`:
#' \describe{
#' \item{fun_spaces}{`fm_list` of function space objects}
#' \item{manifold}{character; manifold type summary, obtained from the
#'   function spaces.}
#' }
#' @family object creation and conversion
#' @examples
#' m <- fm_collect(list(
#'   fmexample$mesh,
#'   fmexample$mesh
#' ))
#' m2 <- fm_as_collection(m)
#' m3 <- fm_as_collection_list(list(m, m))
#' c(fm_dof(m$fun_spaces[[1]]) + fm_dof(m$fun_spaces[[2]]), fm_dof(m))
#' fm_basis(m, loc = tibble::tibble(loc = cbind(0, 0), index = 2), full = TRUE)
#' fm_evaluator(m, loc = tibble::tibble(loc = cbind(0, 0), index = 2))
#' names(fm_fem(m))
#' fm_diameter(m)
fm_collect <- function(x, ...) {
  m <- structure(
    list(
      fun_spaces = fm_as_list(x),
      manifold = ""
    ),
    class = "fm_collect"
  )
  type <- vapply(m$fun_spaces, fm_manifold, character(1))
  type <- unique(type)
  if (length(type) > 1L) {
    stop(
      "All function spaces in a collection need to be of ",
      "the same manifold type, ",
      "but found: ", paste(type, collapse = ", ")
    )
  }
  m$manifold <- type
  m
}

#' @title Convert objects to `fm_collect`
#' @describeIn fm_as_collection Convert an object to `fm_collect`.
#' @param x Object to be converted.
#' @param ... Arguments passed on to submethods
#' @returns An `fm_collect` object
#' @export
#' @family object creation and conversion
#' @export
#' @examples
#' fm_as_collection_list(list(fm_collect(list())))
#'
fm_as_collection <- function(x, ...) {
  if (is.null(x)) {
    return(NULL)
  }
  UseMethod("fm_as_collection")
}
#' @describeIn fm_as_collection Convert each element of a list
#' @export
fm_as_collection_list <- function(x, ...) {
  fm_as_list(x, ..., .class_stub = "collection")
}
#' @rdname fm_as_collection
#' @param x Object to be converted
#' @export
fm_as_collection.fm_collect <- function(x, ...) {
  #  class(x) <- c("fm_collect", setdiff(class(x), "fm_collect"))
  x
}
