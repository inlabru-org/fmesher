#' @include deprecated.R

# fm_list ####

fm_class_stubs <- function() {
  c("segm", "mesh_1d", "mesh_2d", "lattice_2d", "tensor")
}

fm_class_stub <- function(x) {
  if (!is.character(x)) {
    x <- class(x)[1]
  }
  if (!grepl("^fm_", x)) {
    return(NULL)
  }
  x <- sub("^fm_", "", x)
  if (grepl("^list$", x)) {
    return(NULL)
  }
  x <- sub("_list$", "", x)
  if (!(x %in% fm_class_stubs())) {
    return(NULL)
  }
  x
}

#' @title Handle lists of fmesher objects
#'
#' @description
#' Methods for constructing and manipulating `fm_list` objects.
#'
#' @describeIn fm_list Convert each element of a list, or convert a single
#'   non-list object and return in a list
#' @param x list of objects to be converted.
#' @param ... Arguments passed to each individual conversion call.
#' @param .class_stub character; class stub name of class to convert each list
#'   element to. If `NULL`, uses `fm_as_fm` and auto-detects if the resulting
#'   list has consistent class, and then adds that to the class list.  If
#'   non-null, uses `paste0("fm_as_", .class_stub)` for conversion, and verifies
#'   that the resulting list has elements consistent with that class.
#' @returns An `fm_list` object, potentially with `fm_{class_stub}_list`
#'   added.
#' @export
fm_list <- function(x, ..., .class_stub = NULL) {
  fm_as_list(x, ..., .class_stub = NULL)
}
#' @describeIn fm_list Convert each element of a list, or convert a single
#'   non-list object and return in a list
#' @export
fm_as_list <- function(x, ..., .class_stub = NULL) {
  if (is.null(.class_stub)) {
    .method <- "fm_as_fm"
    .class_name <- NULL
    .class_list_name <- NULL
  } else {
    .method <- paste0("fm_as_", .class_stub)
    .class_name <- paste0("fm_", .class_stub)
    .class_list_name <- paste0(.class_name, "_list")
    if (!missing(x) && inherits(x, .class_list_name)) {
      return(x)
    }
  }
  if (missing(x) || is.null(x) || (length(x) == 0)) {
    return(structure(list(),
                     class = c(.class_list_name, "fm_list")
    ))
  }

  if (inherits(x, paste0("fm_", fm_class_stubs(), "_list"))) {
    return(x)
  }

  if (!inherits(x, "fm_list")) {
    m_c <- setdiff(method_classes(.method), "list")
    if (inherits(x, paste0("fm_", fm_class_stubs())) ||
        (!is.null(m_c) && inherits(x, m_c))) {
      # Single element of known or coercible non-list type
      #      y <- do.call(.method, list(x, ...))
      return(fm_as_list(list(x), ..., .class_stub = .class_stub))
    }
  }

  if (!inherits(x, "list")) {
    stop(paste0(
      "'list' object expected. Received '",
      paste0(class(x), collapse = ", "),
      "'."
    ))
  }

  y <- lapply(x, function(xx) do.call(.method, list(xx, ...)))

  if ((length(y) > 0) && is.null(.class_stub)) {
    stubs <- fm_class_stubs()
    is_stub <- vapply(
      stubs, function(stub) {
        all(vapply(y, function(yy) {
          is.null(yy) || inherits(yy, paste0("fm_", stub))
        }, TRUE))
      },
      TRUE
    )
    if (any(is_stub)) {
      .class_stub <- stubs[is_stub][1]
    }
  }
  if (!is.null(.class_stub)) {
    .class_name <- paste0("fm_", .class_stub)
    if (length(y) > 0) {
      is_stub <-
        all(vapply(
          y,
          function(yy) {
            is.null(yy) || inherits(yy, .class_name)
          },
          TRUE
        ))
      if (!is_stub) {
        stop(
          "Inconsistent element classes for 'fm_list' for class '",
          .class_name, "'"
        )
      }
    }
    .class_list_name <- paste0(.class_name, "_list")
    class(y) <- c(.class_list_name, "fm_list")
    return(y)
  }

  class(y) <- "fm_list"
  return(y)
}



#' @export
#' @describeIn fm_list The `...` arguments should be coercible to `fm_list`
#' objects.
`c.fm_list` <- function(...) {
  if (!all(vapply(
    list(...),
    function(xx) is.null(xx) || inherits(xx, "fm_list"),
    TRUE
  ))) {
    y <- lapply(list(...), fm_as_list)
    return(do.call("c", y))
  }
  object <- NextMethod()
  fm_as_list(object)
}

#' @export
#' @param x `fm_list` object from which to extract element(s)
#' @param i indices specifying elements to extract
#' @describeIn fm_list Extract sub-list
`[.fm_list` <- function(x, i) {
  object <- NextMethod()
  class(object) <- class(x)
  object
}


