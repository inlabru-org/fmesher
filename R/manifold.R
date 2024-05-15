#' @include deprecated.R

# fm_manifold ####

#' @title Query the mesh manifold type
#' @description
#' Extract a manifold definition string, or a logical for matching
#' manifold type
#' @param x An object with `manifold` information, or a character string
#' @param type `character`; if `NULL` (the default), returns the manifold definition string.
#' If `character`, returns `TRUE` if the manifold type of `x` matches at least
#' one of the character vector elements.
#' @returns `fm_manifold()`: Either logical (matching manifold type yes/no),
#' or character (the stored manifold, when `is.null(type)` is `TRUE`)
#' @export
#' @examples
#' fm_manifold(fmexample$mesh)
#' fm_manifold_type(fmexample$mesh)
#' fm_manifold_dim(fmexample$mesh)
fm_manifold <- function(x, type = NULL) {
  if (!is.character(x)) {
    x <- x[["manifold"]]
  }
  if (is.null(type)) {
    return(x)
  }
  if (is.null(x)) {
    return(FALSE)
  }
  # Match exact manifold?
  if (x %in% type) {
    return(TRUE)
  }
  m <- fm_manifold_type(type)
  if (!is.null(m) && !identical(m, fm_manifold_type(x))) {
    return(FALSE)
  }
  d <- fm_manifold_dim(type)
  if (!is.null(d) && !identical(d, fm_manifold_dim(x))) {
    return(FALSE)
  }
  if (is.null(m) && is.null(d)) {
    return(TRUE)
  }
  TRUE
}

#' @rdname fm_manifold
#' @returns `fm_manifold_type()`: character or NULL; "M", "R", "S", or "T"
#' @export
fm_manifold_type <- function(x) {
  if (!is.character(x)) {
    x <- x[["manifold"]]
  }
  if (is.null(x)) {
    return(NULL)
  }

  splt <- strsplit(x, "")[[1]]
  splt <- splt[splt %in% c("M", "R", "S", "T")]
  if (length(splt) == 0) {
    return(NULL)
  }
  paste0(splt, collapse = "")
}

#' @rdname fm_manifold
#' @returns `fm_manifold_dim()`: integer or NULL
#' @export
fm_manifold_dim <- function(x) {
  if (!is.character(x)) {
    x <- x[["manifold"]]
  }
  if (is.null(x)) {
    return(NULL)
  }

  figures <- as.character(c(0, seq_len(9)))
  splt <- strsplit(x, "")[[1]]
  splt <- splt[splt %in% figures]
  if (length(splt) == 0) {
    return(NULL)
  }
  d <- paste0(splt, collapse = "")
  as.integer(d)
}
