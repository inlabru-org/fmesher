#' @include deprecated.R

# fm_manifold ####

#' @title Query the mesh manifold type
#' @description
#' Extract a manifold definition string, or a logical for matching
#' manifold type
#' @param x An object with `manifold` information, or a character string
#' @param type `character`; if `NULL` (the default), returns the manifold
#' definition string by calling `fm_manifold_get(x)`.
#' If `character`, returns `TRUE` if the manifold type of `x` matches at least
#' one of the character vector elements.
#' @returns `fm_manifold()`: Either logical (matching manifold type yes/no),
#' or character (the stored manifold, when `is.null(type)` is `TRUE`)
#' @export
#' @examples
#' fm_manifold_get(fmexample$mesh)
#' fm_manifold(fmexample$mesh)
#' fm_manifold(fmexample$mesh, "R2")
#' fm_manifold_type(fmexample$mesh)
#' fm_manifold_dim(fmexample$mesh)
fm_manifold <- function(x, type = NULL) {
  x <- fm_manifold_get(x)
  if (is.null(type)) {
    return(x)
  }
  if (is.null(x)) {
    return(FALSE)
  }
  any(vapply(type, function(t) fm_manifold_match(x, t), logical(1)))
}

#' @describeIn fm_manifold Method for obtaining a text representation of the
#' manifold characteristics, e.g. "R1", "R2", "M2", or "T3". The default
#' method assumes that the manifold is stored as a `character` string in a
#' "manifold" element of the object, so it can be extracted with
#' `x[["manifold"]]`.
#' Object classes that do not store the information in this way need to
#' implement their own method.
#' @returns `fm_manifold_get()`: `character` or `NULL`
#' @export
fm_manifold_get <- function(x) {
  UseMethod("fm_manifold_get")
}

#' @rdname fm_manifold
#' @export
fm_manifold_get.default <- function(x) {
  x[["manifold"]]
}

#' @rdname fm_manifold
#' @export
fm_manifold_get.character <- function(x) {
  x
}

# Check match for a single type
fm_manifold_match <- function(x, type) {
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
#' @returns `fm_manifold_type()`: character or NULL; "M" (curved manifold),
#' "R" (flat space), "S" (generalised spherical space), "T"
#' (general tensor product space), or "G" (metric graph)
#' @export
fm_manifold_type <- function(x) {
  x <- fm_manifold_get(x)
  if (is.null(x)) {
    return(NULL)
  }

  splt <- strsplit(x, "")[[1]]
  splt <- splt[splt %in% c("M", "R", "S", "T", "G")]
  if (length(splt) == 0) {
    return(NULL)
  }
  paste0(splt, collapse = "")
}

#' @rdname fm_manifold
#' @returns `fm_manifold_dim()`: integer or NULL
#' @export
fm_manifold_dim <- function(x) {
  x <- fm_manifold_get(x)
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
