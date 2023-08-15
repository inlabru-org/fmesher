#' @include deprecated.R

# fm_manifold ####

#' @title Query the mesh manifold type
#' @description
#' Extract a manifold definition string, or a logical for matching
#' manifold type
#' @param x A [fm_mesh_1d] or [fm_mesh_2d] object (or other object containing a
#' `manifold` element)
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
  if (is.null(type)) {
    return(x[["manifold"]])
  }
  if (is.null(x) || is.null(x[["manifold"]])) {
    return(FALSE)
  }
  # Match exact manifold?
  if (x[["manifold"]] %in% type) {
    return(TRUE)
  }
  # Match space name or dimension?
  m <- intersect(c("M", "R", "S"), type)
  d <- intersect(as.character(seq_len(3)), type)
  if (length(c(m, d)) == 0) {
    return(FALSE)
  }
  grepl(paste0(c(m, d), collapse = "|"), x[["manifold"]])
}

#' @rdname fm_manifold
#' @returns `fm_manifold_type()`: character or NULL; "M", "R", or "S"
#' @export
fm_manifold_type <- function(x) {
  if (is.null(x) || is.null(x[["manifold"]])) {
    return(NULL)
  }
  # Match space name or dimension?
  m <- intersect(c("M", "R", "S"), strsplit(x[["manifold"]], "")[[1]])
  if (length(m) == 0) {
    return(NULL)
  }
  m
}

#' @rdname fm_manifold
#' @returns `fm_manifold_dim()`: integer or NULL
#' @export
fm_manifold_dim <- function(x) {
  if (is.null(x) || is.null(x[["manifold"]])) {
    return(NULL)
  }
  # Match space name or dimension?
  d <- intersect(as.character(seq_len(3)), strsplit(x[["manifold"]], "")[[1]])
  if (length(d) == 0) {
    return(NULL)
  }
  as.integer(d)
}
