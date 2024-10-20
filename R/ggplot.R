#' @title ggplot2 geomes for fmesher related objects
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' `geom_fm` is a generic function for generating geomes from various kinds of
#' `fmesher` objects, e.g. `fm_segm` and `fm_mesh_2d`.
#' The function invokes particular methods which depend
#' on the [class] of the `data` argument.
#' Requires the `ggplot2` package.
#'
#' Note: `geom_fm` is not yet a "proper" `ggplot2` geom method; the interface
#' may therefore change in the future.
#'
#' @export
#' @param mapping an object for which to generate a geom.
#' @param data an object for which to generate a geom.
#' @param ... Arguments passed on to the geom method.
#' @returns A combination of `ggplot2` geoms.
#'
#' @examplesIf require("ggplot2", quietly = TRUE)
#' ggplot() +
#'   geom_fm(data = fmexample$mesh)
#'
geom_fm <- function(mapping = NULL, data = NULL, ...) {
  UseMethod("geom_fm", data)
}



#' @describeIn geom_fm
#' Converts an [fm_mesh_2d()] object to `sf` with [fm_as_sfc()] and uses
#' `geom_sf` to visualize the triangles and edges.
#'
#' The mesh vertices are only plotted if `mappings$loc` or `defs$loc`
#' is non-`NULL`, e.g. `defs = list(loc = list())`. Default argument settings:
#' ```
#' ... = linewidth = 0.25, color = "grey" # default for triangle mapping
#' defs = list(
#'   int = list(linewidth = 0.5, color = "blue"),
#'   bnd = list(linewidth = 1, color = "black"),
#'   loc = list(size = 1, color = "red")
#' )
#' ```
#'
#' @export
#' @param mappings optional list of `aes` mappings for the non-triangle parts of
#'   the mesh, named "int" for interior constraint edges, "bnd" for boundary
#'   edges, and "loc" for the vertices.
#' @param defs optional list of non-`aes` settings for the non-triangle parts of
#'   the mesh, named "int" for interior constraint edges, "bnd" for boundary
#'   edges, and "loc" for the vertices.
#' @param mapping_int,mapping_bnd,defs_int,defs_bnd
#' `r lifecycle::badge("deprecated")` arguments; see `mappings` and `defs`.
#' @param crs Optional crs to transform the object to before plotting.
#' @importFrom utils modifyList
#' @examplesIf require("ggplot2", quietly = TRUE)
#' m <- fm_mesh_2d(
#'   cbind(10, 20),
#'   boundary = fm_extensions(cbind(10, 20), c(25, 65)),
#'   max.edge = c(4, 10),
#'   crs = fm_crs("+proj=longlat")
#' )
#' ggplot() +
#'   geom_fm(data = m)
#' ggplot() +
#'   geom_fm(data = m, defs = list(loc = list()))
#' ggplot() +
#'   geom_fm(data = m, crs = fm_crs("epsg:27700"))
#' \donttest{
#' # Compute a mesh vertex based function on a different grid
#' px <- fm_pixels(fm_transform(m, fm_crs("mollweide_globe")))
#' px$fun <- fm_evaluate(m,
#'   loc = px,
#'   field = sin(m$loc[, 1] / 5) * sin(m$loc[, 2] / 5)
#' )
#' ggplot() +
#'   geom_tile(aes(geometry = geometry, fill = fun),
#'     data = px,
#'     stat = "sf_coordinates"
#'   ) +
#'   geom_fm(
#'     data = m, alpha = 0.2, linewidth = 0.05,
#'     crs = fm_crs("mollweide_globe")
#'   )
#' }
geom_fm.fm_mesh_2d <- function(mapping = NULL,
                               data = NULL,
                               ...,
                               mappings = NULL,
                               defs = NULL,
                               crs = NULL,
                               mapping_int = deprecated(),
                               mapping_bnd = deprecated(),
                               defs_int = deprecated(),
                               defs_bnd = deprecated()) {
  if (is.null(mappings)) {
    mappings <- list()
  }
  if (is.null(defs)) {
    mappings <- list()
  }
  if (lifecycle::is_present(mapping_int)) {
    lifecycle::deprecate_warn(
      "0.1.7.9009",
      "geom_fm(mapping_int)",
      "geom_fm(mappings = list(int = mapping_int))"
    )
    if (!is.null(mappings$int)) {
      stop("Both mapping_int and mappings$int are provided.")
    }
    mappings$int <- mapping_int
  }
  if (lifecycle::is_present(mapping_bnd)) {
    lifecycle::deprecate_warn(
      "0.1.7.9009",
      "geom_fm(mapping_bnd)",
      "geom_fm(mappings = list(bnd = mapping_bnd))"
    )
    if (!is.null(mappings$bnd)) {
      stop("Both mapping_bnd and mappings$bnd are provided.")
    }
    mappings$bnd <- mapping_bnd
  }
  if (lifecycle::is_present(defs_int)) {
    lifecycle::deprecate_warn(
      "0.1.7.9009",
      "geom_fm(defs_int)",
      "geom_fm(defs = list(int = defs_int))"
    )
    if (!is.null(defs$int)) {
      stop("Both defs_int and mappings$int are provided.")
    }
    mappings$int <- defs_int
  }
  if (lifecycle::is_present(defs_bnd)) {
    lifecycle::deprecate_warn(
      "0.1.7.9009",
      "geom_fm(defs_bnd)",
      "geom_fm(defs = list(bnd = defs_bnd))"
    )
    if (!is.null(defs$bnd)) {
      stop("Both defs_bnd and mappings$bnd are provided.")
    }
    mappings$bnd <- defs_bnd
  }

  if (!is.null(crs)) {
    data <- fm_transform(data, crs = crs)
  }

  if (!fm_manifold(data, "R2")) {
    stop("geom_fm.fm_mesh_2d only implemented for R2 meshes.")
  }

  mesh_sf <- fm_as_sfc(data)
  int <- fm_segm(data, boundary = FALSE)
  bnd <- fm_segm(data, boundary = TRUE)
  loc <- fm_as_sfc(data, format = "loc")

  maps <-
    list(
      mesh = mapping,
      int = mappings$int,
      bnd = mappings$bnd,
      loc = mappings$loc
    )
  maps_def <- list(
    mesh = ggplot2::aes(),
    int = ggplot2::aes(),
    bnd = ggplot2::aes(),
    loc = ggplot2::aes()
  )
  maps <- lapply(
    names(maps_def),
    function(x) {
      if (is.null(maps[[x]])) {
        maps_def[[x]]
      } else {
        modifyList(maps_def[[x]], maps[[x]])
      }
    }
  )
  names(maps) <- names(maps_def)

  defs <-
    list(
      mesh = list(...),
      int = defs$int,
      bnd = defs$bnd,
      loc = defs$loc
    )
  defs_def <- list(
    mesh = list(linewidth = 0.25, color = "grey"),
    int = list(linewidth = 0.5, color = "blue"),
    bnd = list(linewidth = 1, color = "black"),
    loc = list(size = 1, color = "red")
  )
  def <- lapply(
    names(defs_def),
    function(x) {
      if (is.null(defs[[x]])) {
        defs_def[[x]]
      } else {
        modifyList(defs_def[[x]], defs[[x]])
      }
    }
  )
  names(def) <- names(defs_def)

  geoms <-
    c(
      do.call(ggplot2::geom_sf, c(
        list(mapping = maps$mesh, data = mesh_sf), def$mesh
      )),
      do.call(geom_fm, c(
        list(mapping = maps$int, data = int), def$int
      )),
      do.call(geom_fm, c(
        list(mapping = maps$bnd, data = bnd), def$bnd
      ))
    )
  if (!is.null(defs$loc) || !is.null(mappings$loc)) {
    geoms <- c(
      geoms,
      do.call(ggplot2::geom_sf, c(
        list(mapping = maps$loc, data = loc), def$loc
      ))
    )
  }
  geoms
}



#' @describeIn geom_fm
#' Converts an [fm_segm()] object to `sf` with [fm_as_sfc()] and uses
#' `geom_sf` to visualize it.
#' @export
#' @examplesIf require("ggplot2", quietly = TRUE)
#' m <- fm_mesh_1d(c(1, 2, 4, 6, 10), boundary = c("n", "d"), degree = 2)
#' ggplot() +
#'   geom_fm(data = m, weights = c(4, 2, 4, -1))
#'
geom_fm.fm_segm <- function(mapping = NULL,
                            data = NULL,
                            ...,
                            crs = NULL) {
  if (!is.null(crs)) {
    data <- fm_transform(data, crs = crs)
  }

  fm_is_bnd(data) <- FALSE # Avoid warning from fm_as_sfc
  segm_sf <- fm_as_sfc(data)

  maps <-
    list(
      segm = mapping
    )
  maps_def <- list(
    segm = ggplot2::aes()
  )
  maps <- lapply(
    names(maps_def),
    function(x) {
      if (is.null(maps[[x]])) {
        maps_def[[x]]
      } else {
        modifyList(maps_def[[x]], maps[[x]])
      }
    }
  )
  names(maps) <- names(maps_def)

  defs <-
    list(
      segm = list(...)
    )
  defs_def <- list(
    segm = list()
  )
  defs <- lapply(
    names(defs_def),
    function(x) {
      if (is.null(defs[[x]])) {
        defs_def[[x]]
      } else {
        modifyList(defs_def[[x]], defs[[x]])
      }
    }
  )
  names(defs) <- names(defs_def)

  do.call(ggplot2::geom_sf, c(
    list(mapping = maps$segm, data = segm_sf), defs$segm
  ))
}





#' @describeIn geom_fm
#' Evaluates and plots the basis functions defined by an [fm_mesh_1d()] object.
#'
#' @param xlim numeric 2-vector; specifies the interval for which to compute
#' functions. Default is `data$interval`
#' @param basis logical; if `TRUE` (default), show the spline basis functions
#' @param knots logical; if `TRUE` (default), show the spline knot locations
#' @param derivatives logical; if `TRUE` (not default), draw first order
#' derivatives instead of function values
#' @param weights numeric vector; if provided, draw weighted basis functions and
#' the resulting weighted sum.
#' @export
#' @importFrom utils modifyList
#' @examplesIf require("ggplot2", quietly = TRUE)
#' m <- fm_mesh_1d(
#'   c(1, 2, 3, 5, 7),
#'   boundary = c("dirichlet", "neumann"),
#'   degree = 2
#' )
#' ggplot() +
#'   geom_fm(data = m)
#'
geom_fm.fm_mesh_1d <- function(mapping = NULL,
                               data = NULL,
                               ...,
                               xlim = NULL,
                               basis = TRUE,
                               knots = TRUE,
                               derivatives = FALSE,
                               weights = NULL) {
  if (is.null(xlim)) {
    xlim <- data$interval
  }
  x <- seq(xlim[1], xlim[2], length.out = data$n * 100)
  A <- as.matrix(fm_basis(data, loc = x))
  A_d1 <- rbind(
    A[2, , drop = FALSE] - A[1, , drop = FALSE],
    (A[seq_len(nrow(A) - 2) + 2, , drop = FALSE] -
      A[seq_len(nrow(A) - 2), , drop = FALSE]),
    A[nrow(A), ] - A[nrow(A) - 1, , drop = FALSE]
  ) / c(
    x[2] - x[1], x[seq_len(nrow(A) - 2) + 2] - x[seq_len(nrow(A) - 2)],
    x[length(x)] - x[length(x) - 1]
  )

  df <- data.frame(
    x = rep(x, times = data$m),
    basis = factor(rep(seq_len(data$m), each = length(x))),
    value = as.vector(A),
    derivative = as.vector(A_d1)
  )
  if (!is.null(weights)) {
    df$value <- as.vector(A %*% Matrix::Diagonal(ncol(A), weights))
    df$derivative <- as.vector(A_d1 %*% Matrix::Diagonal(ncol(A), weights))
    if (derivatives) {
      df_fun <- data.frame(
        x = x,
        value = as.vector(A_d1 %*% weights)
      )
    } else {
      df_fun <- data.frame(
        x = x,
        value = as.vector(A %*% weights)
      )
    }
  }

  knots_ <- if (data$cyclic) {
    knots_ <- sort(unique((data$loc - data$interval[1]) %% diff(data$interval)))
    c(
      knots_[length(knots_)] - diff(data$interval),
      knots_,
      knots_[1] + diff(data$interval)
    ) + data$interval[1]
  } else {
    data$loc
  }
  df_knots <- data.frame(
    knots = knots_[(knots_ >= xlim[1]) & (knots_ <= xlim[2])]
  )

  maps <-
    list(
      basis = mapping,
      knots = ggplot2::aes(),
      fun = ggplot2::aes()
    )
  if (derivatives) {
    maps_def <- list(basis = ggplot2::aes(
      x = .data[["x"]],
      y = .data[["derivative"]],
      color = .data[["basis"]]
    ))
  } else {
    maps_def <- list(basis = ggplot2::aes(
      x = .data[["x"]],
      y = .data[["value"]],
      color = .data[["basis"]]
    ))
  }
  maps_def$knots <- ggplot2::aes(xintercept = .data[["knots"]])
  maps_def$fun <- ggplot2::aes(x = .data[["x"]], y = .data[["value"]])
  maps <- lapply(
    names(maps_def),
    function(x) {
      if (is.null(maps[[x]])) {
        maps_def[[x]]
      } else {
        modifyList(maps_def[[x]], maps[[x]])
      }
    }
  )
  names(maps) <- names(maps_def)

  defs <-
    list(
      basis = list(...),
      knots = list()
    )
  defs_def <- list(
    basis = list(),
    knots = list(linewidth = 0.5, alpha = 0.5)
    #    basis = list(linewidth = 0.25),
    #    knots = list(linewidth = 0.25)
  )
  defs <- lapply(
    names(defs_def),
    function(x) {
      if (is.null(defs[[x]])) {
        defs_def[[x]]
      } else {
        modifyList(defs_def[[x]], defs[[x]])
      }
    }
  )
  names(defs) <- names(defs_def)

  result <- list()
  if (basis) {
    result <- c(
      result,
      do.call(ggplot2::geom_line, c(
        list(mapping = maps$basis, data = df), defs$basis
      ))
    )
  }
  if (knots) {
    result <- c(
      result,
      do.call(ggplot2::geom_vline, c(
        list(mapping = maps$knots, data = df_knots), defs$knots
      ))
    )
  }
  if (!is.null(weights)) {
    result <- c(
      result,
      do.call(ggplot2::geom_line, c(
        list(mapping = maps$fun, data = df_fun), defs$fun
      ))
    )
  }

  result
}
