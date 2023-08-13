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
#' @examples
#' if (require("ggplot2", quietly = TRUE)) {
#'   ggplot() +
#'     geom_fm(data = fmexample$mesh)
#' }
#'
geom_fm <- function(mapping = NULL, data = NULL, ...) {
  UseMethod("geom_fm", data)
}



#' @describeIn geom_fm
#' Converts an [fm_mesh_2d()] object to `sf` with [fm_as_sfc()] and uses
#' `geom_sf` to visualize the triangles and edges.
#'
#' @export
#' @param mapping_int `aes` for interior constraint edges.
#' @param mapping_bnd `aes` for boundary edges.
#' @param defs_int additional settings for interior constraint edges.
#' @param defs_bnd additional settings for boundary edges.
#' @param crs Optional crs to transform the object to before plotting.
#' @importFrom utils modifyList
#' @examples
#' if (require("ggplot2", quietly = TRUE)) {
#'   m <- fm_mesh_2d(
#'     cbind(10, 20),
#'     boundary = fm_extensions(cbind(10, 20), c(25, 65)),
#'     max.edge = c(2, 5),
#'     crs = fm_crs("+proj=longlat")
#'   )
#'   print(ggplot() +
#'     geom_fm(data = m))
#'   print(ggplot() +
#'     geom_fm(data = m, crs = fm_crs("epsg:27700")))
#' }
#' \donttest{
#' if (require("ggplot2", quietly = TRUE)) {
#'   # Compute a mesh vertex based function on a different grid
#'   px <- fm_pixels(fm_transform(m, fm_crs("mollweide_globe")))
#'   px$fun <- fm_evaluate(m,
#'     loc = px,
#'     field = sin(m$loc[, 1] / 5) * sin(m$loc[, 2] / 5)
#'   )
#'   print(ggplot() +
#'     geom_tile(aes(geometry = geometry, fill = fun),
#'       data = px,
#'       stat = "sf_coordinates"
#'     ) +
#'     geom_fm(
#'       data = m, alpha = 0.2, linewidth = 0.05,
#'       crs = fm_crs("mollweide_globe")
#'     ))
#' }
#' }
geom_fm.fm_mesh_2d <- function(mapping = NULL,
                               data = NULL,
                               mapping_int = NULL,
                               mapping_bnd = NULL,
                               defs_int = NULL,
                               defs_bnd = NULL,
                               ...,
                               crs = NULL) {
  if (!is.null(crs)) {
    data <- fm_transform(data, crs = crs)
  }

  if (!fm_manifold(data, "R2")) {
    stop("geom_fm.fm_mesh_2d only implemented for R2 meshes.")
  }

  mesh_sf <- fm_as_sfc(data)
  int <- fm_segm(data, boundary = FALSE)
  bnd <- fm_segm(data, boundary = TRUE)

  maps <-
    list(
      mesh = mapping,
      int = mapping_int,
      bnd = mapping_bnd
    )
  maps_def <- list(
    mesh = ggplot2::aes(),
    int = ggplot2::aes(),
    bnd = ggplot2::aes()
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
      int = defs_int,
      bnd = defs_bnd
    )
  defs_def <- list(
    mesh = list(linewidth = 0.25, color = "grey"),
    int = list(linewidth = 0.5, color = "blue"),
    bnd = list(linewidth = 1, color = "black")
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

  c(
    do.call(ggplot2::geom_sf, c(
      list(mapping = maps$mesh, data = mesh_sf), defs$mesh
    )),
    do.call(geom_fm, c(
      list(mapping = maps$int, data = int), defs$int
    )),
    do.call(geom_fm, c(
      list(mapping = maps$bnd, data = bnd), defs$bnd
    ))
  )
}



#' @describeIn geom_fm
#' Converts an [fm_segm()] object to `sf` with [fm_as_sfc()] and uses
#' `geom_sf` to visualize it.
#' @export
#' @examples
#' if (require("ggplot2", quietly = TRUE)) {
#'   m <- fm_mesh_1d(c(1, 2, 4, 6, 10), boundary = c("n", "d"), degree = 2)
#'   ggplot() +
#'     geom_fm(data = m, weights = c(4, 2, 4, -1))
#' }
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
#' @param knots logical; if true, show the spline knot locations
#' @param derivatives logical; if true, draw first order derivatives instead of
#' function values
#' @param weights numeric vector; if provided, draw weighted basis functions and
#' the resulting weighted sum.
#' @export
#' @importFrom utils modifyList
#' @examples
#' if (require("ggplot2", quietly = TRUE)) {
#'   m <- fm_mesh_1d(
#'     c(1, 2, 3, 5, 7),
#'     boundary = c("dirichlet", "neumann"),
#'     degree = 2
#'   )
#'   print(ggplot() +
#'     geom_fm(data = m))
#' }
geom_fm.fm_mesh_1d <- function(mapping = NULL,
                               data = NULL,
                               ...,
                               xlim = NULL,
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
    c(data$loc, data$interval[2])
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

  result <-
    c(
      do.call(ggplot2::geom_line, c(
        list(mapping = maps$basis, data = df), defs$basis
      ))
    )
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
