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
#' Extracts converts an [fm_mesh_2d()] object to `sf` with [fm_as_sfc()] and uses
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
#'     max.edge = c(0.1, 0.5),
#'     crs = fm_crs("+proj=longlat")
#'   )
#'   ggplot() +
#'     geom_fm(data = m, crs = fm_crs("epsg:27700"))
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
#' Extracts converts an [fm_segm()] object to `sf` with [fm_as_sfc()] and uses
#' `geom_sf` to visualize it.
#' @export
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
