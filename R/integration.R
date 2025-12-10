#' @include deprecated.R


#' @title (Blockwise) cross product of integration points
#'
#' @description
#' Calculates the groupwise cross product of integration points in different
#' dimensions and multiplies their weights accordingly.
#' If the object defining points in a particular dimension has no
#' weights attached to it all weights are assumed to be 1.
#'
#'
#' @export
#' @keywords internal
#'
#' @param ... `tibble`, `data.frame`, `sf`, or `SpatialPointsDataFrame` objects,
#'   each one usually obtained by a call to an [fm_int()] method.
#' @param na.rm logical; if `TRUE`, the rows with weight `NA` from the
#'   non-overlapping full_join will be removed; if `FALSE`, set the undefined
#'   weights to `NA`. If `NULL` (default), act as `TRUE`, but warn if any
#'   elements needed removing.
#' @param .blockwise logical; if `FALSE`, computes full tensor product
#'   integration. If `TRUE`, computes within-block tensor product integration
#'   (used internally by [fm_int()]). Default `FALSE`
#' @returns A `data.frame`, `sf`, or `SpatialPointsDataFrame` of
#'   multidimensional integration points and their weights
#'
#' @examples
#' if (require("ggplot2")) {
#'   # Create integration points in dimension 'myDim' and 'myDiscreteDim'
#'   ips1 <- fm_int(fm_mesh_1d(1:20),
#'     rbind(c(0, 3), c(3, 8)),
#'     name = "myDim"
#'   )
#'   ips2 <- fm_int(domain = c(1, 2, 4), name = "myDiscreteDim")
#'
#'   # Calculate the cross product
#'   ips <- fm_cprod(ips1, ips2)
#'
#'   # Plot the integration points
#'   ggplot(ips) +
#'     geom_point(aes(myDim, myDiscreteDim, size = weight)) +
#'     scale_size_area()
#' }
#'
#' @importFrom stats na.omit
fm_cprod <- function(..., na.rm = NULL, .blockwise = FALSE) {
  ipl <- list(...)

  # Transform sp to sf
  # TODO make a test. and give a warning for NA non-overlapping outcome?
  # check for each element, or on the subset, change only for sp anonymous
  # function on lapply
  ipl_sp <- vapply(ipl, function(x) inherits(x, "Spatial"), TRUE)
  ipl_sf <- vapply(ipl, function(x) inherits(x, c("sf", "sfc")), TRUE)
  ipl[ipl_sp] <- lapply(ipl[ipl_sp], fm_Spatial_as_int_object)

  ipl <- ipl[!vapply(ipl, is.null, TRUE)]
  if (length(ipl) == 0) {
    return(NULL)
  }

  nms <- names(ipl)
  if (is.null(nms)) {
    nms <- paste0(seq_along(ipl))
    names(ipl) <- nms
  } else if (any(nms == "" | is.na(nms))) {
    nms[nms == "" | is.na(nms)] <-
      paste0(seq_along(ipl))[nms == "" | is.na(nms)]
    names(ipl) <- nms
  }

  if (length(ipl) == 1) {
    ips <- new_fm_int(ipl[[1]], name = names(ipl)[1])
    nms1 <- names(ips)
    nms2 <- c("weight", ".block", ".block_origin")
  } else {
    ips1 <- new_fm_int(ipl[[1]], name = names(ipl)[1])
    if (length(ipl) > 2) {
      ips2 <- do.call(
        fm_cprod,
        c(ipl[-1], list(na.rm = na.rm, .blockwise = .blockwise))
      )
    } else {
      ips2 <- new_fm_int(ipl[[2]], name = names(ipl)[2])
    }

    by <- setdiff(intersect(names(ips1), names(ips2)), "weight")
    if (.blockwise) {
      by <- setdiff(by, c(".block_origin"))
    } else {
      by <- setdiff(by, c(".block", ".block_origin"))
    }

    # `sf::st_join` performs spatial join/filter; `dplyr::*_join` expects `x` of
    # class `sf` and `y` of class `data.frame`. The trick `as.tibble(sf_obj)`
    # allows `dplyr::full_join` and turn it back to `sf` with active geometry as
    # the ips1.
    # Z <- full_join(as_tibble(X), as_tibble(Y), by = "group")
    # st_as_sf(Z)
    # https://stackoverflow.com/questions/64365792/
    #   dplyr-full-join-on-geometry-columns-of-sf-objects
    if (inherits(ips1, c("sf", "sfc")) ||
      inherits(ips2, c("sf", "sfc"))) {
      if (length(by) == 0) {
        ips <-
          sf::st_as_sf(
            dplyr::cross_join(
              tibble::as_tibble(ips2),
              tibble::as_tibble(ips1)
            ),
            sf_column_name = if (inherits(ips1, "sf")) {
              attr(ips1, "sf_column", exact = TRUE)
            } else if (inherits(ips2, "sf")) {
              attr(ips2, "sf_column", exact = TRUE)
            } else {
              NULL
            }
          )
      } else {
        ips <-
          sf::st_as_sf(
            dplyr::full_join(
              tibble::as_tibble(ips2),
              tibble::as_tibble(ips1),
              by = by,
              relationship = "many-to-many"
            ),
            sf_column_name = if (inherits(ips1, "sf")) {
              attr(ips1, "sf_column", exact = TRUE)
            } else if (inherits(ips2, "sf")) {
              attr(ips2, "sf_column", exact = TRUE)
            } else {
              NULL
            }
          )
      }
    } else {
      # equivalent to base::merge(ips2, ips1, by = by, all = TRUE)
      if (length(by) == 0) {
        ips <-
          dplyr::cross_join(ips2, ips1)
      } else {
        ips <-
          dplyr::full_join(ips2, ips1,
            by = by,
            relationship = "many-to-many"
          )
      }
    }

    ips$weight <- ips$weight.x * ips$weight.y
    ips[["weight.x"]] <- NULL
    ips[["weight.y"]] <- NULL
    tibble::remove_rownames(ips)

    if (!.blockwise) {
      # Order by ips2, then ips1, so that ips1-indices change fastest
      levs <- paste(format(ips$.block.x), format(ips$.block.y), sep = "|")
      ordered_levs <- sort(unique(levs))
      ips$.block <- as.integer(factor(levs, levels = ordered_levs))
      ips[[".block.x"]] <- NULL
      ips[[".block.y"]] <- NULL
    }
    ips$.block_origin <- cbind(ips$.block_origin.y, ips$.block_origin.x)
    ips[[".block_origin.x"]] <- NULL
    ips[[".block_origin.y"]] <- NULL

    nms1 <- names(ips1)
    nms2 <- names(ips2)
  }

  # Reorder the columns to have ips1-only first, ips2-only second,
  # and joint last
  nms <- names(ips)
  nms_joint <- setdiff(
    nms,
    union(
      setdiff(nms1, nms2),
      setdiff(nms2, nms1)
    )
  )
  nms1 <- intersect(nms, setdiff(nms1, nms_joint))
  nms2 <- intersect(nms, setdiff(nms2, nms_joint))
  ips <- ips[, c(nms1, nms2, nms_joint), drop = FALSE]

  if (any(is.na(ips$weight)) && !isFALSE(na.rm)) {
    if (is.null(na.rm)) {
      warning(
        paste0(
          "Block information mismatch resulting in NA weights,",
          " and 'na.rm' was not supplied.",
          " These rows will be removed."
        )
      )
    }
    ips <- na.omit(ips)
  }

  ips <- dplyr::arrange(ips, .data$.block)

  # TODO Transform back to sp only if they are required. ips is a tibble sf tbl
  # data.frame.
  # It does not make sense to revert certain indices back after merging. Hence,
  # I revert the entire object back to sp.
  if (any(ipl_sp)) {
    if (any(ipl_sf)) {
      lifecycle::deprecate_stop(
        when = "0.0.1",
        what = "fm_cprod('...'='should not mix `sp` and `sf` objects')"
      )
    }
    ips <- fm_int_object_as_Spatial(ips)
  }
  ips
}

#' @title Construct integration scheme objects
#' @description Constructor method for integration scheme objects, allowing
#'   default construction of `.block` information. Primarily meant for internal
#'   use, but can be used to manually create data of the same structure as
#'   [fm_int()] output.
#'
#' @param object An object representing integration points; either a
#' data.frame-like object, or a vector/list of coordinates or other location
#' reference objects.
#' @param blocks logical; if `TRUE`, set per-element `.block` indices.
#' If `FALSE` (default), set a common block, `1L`.
#' @param weight Optional weight variable; if `NULL`, all weights are set to 1.
#' @param name character; name of the integration domain.
#' @param override logical; If `name` is non-NULL and `override=TRUE` for sf
#'   object, the current `sf_column` is renamed to `name`.
#' @returns A tibble or sf/tibble object. May acquire additional class
#' attributes in the future.
#' @seealso [fm_int()]
#' @export
#' @examples
#' new_fm_int(1:4, blocks = TRUE, weight = c(1, 2, 1, 3), name = "z")
new_fm_int <- function(object, blocks = FALSE, weight = NULL,
                       name = NULL, override = FALSE) {
  if (!is.data.frame(object)) {
    if (is.null(name) || (nzchar(name) == 0)) {
      stop("A dimension name must be provided for the integration points.")
    }
    .block <- if (blocks) {
      seq_len(NROW(object))
    } else {
      1L
    }
    if (inherits(object, "sfg")) {
      stop("fm_int_object() sfg input should be converted to sfc first.")
    }
    object <- tibble::tibble(
      "{name}" := object,
      weight = if (is.null(weight)) 1 else weight,
      .block = .block,
      .block_origin = matrix(.block, NROW(object), 1,
        dimnames = list(NULL, name)
      )
    )
    if (inherits(object[[name]], "sfc")) {
      object <- sf::st_as_sf(object, sf_column_name = name)
    }
  } else {
    if (!tibble::is_tibble(object)) {
      if (inherits(object, "sf")) {
        # Convert data.frame/sf to tibble/sf
        geometry_name <- attr(object, "sf_column", exact = TRUE)
        object <- tibble::as_tibble(object)
        object <- sf::st_as_sf(object, sf_column_name = geometry_name)
      } else {
        object <- tibble::as_tibble(object)
      }
    }
    if (override) {
      if (!is.null(name) && inherits(object, "sf")) {
        if (name %in% names(object)) {
          if (inherits(object[[name]], "sfc")) {
            sf::st_geometry(object) <- name
          }
        } else if (!is.null(attr(object, "sf_column", exact = TRUE))) {
          sf::st_geometry(object) <- name
        }
      }
    }
    if (is.null(object[["weight"]])) {
      if (is.null(weight)) {
        object[["weight"]] <- 1
      } else {
        object[["weight"]] <- weight
      }
    }
    if (is.null(object[[".block"]])) {
      if (blocks) {
        object[[".block"]] <- seq_len(NROW(object))
      } else {
        object[[".block"]] <- 1L
      }
    }
    if (is.null(object[[".block_origin"]])) {
      object[[".block_origin"]] <-
        matrix(object[[".block"]], NROW(object), 1,
          dimnames = list(NULL, name)
        )
    } else if (length(colnames(object[[".block_origin"]])) == 0) {
      colnames(object[[".block_origin"]]) <- name
    }
  }

  object
}

#' @describeIn fmesher-deprecated Deprecated function since `0.5.0.9013`;
#'   use [new_fm_int()] instead.
fm_int_object <- function(...) {
  lifecycle::deprecate_warn(
    when = "0.5.0.9013",
    what = "fm_int_object()",
    with = "new_fm_int()",
    always = TRUE
  )
  new_fm_int(...)
}

fm_int_object_as_Spatial <- function(ips) {
  geom <- sf::st_geometry(ips)
  sf_column <- attr(ips, "sf_column")
  ips <- tibble::as_tibble(ips)
  ips[[sf_column]] <- NULL
  sp::addAttrToGeom(sf::as_Spatial(geom, cast = TRUE, IDs = row.names(ips)),
    ips,
    match.ID = FALSE
  )
}
fm_Spatial_as_int_object <- function(ips) {
  sf::st_as_sf(ips)
  ##  crs <- fm_crs(ips)
  ##  coord_names <- sp::coordnames(ips)
  ##  ips <- tibble::as_tibble(ips)
  ##  sf::st_as_sf(ips, coords = coord_names, crs = crs,
  ##               sf_column_name = "geometry")
}

#' @title Multi-domain integration
#'
#' @description Construct integration points on tensor product spaces
#' @param domain Functional space specification; single domain or a named list
#' of domains
#' @param samplers For single domain `fm_int` methods, an object specifying one
#'   or more subsets of the domain, and optional weighting in a `weight`
#'   variable. For `fm_int.list`, a list of sampling definitions, where data
#'   frame elements may contain information for multiple domains, in which case
#'   each row represent a separate tensor product integration subspace.
#' @param name For single-domain methods, the variable name to use for the
#' integration points. Default 'x'
#' @param \dots Additional arguments passed on to other methods
#'
#' @returns A `tibble`, `sf`, or `SpatialPointsDataFrame` of 1D
#'   and 2D integration points, including a `weight` column, a`.block` column,
#'   and a matrix column `.block_origin`.
#'   The `.block` column is used to identify the integration
#'   blocks defined by the samplers. The `.block_origin` collects the original
#'   subdomain block information for tensor product blocks.
#' @export
#' @examples
#' # Integration on the interval (2, 3.5) with Simpson's rule
#' ips <- fm_int(fm_mesh_1d(0:4), samplers = cbind(2, 3.5))
#' plot(ips$x, ips$weight)
#'
#' # Create integration points for the two intervals [0,3] and [5,10]
#' ips <- fm_int(
#'   fm_mesh_1d(0:10),
#'   rbind(c(0, 3), c(5, 10))
#' )
#' plot(ips$x, ips$weight)
#'
#' # Convert a 1D mesh into integration points
#' mesh <- fm_mesh_1d(seq(0, 10, by = 1))
#' ips <- fm_int(mesh, name = "time")
#' plot(ips$time, ips$weight)
#'
#' if (require("ggplot2", quietly = TRUE)) {
#'   #' Integrate on a 2D mesh with polygon boundary subset
#'   ips <- fm_int(fmexample$mesh, fmexample$boundary_sf[[1]])
#'   ggplot() +
#'     geom_sf(data = fm_as_sfc(fmexample$mesh, multi = TRUE), alpha = 0.5) +
#'     geom_sf(data = fmexample$boundary_sf[[1]], fill = "red", alpha = 0.5) +
#'     geom_sf(data = ips, aes(size = weight)) +
#'     scale_size_area()
#' }
#'
fm_int <- function(domain, samplers = NULL, ...) {
  UseMethod("fm_int")
}

#' Multi-domain sampler integration
#'
#' Combine integration over different domains
#'
#' @param domain A list of named domains
#' @param samplers A named list of samplers
#' @param ... Passed on to each [fm_int()] call.
#' @param extra Optional character vector with names of variables other than the
#'   integration domains to be included from the samplers. If `NULL` (default),
#'   all additional variables are included.
#' @returns An object with integration points and weights
#' @export
#' @keywords internal
#' @examples
#' fm_int_multi_sampler(
#'   domain = list(x = fm_mesh_1d(1:4), y = 11:12),
#'   samplers = tibble::tibble(
#'     x = rbind(c(1, 3), c(2, 4)),
#'     y = c(12, 11)
#'   )
#' )
fm_int_multi_sampler <- function(domain, samplers, ..., extra = NULL) {
  if (is.null(names(domain))) {
    stop("For 'fm_int_multi_sampler', the domain must be a named list.")
  }

  names_domain <- names(domain)
  names_samplers <- names(samplers)
  names_reserved <- c("weight", ".block", ".block_origin")

  if (length(intersect(names_domain, names_reserved)) > 0) {
    stop(paste0(
      "The reserved name(s) ",
      paste0("'",
        intersect(names_domain, names_reserved),
        "'",
        collapse = ", "
      ),
      " cannot be used as domain names."
    ))
  }

  names_intersect <- intersect(names_samplers, names_domain)
  ips_list <- lapply(
    names_intersect,
    function(nm) {
      fm_int(
        domain = domain[[nm]],
        samplers = samplers[[nm]],
        name = nm,
        ...
      )
    }
  )
  ips <- do.call(fm_cprod, c(ips_list, list(.blockwise = TRUE)))

  if ("weight" %in% names_samplers) {
    ips$weight <- ips[["weight"]] * samplers[["weight"]][ips$.block]
  }

  names_extra <- setdiff(
    names_samplers,
    c(names_domain, names_reserved, names(ips))
  )
  if (!is.null(extra)) {
    names_extra <- intersect(names_extra, extra)
  }
  if (length(names_extra) > 0) {
    for (nm in names_extra) {
      ips[[nm]] <- samplers[[nm]][ips$.block]
    }
  }

  ips
}


#' @param extra Optional character vector with names of variables other than the
#'   integration domains to be included from the samplers. If `NULL` (default),
#'   all additional variables are included.
#' @export
#' @describeIn fm_int Multi-domain integration
fm_int.list <- function(domain, samplers = NULL, ..., extra = NULL) {
  weight_name <- "weight"

  if (is.null(names(domain))) {
    stop("For 'fm_int.list', the domain must be a named list.")
  }

  if (!is.null(samplers) && !inherits(samplers, "list")) {
    samplers <- list(samplers)
  }

  # Change a mix of sp, sfc, and sf objects to sf
  segm_samplers <- unlist(lapply(samplers, function(x) inherits(x, "fm_segm")))
  sfc_samplers <- unlist(lapply(samplers, function(x) inherits(x, "sfc")))
  sf_samplers <- unlist(lapply(samplers, function(x) inherits(x, "sf")))
  sp_samplers <- unlist(lapply(samplers, function(x) inherits(x, "Spatial")))
  if (any(sp_samplers)) {
    if (any(sf_samplers) || any(sfc_samplers) || any(segm_samplers)) {
      warning(paste0(
        "Both `sf` and `sp` objects in the samplers are detected.",
        " Output will be `sf`."
      ))
    }
    samplers[sp_samplers] <- lapply(samplers[sp_samplers], sf::st_as_sf)
    if (!("coordinates" %in% names(domain))) {
      stop("`sp` input detected but no `coordinates` domain present.")
    }
    names(domain)[names(domain) %in% "coordinates"] <- "geometry"
  }

  names_domain <- names(domain)
  names_lsamplers <- names(samplers)
  if (is.null(names_lsamplers)) {
    names_lsamplers <- rep("", length(samplers))
  }
  index_single_samplers <- which(names_lsamplers != "")
  index_multi_samplers <- which(names_lsamplers == "")
  names_samplers <- as.list(names_lsamplers)
  names_samplers[index_multi_samplers] <-
    lapply(
      samplers[index_multi_samplers],
      function(x) {
        if (inherits(x, c("fm_segm", "sfc"))) {
          stop(
            paste0(
              "Unnamed sampler in the samplers is an 'fm_segm' or",
              " 'sfc' object.\n",
              "  Name them explicitly, or convert to an 'sf' object with the ",
              "appropriate geometry column name,\n",
              "  or use other supported multi-sampler class instead."
            )
          )
          NULL
        } else {
          names(x)
        }
      }
    )
  # coordinate and geometry are not required here
  names_reserved <- c("weight", ".block", ".block_origin")

  if (length(intersect(names_domain, names_reserved)) > 0) {
    stop(paste0(
      "The reserved names ",
      paste0(intersect(names_domain, names_reserved), collapse = ", "),
      " cannot be used as domain names."
    ))
  }

  lips_samplers <- list()

  #######################
  # multidomain samplers, ie unnamed element(s) in samplers, for each sampler
  # and then for each domain(lapply)
  # TODO still have to deal with secondary geometry
  for (i in index_multi_samplers) {
    if (is.null(names_samplers[[i]])) {
      stop(paste0(
        "The unnamed sampler #", i, " in the samplers has no sub-names."
      ))
    }
    if (!any(names_samplers[[i]] %in% names_domain)) {
      stop(
        paste0(
          "Sampler #", i, " with names (",
          paste0(names_samplers[[i]], collapse = ","),
          ") has no matching domains (",
          paste0(names(domain), collapse = ","),
          ")."
        )
      )
    }
    lips_samplers[[i]] <-
      fm_int_multi_sampler(
        domain = domain,
        samplers = samplers[[i]],
        ...
      )
  }

  #######################
  # singledomain samplers, ie named element(s) in samplers
  for (i in index_single_samplers) {
    nm <- intersect(names_samplers[[i]], names_domain)
    if (length(nm) == 0) {
      stop(
        paste0(
          "Named sampler #", i, " (",
          names_lsamplers[[i]],
          ") has no corresponding domain (",
          paste0(names_domain, collapse = ","),
          ")."
        )
      )
    }
    lips_samplers[[i]] <-
      fm_int(
        domain = domain[[nm]],
        samplers = samplers[[i]],
        name = nm,
        ...
      )
  }

  # Full domain samplers
  names_full_domain_samplers <- setdiff(names_domain, unlist(names_samplers))
  lips_full_domain_samplers <-
    lapply(
      names_full_domain_samplers,
      function(nm) {
        fm_int(
          domain = domain[[nm]],
          name = nm,
          ...
        )
      }
    )

  ips <- do.call(fm_cprod, c(
    lips_samplers,
    lips_full_domain_samplers,
    list(.blockwise = FALSE)
  ))

  if (any(sp_samplers) && !any(sf_samplers)) {
    ips <- fm_int_object_as_Spatial(ips)
    cnames <- sp::coordnames(ips)
    sp::coordnames(ips) <- c("x", "y", "z")[seq_along(cnames)]
  }

  ips
}


# Helper for blockwise integration; used when a sampler is a list of valid
# (sub)samplers.
fm_int_block_sampler <- function(domain,
                                 sampler_row,
                                 name,
                                 ...) {
  subsampler <- sampler_row[1L, name, drop = TRUE]
  weight <- sampler_row[1L, "weight", drop = TRUE]
  block <- sampler_row[1L, ".block", drop = TRUE]
  block_origin <- sampler_row[1L, ".block_origin", drop = TRUE]

  ips_list <- list()
  for (k in seq_along(subsampler)) {
    loc_subsampler <- subsampler[[k]]
    ips_list[[k]] <- fm_int(
      domain = domain,
      samplers = loc_subsampler,
      name = name,
      ...
    )
    ips_list[[k]]$weight <- ips_list[[k]]$weight * weight
    ips_list[[k]]$.block <- block
    ips_list[[k]]$.block_origin <-
      matrix(
        block_origin,
        NROW(ips_list[[k]]),
        ncol(block_origin),
        byrow = TRUE,
        dimnames = list(
          NULL,
          colnames(block_origin)
        )
      )
  }
  dplyr::bind_rows(ips_list)
}

# Wrapper for integration over samplers that are potentially nested lists
fm_int_wrapper <- function(domain, samplers, name, ..., int_fun) {
  if (is.list(samplers[[name]])) {
    ips <- list()
    for (j in seq_along(samplers[[name]])) {
      subsampler <- samplers[j, , drop = TRUE]
      if (is.list(subsampler)) {
        ips_list <- fm_int_block_sampler(
          domain,
          sampler_row = samplers[j, , drop = FALSE],
          name = name,
          ...
        )
        ips[[j]] <- dplyr::bind_rows(ips_list)
      } else {
        ips[[j]] <- int_fun(
          domain,
          samplers = samplers[j, , drop = FALSE],
          name = name,
          ...
        )
      }
    }
    ips <- new_fm_int(
      do.call(dplyr::bind_rows, ips),
      name = name
    )
  } else {
    ips <- int_fun(domain, samplers = samplers, name = name, ...)
  }

  if (NROW(ips) == 0) {
    ips <- new_fm_int(tibble::tibble(
      "{name}" := numeric(0),
      weight = numeric(0),
      .block = integer(0),
      .block_origin = matrix(
        integer(0),
        nrow = 0,
        ncol = NCOL(ips[[".block_origin"]]),
        dimnames = list(NULL, colnames(ips[[".block_origin"]]))
      )
    ))
  }

  ips
}


#' @export
#' @describeIn fm_int Discrete double or integer space integration
#' @examples
#' # Individual sampling points:
#' (ips <- fm_int(0:10, c(0, 3, 5, 6, 10)))
#' # Sampling blocks:
#' (ips <- fm_int(0:10, list(c(0, 3), c(5, 6, 10))))
#'
fm_int.numeric <- function(domain, samplers = NULL, name = "x", ...) {
  if (is.null(samplers)) {
    ips <- new_fm_int(as.vector(domain), name = name)
    return(ips)
  }

  samplers <- new_fm_int(samplers, blocks = TRUE, name = name)

  fm_int_numeric <- function(domain, samplers, name, ...) {
    storage.mode(samplers[[name]]) <- storage.mode(domain)

    ok <- samplers[[name]] %in% domain
    ips <- samplers[ok, , drop = FALSE]
    ips
  }

  ips <- fm_int_wrapper(
    domain = domain, samplers = samplers, name = name, ...,
    int_fun = fm_int_numeric
  )

  ips
}

#' @export
#' @describeIn fm_int Discrete character space integration
fm_int.character <- function(domain, samplers = NULL, name = "x", ...) {
  if (is.null(samplers)) {
    ips <- new_fm_int(as.vector(domain), name = name)
    return(ips)
  }

  samplers <- new_fm_int(samplers, blocks = TRUE, name = name)

  fm_int_character <- function(domain, samplers, name, ...) {
    storage.mode(samplers[[name]]) <- storage.mode(domain)

    ok <- samplers[[name]] %in% domain
    ips <- samplers[ok, , drop = FALSE]
    ips
  }

  ips <- fm_int_wrapper(
    domain = domain, samplers = samplers, name = name, ...,
    int_fun = fm_int_character
  )

  ips
}

#' @export
#' @describeIn fm_int Discrete factor space integration
fm_int.factor <- function(domain, samplers = NULL, name = "x", ...) {
  if (is.null(samplers)) {
    ips <- new_fm_int(as.vector(domain), name = name)
    return(ips)
  }

  samplers <- new_fm_int(samplers, blocks = TRUE, name = name)

  fm_int_factor <- function(domain, samplers, name, ...) {
    samplers[[name]] <- factor(as.vector(samplers[[name]]),
      levels = levels(domain)
    )

    ok <- samplers[[name]] %in% domain
    ips <- samplers[ok, , drop = FALSE]
    ips
  }

  ips <- fm_int_wrapper(
    domain = domain, samplers = samplers, name = name, ...,
    int_fun = fm_int_factor
  )

  ips
}


#' @export
#' @describeIn fm_int `SpatRaster` integration. Not yet implemented.
fm_int.SpatRaster <- function(domain, samplers = NULL, name = "x", ...) {
  stop("'SpatRaster' integration is not yet implemented.")
}

#' @export
#' @describeIn fm_int `fm_lattice_2d` integration. Not yet implemented.
fm_int.fm_lattice_2d <- function(domain, samplers = NULL, name = "x", ...) {
  stop("'fm_lattice_2d' integration is not yet implemented.")
}


# fm_mesh_1d integration ####

#' @param int.args List of arguments passed to line and integration methods.
#' * `method`: "stable" (to aggregate integration weights onto mesh nodes)
#'   or "direct" (to construct a within triangle/segment integration scheme
#'   without aggregating onto mesh nodes)
#' * `nsub1`, `nsub2`: integers controlling the number of internal integration
#'   points before aggregation. Points per triangle: `(nsub2+1)^2`.
#'   Points per knot segment: `nsub1`
#' @export
#' @describeIn fm_int `fm_mesh_1d` integration. Supported samplers:
#' * `NULL` for integration over the entire domain;
#' * A vector defining points for summation (up to `0.5.0`, length 2 vectors
#'   were interpreted as intervals. From 0.6.0 intervals must be specified as
#'   rows of a 2-column matrix);
#' * A 2-column matrix with a single interval in each row;
#' * A list of such vectors or matrices
#' * A tibble with a named column containing a vector/matrix/list as above,
#'   and optionally a `weight` column.
#' @examples
#' # Continuous integration on intervals
#' ips <- fm_int(
#'   fm_mesh_1d(0:10, boundary = "cyclic"),
#'   rbind(c(0, 3), c(5, 10))
#' )
#' plot(ips$x, ips$weight)
#'
fm_int.fm_mesh_1d <- function(domain,
                              samplers = NULL,
                              name = "x",
                              int.args = NULL,
                              format = NULL,
                              ...) {
  format <- match.arg(format, c("numeric", "bary"))

  int.args.default <- list(method = "stable", nsub1 = 30, nsub2 = 9)
  if (is.null(int.args)) {
    int.args <- list()
  }
  missing.args <- setdiff(names(int.args.default), names(int.args))
  int.args[missing.args] <- int.args.default[missing.args]
  if (!is.null(int.args[["nsub"]])) {
    int.args[["nsub1"]] <- int.args[["nsub"]]
  }

  if (is.null(samplers)) {
    samplers <- new_fm_int(
      cbind(domain$interval[1], domain$interval[2]),
      blocks = FALSE,
      name = name
    )
  } else {
    samplers <- new_fm_int(samplers, blocks = TRUE, name = name)
    if (is.data.frame(samplers) && !(name %in% colnames(samplers))) {
      stop(paste0("Domain name '", name, "' missing from `samplers`."))
    }
  }

  fm_int_mesh_1d <- function(domain, samplers, name, int.args, format, ...) {
    fm_int_1d_points <- function(domain, sampler_row) {
      subsampler <- sampler_row[1L, name, drop = TRUE]
      theweight <- sampler_row[1L, "weight", drop = TRUE]
      the.block <- sampler_row[1L, ".block", drop = TRUE]
      the.block_origin <- sampler_row[1L, ".block_origin", drop = TRUE]

      # Pointwise summation
      loc_point <- subsampler
      weight_point <- rep(1.0, length(loc_point))

      ips <- new_fm_int(tibble::tibble(
        "{name}" := loc_point[weight_point > 0],
        weight = weight_point[weight_point > 0] * theweight,
        .block = the.block,
        .block_origin = matrix(
          the.block_origin,
          sum(weight_point > 0),
          ncol(the.block_origin),
          byrow = TRUE,
          dimnames = list(NULL, colnames(the.block_origin))
        )
      ))
      ips
    }

    fm_int_1d_interval <- function(domain, sampler_row) {
      subsampler <- sampler_row[1L, name, drop = TRUE]
      theweight <- sampler_row[1L, "weight", drop = TRUE]
      the.block <- sampler_row[1L, ".block", drop = TRUE]
      the.block_origin <- sampler_row[1L, ".block_origin", drop = TRUE]

      # Interval integration
      subsampler <- as.vector(subsampler)
      if (isTRUE(domain$cyclic)) {
        if (diff(subsampler) >= diff(domain$interval)) {
          subsampler <- domain$interval
        } else {
          subsampler[1] <- domain$interval[1] +
            (subsampler[1] - domain$interval[1]) %% diff(domain$interval)
          subsampler[2] <- subsampler[1] +
            diff(subsampler) %% diff(domain$interval)
          if (diff(subsampler) == 0.0) {
            subsampler <- domain$interval
          } else if (subsampler[2] > domain$interval[2]) {
            subsampler[2] <- domain$interval[1] +
              (subsampler[2] - domain$interval[1]) %% diff(domain$interval)
          }
        }
      } else if (diff(subsampler) <= 0.0) {
        # Empty interval, skip to next subsampler
        return(NULL)
      }

      if (identical(int.args[["method"]], "stable")) {
        if (isTRUE(domain$cyclic)) {
          loc_trap <- sort(unique(pmin(
            domain$interval[2],
            pmax(
              domain$interval[1],
              c(
                domain$loc,
                as.vector(subsampler),
                domain$interval
              )
            )
          )))
        } else {
          loc_trap <- sort(unique(c(
            domain$interval,
            domain$loc,
            as.vector(subsampler)
          )))
        }

        # Simpson's rule integration
        loc_mid <- (loc_trap[-1] + loc_trap[-length(loc_trap)]) / 2
        # Detect mid-points inside the samplers
        if (isTRUE(domain$cyclic) && (subsampler[1] > subsampler[2])) {
          inside <- (loc_mid < min(subsampler)) |
            (loc_mid > max(subsampler))
        } else {
          inside <- (loc_mid >= min(subsampler)) &
            (loc_mid <= max(subsampler))
        }
        weight_mid <- diff(loc_trap)
        weight_mid[!inside] <- 0.0

        weight_trap <- c(weight_mid / 2, 0) + c(0, weight_mid / 2)
        loc_simpson <- c(loc_trap, loc_mid)
        weight_simpson <- c(weight_trap / 3, weight_mid * 2 / 3)

        ips <- new_fm_int(tibble::tibble(
          "{name}" := loc_simpson[(weight_simpson > 0)],
          weight = weight_simpson[(weight_simpson > 0)] * theweight,
          .block = the.block,
          .block_origin = matrix(
            the.block_origin,
            sum(weight_simpson > 0),
            ncol(the.block_origin),
            byrow = TRUE,
            dimnames = list(NULL, colnames(the.block_origin))
          )
        ))
      } else {
        nsub <- int.args[["nsub1"]]
        u <- rep(
          (seq_len(nsub) - 0.5) / nsub,
          domain$n - 1
        )
        int_loc <-
          domain$loc[rep(seq_len(domain$n - 1), each = nsub)] * (1 - u) +
          domain$loc[rep(seq_len(domain$n - 1) + 1, each = nsub)] * u
        int_w <-
          (domain$loc[rep(seq_len(domain$n - 1) + 1, each = nsub)] -
            domain$loc[rep(seq_len(domain$n - 1), each = nsub)]) /
            nsub

        if (isTRUE(domain$cyclic) && (subsampler[1] > subsampler[2])) {
          inside <- (int_loc < min(subsampler)) |
            (int_loc > max(subsampler))
        } else {
          inside <- (int_loc >= min(subsampler)) &
            (int_loc <= max(subsampler))
        }

        ips <- new_fm_int(tibble::tibble(
          "{name}" := int_loc[inside],
          weight = int_w[inside] * theweight,
          .block = the.block,
          .block_origin = matrix(the.block_origin,
            sum(inside),
            ncol(the.block_origin),
            byrow = TRUE,
            dimnames = list(NULL, name)
          )
        ))
      }

      ips
    }

    ips <- list()
    for (j in seq_len(NROW(samplers))) {
      subsampler <- samplers[j, name, drop = TRUE]

      if (!is.matrix(subsampler)) {
        ips[[j]] <- fm_int_1d_points(
          domain = domain,
          sampler_row = samplers[j, , drop = FALSE]
        )
      } else {
        ips[[j]] <- fm_int_1d_interval(
          domain = domain,
          sampler_row = samplers[j, , drop = FALSE]
        )
      }
    }

    ips <- new_fm_int(
      do.call(dplyr::bind_rows, ips),
      name = name
    )

    if (NROW(ips) == 0) {
      ips <- new_fm_int(tibble::tibble(
        "{name}" := numeric(0),
        weight = numeric(0),
        .block = integer(0),
        .block_origin = matrix(
          integer(0),
          nrow = 0,
          ncol = NCOL(ips[[".block_origin"]]),
          dimnames = list(NULL, colnames(ips[[".block_origin"]]))
        )
      ))
    }

    if (identical(format, "bary")) {
      # TODO: Reverse the logic above, and construct barycentric coordinates
      # directly
      ips[[name]] <- fm_bary(domain, ips[[name]])
    }
    ips
  }

  ips <- fm_int_wrapper(
    domain = domain, samplers = samplers, name = name,
    int.args = int.args, format = format,
    int_fun = fm_int_mesh_1d
  )

  ips
}


# fm_mesh_2d integration ####

#' @export
#' @describeIn fm_int `fm_mesh_2d` integration. Any sampler class with an
#' associated [fm_int_mesh_2d()] method is supported.
#' @param format character; determines the output format, as either "sf"
#'   (default for `fm_mesh_2d` when the sampler is `NULL`),
#'   "numeric" (default for `fm_mesh_1d`), "bary", or "sp".
#'   When `NULL`, determined by the domain and sampler types.
fm_int.fm_mesh_2d <- function(domain,
                              samplers = NULL,
                              name = NULL,
                              int.args = NULL,
                              format = NULL,
                              ...) {
  int.args.default <- list(method = "stable", nsub1 = 30, nsub2 = 9)
  if (is.null(int.args)) {
    int.args <- list()
  }
  missing.args <- setdiff(names(int.args.default), names(int.args))
  int.args[missing.args] <- int.args.default[missing.args]
  if (!is.null(int.args[["nsub"]])) {
    int.args[["nsub2"]] <- int.args[["nsub"]]
  }

  ips <- fm_int_mesh_2d(
    samplers,
    domain = domain,
    name = name,
    int.args = int.args,
    ...
  )

  if (is.null(format) && inherits(samplers, "Spatial")) {
    format <- "sp"
  }
  if (!is.null(format)) {
    if (identical(format, "bary")) {
      # TODO: Reverse the logic of fm_int_mesh_2d() to generate fm_bary directly
      ips <- fm_bary(domain, ips)
    } else if (identical(format, "sf") && !inherits(ips, "sf")) {
      ips <- sf::st_as_sf(ips)
      if (!is.null(name) && (name != attr(ips, "sf_column"))) {
        ips <- dplyr::rename(ips, "{name}" := attr(ips, "sf_column"))
      }
    } else if (identical(format, "sp") && !inherits(ips, "Spatial")) {
      ips <- fm_int_object_as_Spatial(ips)
      cnames <- sp::coordnames(ips)
      sp::coordnames(ips) <- c("x", "y", "z")[seq_along(cnames)]
    }
  }

  ips
}

#' @title Project integration points to mesh vertices
#'
#' @description
#' Compute information for assigning points to the vertices of the covering
#' triangle
#'
#' @param points A `SpatialPointsDataFrame`, `sf`, `tibble`, or `list` object
#' @param mesh An `fm_mesh_2d` object
#' @returns `SpatialPointsDataFrame`, `sf`, `tibble`, or `list` of mesh
#' vertices with projected data attached
#' @importFrom rlang .data
#' @keywords internal
#' @export
#' @examples
#' head(fm_vertex_projection(list(loc = fmexample$loc), fmexample$mesh))
#' head(fm_vertex_projection(fmexample$loc_sf, fmexample$mesh))
#'
fm_vertex_projection <- function(points, mesh) {
  if (inherits(points, c("sf", "sfc")) ||
    inherits(points, "Spatial")) {
    n_points <- NROW(points)
    res <- fm_bary(mesh, points)
  } else if (inherits(points, "fm_bary")) {
    n_points <- NROW(points)
    res <- points
  } else {
    n_points <- NROW(points$loc)
    res <- fm_bary(mesh, points$loc)
  }
  tri <- res$index
  bary <- res$where

  if (inherits(points, "sfc")) {
    points <- sf::st_sf(geometry = points)
  }

  if (is.null(points[["weight"]])) {
    points$weight <- rep(1L, n_points)
  }
  if (is.null(points[[".block"]])) {
    points$.block <- rep(1L, n_points)
  }
  if (is.null(points[[".block_origin"]])) {
    points$.block_origin <- matrix(
      points$.block,
      n_points,
      1L,
      dimnames = NULL
    )
  }

  ok <- !is.na(tri)
  ok[ok] <- (tri[ok] > 0)
  if (any(!ok)) {
    warning(
      paste0(
        "Some integration points were outside the mesh;",
        " check your coordinate systems."
      )
    )
  }

  data <-
    new_fm_int(tibble::tibble(
      .vertex = as.vector(mesh$graph$tv[tri[ok], ]),
      weight = as.vector(points$weight[ok] * bary[ok, ]),
      .block = rep(points$.block[ok], times = 3),
      .block_origin = rbind(
        points$.block_origin[ok, , drop = FALSE],
        points$.block_origin[ok, , drop = FALSE],
        points$.block_origin[ok, , drop = FALSE]
      )
    ))

  data <-
    dplyr::summarise(
      dplyr::group_by(data, .data$.vertex, .data$.block, .data$.block_origin),
      weight = sum(.data$weight),
      .groups = "drop"
    )
  coords <- mesh$loc[data$.vertex, , drop = FALSE]
  data <- dplyr::select(data, c("weight", ".block", ".block_origin", ".vertex"))

  if (inherits(points, "Spatial")) {
    fm_safe_sp(force = TRUE)
    # sp cannot handle matrix columns, so we convert to a plain vector
    data$.block_origin <- data$.block
    ret <- sp::SpatialPointsDataFrame(
      coords[, seq_len(min(
        ncol(coords),
        ncol(sp::coordinates(points))
      )), drop = FALSE],
      proj4string = fm_CRS(mesh),
      data = data,
      match.ID = FALSE
    )
    sp::coordnames(ret) <- sp::coordnames(points)
  } else if (inherits(points, "sf")) {
    colnames(coords) <- c("X", "Y", "Z")[seq_len(ncol(coords))]
    d <- min(ncol(coords), length(intersect(
      colnames(sf::st_coordinates(points)), c("X", "Y", "Z")
    )))
    data <- dplyr::bind_cols(
      tibble::as_tibble(coords[, seq_len(d), drop = FALSE]),
      tibble::as_tibble(data)
    )
    ret <- sf::st_as_sf(
      data,
      coords = seq_len(d),
      crs = fm_crs(mesh)
    )
  } else {
    colnames(coords) <- c("X", "Y", "Z")[seq_len(ncol(coords))]
    data$loc <- coords
    ret <- data
  }

  if (inherits(ret, "sf") && inherits(points, "sf")) {
    if (!is.null(attr(points, "sf_column")) &&
      (attr(points, "sf_column") != attr(ret, "sf_column"))) {
      ret <- dplyr::rename(ret, "{attr(points, 'sf_column')}" := "geometry")
    }
  }

  ret
}

#' Subset integration on a mesh
#'
#' Integration methods for spatial samplers on `fm_mesh_2d` meshes.
#'
#' @returns A `list`, `sf`, or `Spatial` object with
#' point coordinate information and additional columns `weight` and `.block`
#' @inheritParams fm_int
#' @export
#' @keywords internal
#' @examples
#' str(fm_int_mesh_2d(samplers = NULL, domain = fmexample$mesh))
fm_int_mesh_2d <- function(samplers,
                           domain,
                           name = NULL,
                           int.args = NULL,
                           ...) {
  stopifnot(inherits(domain, "fm_mesh_2d"))

  if (missing(samplers) || is.null(samplers)) {
    return(
      fm_int_mesh_2d_NULL(
        samplers = NULL,
        domain = domain,
        name = name,
        int.args = int.args,
        ...
      )
    )
  }

  UseMethod("fm_int_mesh_2d")
}

#' @describeIn fm_int_mesh_2d Full domain integration
fm_int_mesh_2d_NULL <- function(samplers,
                                domain,
                                name = NULL,
                                int.args = NULL,
                                ...) {
  stopifnot(is.null(samplers))

  ips <- fm_int_mesh_2d_polygon(
    domain = domain,
    samplers = NULL,
    int.args = int.args,
    name = name
  )

  if (!is.null(name) && (name != attr(ips, "sf_column"))) {
    ips <- dplyr::rename(ips, "{name}" := "geometry")
  }

  ips
}

#' @export
#' @describeIn fm_int_mesh_2d `sf` integration
fm_int_mesh_2d.sf <- function(samplers,
                              domain,
                              name = NULL,
                              int.args = NULL,
                              ...) {
  if (is.null(name)) {
    name <- attr(samplers, "sf_column")
  }
  if (!("weight" %in% names(samplers))) {
    weight <- rep(1, NROW(samplers))
  } else {
    weight <- samplers$weight
  }

  fm_int_mesh_2d(
    sf::st_geometry(samplers),
    domain,
    name = name,
    int.args = int.args,
    .weight = weight,
    ...
  )
}

#' @export
#' @describeIn fm_int_mesh_2d `sfc_POINT` integration
fm_int_mesh_2d.sfc_POINT <- function(samplers,
                                     domain,
                                     name = NULL,
                                     int.args = NULL,
                                     .weight = rep(1, NROW(samplers)),
                                     ...) {
  if (is.null(name)) {
    name <- "geometry"
  }
  ips <- new_fm_int(
    samplers,
    blocks = TRUE,
    weight = .weight,
    name = name,
    override = TRUE
  )

  # TODO: remove points outside the domain

  ips
}


#' @export
#' @describeIn fm_int_mesh_2d `sfc_MULTIPOINT` integration
#' @importFrom rlang :=
fm_int_mesh_2d.sfc_MULTIPOINT <- function(samplers,
                                          domain,
                                          name = NULL,
                                          int.args = NULL,
                                          .weight = rep(1, NROW(samplers)),
                                          ...) {
  coords <- tibble::as_tibble(sf::st_coordinates(samplers))
  coords <- dplyr::rename(coords, .block = "L1")
  coords$weight <- .weight[coords$.block]
  ips <- sf::st_as_sf(
    coords,
    coords = intersect(names(coords), c("X", "Y", "Z", "M")),
    crs = fm_crs(samplers)
  )

  ips <- new_fm_int(
    ips,
    name = name,
    override = TRUE
  )

  # TODO: remove points outside the domain

  ips
}


fm_int_mesh_2d_lines <- function(samplers,
                                 domain,
                                 name = NULL,
                                 int.args = NULL,
                                 .weight = rep(1, NROW(samplers)),
                                 ...) {
  project <- identical(int.args$method, "stable")

  weight <- .weight
  .block <- seq_len(NROW(samplers))
  .block_origin <- matrix(.block, NROW(samplers), 1,
    dimnames = list(NULL, name)
  )

  # Extract start and end coordinates
  coords <- sf::st_coordinates(samplers)
  if ("L2" %in% colnames(coords)) {
    # MULTILINESTRING
    feature <- coords[, "L2"]
    part <- coords[, "L1"]
    L_idx <- which(colnames(coords) %in% c("L1", "L2"))
  } else {
    feature <- coords[, "L1"]
    part <- rep(1, nrow(coords))
    L_idx <- which(colnames(coords) %in% "L1")
  }

  segment <- which((diff(part) > 0) | (diff(feature) > 0))
  coordnames <- intersect(colnames(coords), c("X", "Y", "Z", "M"))

  sp <- coords[-c(segment, nrow(coords)), coordnames, drop = FALSE]
  ep <- coords[-c(1L, 1L + segment), coordnames, drop = FALSE]
  origin <- feature[-c(segment, nrow(coords))]
  segm <- fm_segm(
    rbind(sp, ep),
    idx = cbind(seq_len(nrow(sp)), seq_len(nrow(sp)) + nrow(sp)),
    grp = origin,
    crs = fm_crs(samplers)
  )

  sampler_crs <- fm_crs(samplers)
  target_crs <- fm_crs(domain)
  if (!fm_crs_is_null(sampler_crs) &&
    fm_crs_is_null(target_crs)) {
    target_crs <- sampler_crs
  }

  # Filter out points outside the mesh...
  segm <- fm_transform(segm, crs = target_crs, passthrough = TRUE)

  # Split at mesh edges
  segm <- fm_split_lines(domain, segm)
  origin <- origin[segm$origin]

  # At this point, segm is in the target_crs

  # Determine integration points along lines

  if (fm_crs_is_null(sampler_crs) || !fm_crs_is_geocent(fm_crs(domain))) {
    sp <- segm$loc[segm$idx[, 1], , drop = FALSE]
    ep <- segm$loc[segm$idx[, 2], , drop = FALSE]
    if (fm_manifold(domain, "S2")) {
      sp <- sp / rowSums(sp^2)^0.5
      ep <- ep / rowSums(ep^2)^0.5
      ips <- (sp + ep) / 2
      w <- rowSums((ep - sp)^2)^0.5

      radius <- mean(rowSums(domain$loc^2)^0.5)
      w <- asin(pmin(1, pmax(-1, w / 2))) * radius
      ips <- ips * (radius / rowSums(ips^2)^0.5)
    } else {
      ips <- (sp + ep) / 2
      w <- rowSums((ep - sp)^2)^0.5
    }
  } else {
    # Has CRS
    longlat.crs <- fm_crs("longlat_globe")
    geocentric.crs <- fm_crs("sphere")
    segm3d <- fm_transform(segm, crs = geocentric.crs, crs0 = target_crs)
    sp3d <- segm$loc[segm$idx[, 1], , drop = FALSE]
    ep3d <- segm$loc[segm$idx[, 2], , drop = FALSE]
    mp3d <- (sp3d + ep3d) / rowSums((sp3d + ep3d)^2)^0.5

    ips <- fm_transform(mp3d, crs = target_crs, crs0 = geocentric.crs)
    w <- sp::spDists(
      fm_transform(sp3d,
        crs = longlat.crs,
        crs0 = geocentric.crs
      )[, 1:2, drop = FALSE],
      fm_transform(ep3d,
        crs = longlat.crs,
        crs0 = geocentric.crs
      )[, 1:2, drop = FALSE],
      diagonal = TRUE,
      longlat = TRUE
    )
  }

  # Wrap everything up and perform projection according to distance and given
  # .block argument
  ips <- data.frame(ips)
  d_ips <- ncol(ips)
  # Temporary names
  colnames(ips) <- c("x", "y", "z")[seq_len(d_ips)]

  ips <- sf::st_as_sf(
    tibble::as_tibble(ips),
    coords = seq_len(d_ips),
    crs = target_crs
  )
  ips$weight <- w * weight[origin]
  ips$.block <- .block[origin]
  ips$.block_origin <- matrix(
    .block_origin[origin, , drop = FALSE],
    nrow = NROW(ips),
    ncol = 1,
    byrow = TRUE,
    dimnames = list(NULL, colnames(.block_origin))
  )
  ips <- new_fm_int(ips, name = name, override = TRUE)

  # Project to mesh vertices
  if (project) {
    ips <- fm_vertex_projection(ips, domain)
  }

  ips
}


#' @export
#' @describeIn fm_int_mesh_2d `sfc_LINESTRING` integration
fm_int_mesh_2d.sfc_LINESTRING <- function(samplers,
                                          domain,
                                          name = NULL,
                                          int.args = NULL,
                                          .weight = rep(1, NROW(samplers)),
                                          ...) {
  ips <- fm_int_mesh_2d_lines(samplers, domain, name, int.args, .weight, ...)

  if (!is.null(name) && (name != attr(ips, "sf_column"))) {
    ips <- dplyr::rename(ips, "{name}" := "geometry")
  }

  ips
}

#' @export
#' @describeIn fm_int_mesh_2d `sfc_MULTILINESTRING` integration
fm_int_mesh_2d.sfc_MULTILINESTRING <- function(samplers,
                                               domain,
                                               name = NULL,
                                               int.args = NULL,
                                               .weight = rep(1, NROW(samplers)),
                                               ...) {
  ips <- fm_int_mesh_2d_lines(samplers, domain, name, int.args, .weight, ...)

  if (!is.null(name) && (name != attr(ips, "sf_column"))) {
    ips <- dplyr::rename(ips, "{name}" := "geometry")
  }

  ips
}


#' Integration scheme for mesh triangle interiors
#'
#' @param mesh Mesh on which to integrate
#' @param tri_subset Optional triangle index vector for integration on a subset
#' of the mesh triangles (Default `NULL`)
#' @param nsub number of subdivision points along each triangle edge, giving
#'    `(nsub + 1)^2` proto-integration points used to compute
#'   the vertex weights
#'   (default `nsub=9`, giving 100 integration points for each triangle)
#' @returns `tibble` with columns `loc` and `weight` with
#'   integration points for the mesh
#' @author Finn Lindgren <Finn.Lindgren@@gmail.com>
#' @keywords internal
#' @export
#' @examples
#' str(fm_int_mesh_2d_core(fmexample$mesh))
#'
fm_int_mesh_2d_core <- function(mesh, tri_subset = NULL, nsub = NULL) {
  # Construct a barycentric grid of subdivision triangle midpoints
  if (is.null(nsub)) {
    nsub <- 9
  }
  stopifnot(nsub >= 0)
  nB <- (nsub + 1)^2

  nT <- nrow(mesh$graph$tv)
  if (is.null(tri_subset)) {
    tri_subset <- seq_len(nT)
  }

  is_spherical <- fm_manifold(mesh, "S2")

  # Barycentric integration coordinates
  b <- seq(1 / 3, 1 / 3 + nsub, length.out = nsub + 1) / (nsub + 1)
  bb <- as.matrix(expand.grid(b, b))
  # Points above the diagonal should be reflected into the lower triangle:
  refl <- rowSums(bb) > 1
  if (any(refl)) {
    bb[refl, ] <- cbind(1 - bb[refl, 2], 1 - bb[refl, 1])
  }
  # Construct complete barycentric coordinates:
  barycentric_grid <- cbind(1 - rowSums(bb), bb)

  # Construct integration points
  loc <- matrix(0.0, length(tri_subset) * nB, ncol(mesh$loc))
  idx_end <- 0
  for (tri in tri_subset) {
    idx_start <- idx_end + 1
    idx_end <- idx_start + nB - 1
    loc[seq(idx_start, idx_end, length.out = nB), ] <-
      as.matrix(barycentric_grid %*%
        mesh$loc[mesh$graph$tv[tri, ], , drop = FALSE])
  }

  if (is_spherical) {
    # Normalise
    radius <- sum(mesh$loc[1, ]^2)^0.5
    mesh$loc <- mesh$loc / radius
    loc <- loc / rowSums(loc^2)^0.5
  }

  # Construct integration weights
  tri_area <- fm_fem(mesh, order = 1)$ta[tri_subset]

  if (is_spherical) {
    tri_area <- tri_area * radius^2
    loc <- loc * radius
  }

  tibble::tibble(
    loc = loc,
    weight = rep(tri_area / nB, each = nB)
  )
}


fm_int_mesh_2d_polygon <- function(samplers,
                                   domain,
                                   name = NULL,
                                   int.args = NULL,
                                   .weight = rep(1, NROW(samplers)),
                                   ...) {
  method <- match.arg(int.args[["method"]], c("stable", "direct"))

  ipsl <- list()

  # Compute direct integration points
  # TODO: Allow blockwise construction to avoid
  # overly large temporary coordinate matrices (via tri_subset)
  integ <- fm_int_mesh_2d_core(domain, nsub = int.args[["nsub2"]])

  # Keep points with positive weights (This should be all,
  # but if there's a degenerate triangle, this gets rid of it)
  ok <- (integ$weight > 0)
  integ <- integ[ok, , drop = FALSE]

  domain_crs <- fm_crs(domain)

  if (!is.null(samplers)) {
    samplers_crs <- fm_crs(samplers)
    integ_sf <- sf::st_as_sf(
      as.data.frame(integ$loc),
      coords = seq_len(ncol(integ$loc)),
      crs = domain_crs
    )
    if (!identical(domain_crs, samplers_crs) &&
      !fm_crs_is_null(domain_crs) &&
      !fm_crs_is_null(samplers_crs)) {
      integ_sf <- fm_transform(integ_sf,
        crs = samplers_crs,
        passthrough = TRUE
      )
    }

    idx <- sf::st_contains(samplers, integ_sf, sparse = TRUE)

    if (method %in% c("stable")) {
      integ_bary_ <- fm_bary(domain, integ_sf)
      integ_bary_$weight <- integ$weight
    }

    for (g in seq_along(idx)) {
      if (length(idx[[g]]) > 0) {
        integ_ <- integ[idx[[g]], , drop = FALSE]

        if (method %in% c("stable")) {
          integ_ <- integ_bary_[idx[[g]], , drop = FALSE]
          # Project integration points and weights to mesh nodes
          integ_ <- fm_vertex_projection(integ_, domain)
        }

        if (ncol(integ_$loc) > 2) {
          ips <- new_fm_int(
            sf::st_as_sf(
              tibble::tibble(
                x = integ_$loc[, 1],
                y = integ_$loc[, 2],
                z = integ_$loc[, 3],
                weight = integ_$weight,
                .block = g
              ),
              coords = c("x", "y", "z"),
              crs = domain_crs
            ),
            name = name,
            override = TRUE
          )
        } else {
          ips <- new_fm_int(
            sf::st_as_sf(
              tibble::tibble(
                x = integ_$loc[, 1],
                y = integ_$loc[, 2],
                weight = integ_$weight,
                .block = g
              ),
              coords = c("x", "y"),
              crs = domain_crs
            ),
            name = name,
            override = TRUE
          )
        }

        ipsl <- c(ipsl, list(ips))
      }
    }
  } else {
    if (method %in% c("stable")) {
      # Project integration points and weights to mesh nodes
      integ <- fm_vertex_projection(integ, domain)
    }

    if (ncol(integ$loc) > 2) {
      ipsl <- list(new_fm_int(
        sf::st_as_sf(
          tibble::tibble(
            x = integ$loc[, 1],
            y = integ$loc[, 2],
            z = integ$loc[, 3],
            weight = integ$weight,
            .block = 1L
          ),
          coords = c("x", "y", "z"),
          crs = domain_crs
        ),
        name = name,
        override = TRUE
      ))
    } else {
      ipsl <- list(new_fm_int(
        sf::st_as_sf(
          tibble::tibble(
            x = integ$loc[, 1],
            y = integ$loc[, 2],
            weight = integ$weight,
            .block = 1L
          ),
          coords = c("x", "y"),
          crs = domain_crs
        ),
        name = name,
        override = TRUE
      ))
    }
  }

  ips <- new_fm_int(
    do.call(dplyr::bind_rows, ipsl),
    name = name,
    override = TRUE
  )

  ips
}


#' @export
#' @describeIn fm_int_mesh_2d `sfc_POLYGON` integration
fm_int_mesh_2d.sfc_POLYGON <- function(samplers,
                                       domain,
                                       name = NULL,
                                       int.args = NULL,
                                       .weight = rep(1, NROW(samplers)),
                                       ...) {
  weight <- .weight
  .block <- seq_len(NROW(samplers))

  ips <- fm_int_mesh_2d_polygon(
    domain = domain,
    int.args = int.args,
    samplers = samplers,
    name = name
  )

  ips$weight <- ips$weight * .weight[ips$.block]
  ips$.block <- .block[ips$.block]
  ips <- new_fm_int(ips, name = name, override = TRUE)

  ips
}

#' @export
#' @describeIn fm_int_mesh_2d `sfc_MULTIPOLYGON` integration
fm_int_mesh_2d.sfc_MULTIPOLYGON <- function(samplers,
                                            domain,
                                            name = NULL,
                                            int.args = NULL,
                                            .weight = rep(1, NROW(samplers)),
                                            ...) {
  weight <- .weight
  .block <- seq_len(NROW(samplers))

  ips <- fm_int_mesh_2d_polygon(
    domain = domain,
    int.args = int.args,
    samplers = samplers,
    name = name
  )

  ips$weight <- ips$weight * .weight[ips$.block]
  ips$.block <- .block[ips$.block]
  ips <- new_fm_int(ips, name = name, override = TRUE)

  ips
}


#' @export
#' @describeIn fm_int_mesh_2d `sfc_GEOMERY` integration
fm_int_mesh_2d.sfc_GEOMETRY <- function(samplers,
                                        domain,
                                        name = NULL,
                                        int.args = NULL,
                                        .weight = rep(1, NROW(samplers)),
                                        ...) {
  geometry_class <- vapply(
    seq_along(samplers),
    function(x) {
      class(samplers[x])[1]
    },
    ""
  )
  .block <- seq_len(NROW(samplers))

  ips <- list()
  for (g_class in unique(geometry_class)) {
    subset <- geometry_class == g_class
    ips[[g_class]] <-
      fm_int_mesh_2d(samplers[subset],
        domain = domain,
        name = name,
        int.args = int.args,
        .weight = .weight[subset]
      )
    ips[[g_class]][[".block"]] <- .block[subset][ips[[g_class]][[".block"]]]
    ips[[g_class]] <- new_fm_int(ips[[g_class]], name = name)
  }
  ips <- new_fm_int(
    do.call(dplyr::bind_rows, ips),
    name = name,
    override = TRUE
  )

  ips
}


#' @export
#' @describeIn fm_int_mesh_2d `Spatial` integration
fm_int_mesh_2d.Spatial <- function(samplers,
                                   domain,
                                   name = NULL,
                                   int.args = NULL,
                                   format = NULL,
                                   ...) {
  samplers <- sf::st_as_sf(samplers)

  ips <-
    fm_int_mesh_2d(
      samplers,
      domain = domain,
      name = name,
      int.args = int.args,
      ...
    )

  ips
}

#' @export
#' @describeIn fm_int_mesh_2d `fm_segm` integration
fm_int_mesh_2d.fm_segm <- function(samplers,
                                   domain,
                                   name = NULL,
                                   int.args = NULL,
                                   format = NULL,
                                   ...) {
  samplers <- fm_as_sfc(samplers)

  ips <-
    fm_int_mesh_2d(
      samplers,
      domain = domain,
      name = name,
      int.args = int.args,
      ...
    )

  ips
}
