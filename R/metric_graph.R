# inlabru MAPPER ----

# required packages/suggest packages
# @rawNamespace S3method(inlabru::bru_get_mapper,rspde_metric_graph)
# @rawNamespace S3method(inlabru::bru_mapper, metric_graph)
# @rawNamespace S3method(inlabru::bru_mapper, fm_MG)

#' @title Wrapper that calls bru_mapper with correct input
#' @description Automatically generate a `bru_mapper` for `rspde_metric_graph`
#'   models.
#' @param model Model class (contains a metric graph object)
#' @param \dots arguments passed to sub-methods
#' @rdname bru_get_mapper_rspde_metric_graph
bru_get_mapper_rspde_metric_graph <- function(model, ...) {
  if ((model[["f"]]$n) %% (fm_dof(model[["mesh"]])) != 0) {
    stop(paste0(
      "Incompatible degrees of freedom. SPDE: ",
      model[["f"]]$n, " and mesh: ", fm_dof(model["mesh"])
    ))
  }
  bru_mapper_MG(
    model[["mesh"]],
    n_rep = (model[["f"]]$n) / (fm_dof(model[["mesh"]]))
  )
}


#' @title bru_mapper for the fm_MG and metric_graph classes
#' @description Constructs a basic `bru_mapper` mapper object for metric
#'   graphs with internal mesh, for simple FEM model mapping.
#' @param mesh A `metric_graph` or `fm_MG` object.
#' @param n_rep number of components in linear predictor
#' @param \dots arguments passed to sub-methods
#' @returns A `bru_mapper_fmesher` or `bru_mapper_repeat` object
#' @rdname bru_mapper_MG
bru_mapper_MG <- function(mesh, n_rep = 1, ...) {
  UseMethod("bru_mapper_MG")
}

#' @rdname bru_mapper_MG
#' @export
bru_mapper_MG.fm_MG <- function(mesh, n_rep = 1, ...) {
  mapper <- inlabru::bru_mapper_fmesher(mesh)
  if (n_rep > 1) {
    mapper <- inlabru::bru_mapper_repeat(mapper, n_rep = n_rep)
  }
  mapper
}

#' @rdname bru_mapper_MG
#' @export
bru_mapper_MG.metric_graph <- function(mesh, n_rep = 1, ...) {
  bru_mapper_MG(fm_as_MG(mesh), n_rep = n_rep, ...)
}


# fmesher functions ----

#' @title fmesher interface methods for `metric_graph` objects
#' @description Interface for `metric_graph` objects
#' @name fm_MG
NULL

#' @describeIn fm_MG Wrap a `metric_graph` R6 object in an S3 container object
#'   with class `fm_MG` and subclass `fm_MGG` or `fm_MGM`. The graph is
#'   accessible via `x$graph`.
#' @param MGG indicator for `fm_MGG` (`TRUE`) or `fm_MGM` (`FALSE`), or NULL
#' (default) to auto-determine based on the object contents.
#' @export
fm_as_MG <- function(x, MGG = NULL) {
  UseMethod("fm_as_MG")
}

#' @rdname fm_MG
#' @export
fm_as_MG.fm_MG <- function(x, MGG = NULL) {
  if (is.null(MGG)) {
    if (inherits(x, c("fm_MGG", "fm_MGM"))) {
      return(x)
    }
    MGG <- is.null(x[["graph"]][["mesh"]])
  }

  if ((MGG && inherits(x, "fm_MGG")) ||
    (!MGG && inherits(x, "fm_MGM"))) {
    return(x)
  }
  cl <- setdiff(class(x), c("fm_MGG", "fm_MGM"))
  if (isTRUE(MGG)) {
    class(x) <- c("fm_MGG", cl)
  } else {
    stopifnot(!is.null(x[["mesh"]]))
    class(x) <- c("fm_MGM", cl)
  }
  x
}

#' @rdname fm_MG
#' @export
fm_as_MG.metric_graph <- function(x, MGG = NULL) {
  if (is.null(MGG)) {
    MGG <- is.null(x[["mesh"]])
  }
  if (isTRUE(MGG)) {
    class_names <- c("fm_MGG", "fm_MG")
  } else {
    stopifnot(!is.null(x[["mesh"]]))
    class_names <- c("fm_MGM", "fm_MG")
  }
  structure(list(graph = x), class = class_names)
}

#' @describeIn fm_as_fm Wrap a `metric_graph` object with class `fm_MGG` or
#' `fm_MGM`, with [fm_as_MG()].
#' @param MGG indicator for `fm_MGG` (TRUE) or `fm_MGM` (FALSE), or NULL
#' (default) to auto-determine based on the object contents.
#' Passed on to [fm_as_MG()].
#' @export
fm_as_fm.fm_MG <- function(x, ..., MGG = NULL) {
  fm_as_MG(x, MGG = MGG)
}

#' @describeIn fm_as_fm Wrap a `metric_graph` object with class `fm_MGG` or
#' `fm_MGM`, with [fm_as_MG()].
#' @export
fm_as_fm.metric_graph <- function(x, ..., MGG = NULL) {
  fm_as_MG(x, MGG = MGG)
}

#' @describeIn fm_MG Extract the graph R6 object from an `fm_MG` object.
#' If the input is already a `metric_graph` it is returned unchanged.
#' @export
fm_MG_graph <- function(x) {
  UseMethod("fm_MG_graph")
}

#' @rdname fm_MG
#' @export
fm_MG_graph.fm_MG <- function(x) {
  x[["graph"]]
}

#' @rdname fm_MG
#' @export
fm_MG_graph.metric_graph <- function(x) {
  x
}


#' @describeIn fm_MG Construct an interpolation/basis matrix
#' @param x metric_graph or fm_MG object
#' @param loc Observation locations, can be either MGG coordinates, MGM
#'   coordinates or Euclidean coordinates (passed to [fm_bary()])
#' @param weights Optional weight vector, one weight for each location
#' @inheritParams fm_basis
#' @export
#' @returns A `fm_basis` object; a list of evaluator information objects,
#' at least a matrix `A` and logical vector `ok`.
fm_basis.fm_MG <- function(x,
                           loc,
                           weights = NULL,
                           ...,
                           full = FALSE) {
  # Ensure we have fm_MGG or fm_MGM. The code after works for both.
  x <- fm_as_MG(x)
  if (is.null(weights)) {
    weights <- rep(1.0, NROW(loc))
  } else if (length(weights) == 1L) {
    weights <- rep(weights, NROW(loc))
  }

  info <- list()
  # use metric graph function to get basis functions
  # obtain bary wrt to MGG or MGM, as appropriate
  barys <- fm_bary(x, loc)
  simplex <- fm_bary_simplex(x, barys)
  n <- NROW(barys)
  info$A <- Matrix::sparseMatrix(
    i = c(seq_len(n), seq_len(n)),
    j = as.vector(simplex),
    x = c(weights * barys$where[, 1], weights * barys$where[, 2]),
    dims = c(n, fm_dof(x))
  )
  info[["ok"]] <- !is.na(barys$index)

  fm_basis(
    structure(
      info,
      class = "fm_basis"
    ),
    full = full
  )
}

#' @rdname fm_MG
#' @export
fm_basis.metric_graph <- function(x,
                                  loc,
                                  weights = NULL,
                                  ...,
                                  full = FALSE) {
  x <- fm_as_MG(x)
  fm_basis(x, loc, weights = weights, ..., full = full)
}


#' @describeIn fm_MG Compute an `fm_MGG_bary` or `fm_MGM_bary` object
#' @inheritParams fm_bary
#' @param MGG indicator for the barycentric coordinates related to the graph
#'   (MGG) or mesh (MGM), or NULL. Passed on to [fm_as_MG()]
#' @export
fm_bary.metric_graph <- function(mesh,
                                 loc,
                                 ...) {
  fm_bary(fm_as_MG(mesh), loc, ...)
}

#' @describeIn fm_MG Compute an `fm_MGG_bary` object
#' @export
fm_bary.fm_MGG <- function(mesh,
                           loc,
                           ...) {
  mesh <- fm_as_MG(mesh, MGG = TRUE)
  if (inherits(loc, "fm_bary")) {
    bary <- fm_as_MGG_bary(loc, graph = mesh)
  } else if (inherits(loc, "sfg") || inherits(loc, "sf") ||
    inherits(loc, "sfc")) {
    # check the crs of point and convert to the same crs as graph (or
    # coordinates handles this)
    if(!is.null(fm_MG_graph(mesh)$.__enclose_env__$private$crs)){
      loc <- sf::st_transform(loc, fm_MG_graph(mesh)$.__enclose_env__$private$crs)
    }
    res <- fm_MG_graph(mesh)$coordinates(XY = sf::st_coordinates(loc))
    bary <- fm_as_MGG_bary(loc = res)

  } else { # assume loc is a matrix with XY-coordinates in the same coordinate
           # system as the mesh
    res <- fm_MG_graph(mesh)$coordinates(XY = loc)
    bary <- fm_as_MGG_bary(loc = res)
  }
  bary
}

#' @describeIn fm_MG Compute an `fm_MGM_bary` object
#' @export
fm_bary.fm_MGM <- function(mesh,
                           loc,
                           ...) {
  mesh <- fm_as_MG(mesh, MGG = FALSE)
  if (inherits(loc, "fm_bary")) {
    bary_coord <- fm_as_MGM_bary(loc, graph = mesh)
  } else if (inherits(loc, "sfg") || inherits(loc, "sf") ||
    inherits(loc, "sfc")) {
    if(!is.null(fm_MG_graph(mesh)$.__enclose_env__$private$crs)){
      loc <- sf::st_transform(loc, fm_MG_graph(mesh)$.__enclose_env__$private$crs)
    }
    res <- fm_MG_graph(mesh)$coordinates(XY = sf::st_coordinates(loc))
    bary <- fm_as_MGG_bary(loc = res)
    bary_coord <- MGG_to_MGM(bary, mesh)
  } else {
    # warning(paste0("fm_bary: loc was interpreted as a matrix of Euclidean ",
    #                "coordinates in the same reference system as mesh."))
    res <- fm_MG_graph(mesh)$coordinates(XY = loc)
    bary <- fm_as_MGG_bary(loc = res)
    bary_coord <- MGG_to_MGM(bary, mesh)
  }
  bary_coord
}

#' @describeIn fm_MG Returns a `metric_graph` graph or mesh edge index matrix
#' @export
fm_bary_simplex.metric_graph <- function(mesh, bary = NULL, ...) {
  mesh <- fm_as_MG(mesh)
  simplex <- fm_bary_simplex(mesh)
  return(simplex)
}

#' @describeIn fm_MG Returns a `metric_graph` graph edge index matrix
#' @export
fm_bary_simplex.fm_MGG <- function(mesh, bary = NULL, ...) {
  simplex <- fm_MG_graph(mesh)[["E"]]
  if (!is.null(bary)) {
    bary <- fm_as_MGG_bary(bary, graph = mesh)
    simplex <- simplex[bary$index, , drop = FALSE]
  }
  return(simplex)
}

#' @describeIn fm_MG Returns a `metric_graph` mesh edge index matrix
#' @export
fm_bary_simplex.fm_MGM <- function(mesh, bary = NULL, ...) {
  simplex <- fm_MG_graph(mesh)[["mesh"]][["E"]]
  if (!is.null(bary)) {
    bary <- fm_as_MGM_bary(bary, graph = mesh)
    simplex <- simplex[bary$index, , drop = FALSE]
  }
  return(simplex)
}



#' @describeIn fm_MG Return manifold type, always "G1" for `metric_graph`
#' @export
fm_manifold_get.metric_graph <- function() {
  "G1"
}


#' @describeIn fm_MG Returns the degrees of freedom (number of vertices in the
#' graph or mesh, depending on the object). `fm_as_MG(x)` is used to ensure
#' `fm_MGG` or `fm_MGM` class.
#' @export
fm_dof.metric_graph <- function(x) {
  fm_dof(fm_as_MG(x))
}

#' @describeIn fm_MG `fm_MGG` method for `fm_dof()`
#' @export
fm_dof.fm_MGG <- function(x) {
  fm_MG_graph(x)[["nV"]]
}

#' @describeIn fm_MG `fm_MGM` method for `fm_dof()`
#' @export
fm_dof.fm_MGM <- function(x) {
  NROW(fm_MG_graph(x)[["mesh"]][["VtE"]])
}


#' @export
#' @describeIn fm_MG `metric_graph` integration. Supported samplers:
#'   * `NULL` for integration over the entire domain;
#'   * A tibble with a named column containing a matrix with single edge
#'     intervals (ordered), and optionally a `weight` column.
#' @examples
#' if (requireNamespace("MetricGraph")) {
#'   edge1 <- rbind(c(0, 0), c(1, 0))
#'   edge2 <- rbind(c(0, 0), c(0, 1))
#'   edge3 <- rbind(c(0, 1), c(-1, 1))
#'   theta <- seq(from = pi, to = 3 * pi / 2, length.out = 20)
#'   edge4 <- cbind(sin(theta), 1 + cos(theta))
#'   edges <- list(edge1, edge2, edge3, edge4)
#'   graph <- MetricGraph::metric_graph$new(edges = edges)
#'   graph$build_mesh(h = 0.01)
#'   p1 <- simple_path_MGG(
#'     graph = graph,
#'     start = cbind(1, 0.2),
#'     edges = c(2),
#'     end = cbind(3, 0.8)
#'   )
#'   samplers <- tibble::tibble(x = list(p1), weight = c(1))
#'   ips <- fm_int(
#'     graph,
#'     samplers
#'   )
#' }
#' @inheritParams fm_int
fm_int.metric_graph <- function(domain,
                                samplers = NULL,
                                name = "x",
                                ...) {
  fm_int(
    domain = fm_as_MG(domain),
    samplers = samplers,
    name = name,
    ...
  )
}


#' @describeIn fm_MG Integrate on an `fm_MGM` object
#' @export
fm_int.fm_MGM <- function(domain,
                          samplers = NULL,
                          name = "x",
                          ...) {
  ips <- fm_int(
    domain = fm_as_MG(domain, MGG = TRUE),
    samplers = samplers,
    name = name,
    ...
  )
  ips[[name]] <- fm_as_MGM_bary(ips[[name]], graph = domain)
  return(ips)
}

#' @describeIn fm_MG Integration on `metric_graph`; requires a mesh in the
#'   graph.
#' @export
fm_int.fm_MGG <- function(domain,
                          samplers = NULL,
                          name = "x",
                          ...) {
  ips <- list()
  if (is.null(fm_MG_graph(domain)[["mesh"]])) {
    stop("There is no mesh")
  }

  # All mesh locations in MGG format
  mesh_MGG <- fm_as_MGG_bary(fm_MG_graph(domain)$mesh$VtE)
  # if(is.null(samplers)){
  #   samplers <- tibble::tibble(
  #     x = cbind(domain[,1], domain[,2]),
  #     weight = 1,
  #     .block = 1L
  #   )
  # }
  if (is.data.frame(samplers)) {
    .block <- seq_len(NROW(samplers))
  }

  for (j in seq_len(NROW(samplers))) {
    # for a graph interval
    subsampler <- samplers[[name]][[j]]
    theweight <- samplers[["weight"]][[j]]
    ips_edge <- list()
    for (k in seq_len(nrow(subsampler))) {
      interedge <- subsampler[k, , drop = TRUE]
      if (!inherits(interedge, "fm_MGG_interval")) {
        interedge <- fm_as_MGG_intervals(
          graph = domain,
          start = fm_as_MGG_bary(list(
            interedge$start$index,
            interedge$start$where[1, 2]
          )),
          end = fm_as_MGG_bary(list(
            interedge$end$index,
            interedge$end$where[1, 2]
          ))
        )
      }
      the.block <- .block[j]
      # Simpson's rule integration
      # mesh vertices on edge of interest (+ end points)
      # make sure we have only inter edges:
      if (interedge$start$index != interedge$end$index) {
        stop("samplers contain intervals which are not interedges.")
      }
      loc_trap <- sort(unique(c(
        0,
        mesh_MGG$where[mesh_MGG$index == interedge$start$index, 2],
        1
      )))
      loc_mid <- (loc_trap[-1] + loc_trap[-length(loc_trap)]) / 2
      # Detect mid-points inside the interval
      if (interedge$start$where[1, 2] > interedge$end$where[1, 2]) {
        inside <- (loc_mid <= interedge$start$where[1, 2]) &
          (loc_mid >= interedge$end$where[1, 2])
      } else {
        inside <- (loc_mid >= interedge$start$where[1, 2]) &
          (loc_mid <= interedge$end$where[1, 2])
      }
      # convert to MGM (call outside of for-loop and only get the desired rows
      # in this step)
      loc_mid_MGM <- MGG_to_MGM(
        coord = fm_as_MGG_bary(cbind(interedge$start$index, loc_mid)),
        graph = domain
      )
      # get the edge lengths for each mesh
      weight_mid <- fm_MG_graph(domain)$mesh$h_e[loc_mid_MGM$index]
      weight_mid[!inside] <- 0.0

      weight_trap <- c(weight_mid / 2, 0) + c(0, weight_mid / 2)
      loc_simpson <- c(loc_trap, loc_mid)
      weight_simpson <- c(weight_trap / 3, weight_mid * 2 / 3)

      m_ips <- sum(weight_simpson > 0)
      if (m_ips == 0) {
        ips_edge[[k]] <- tibble::tibble(
          x = fm_as_MGG_bary(tibble::tibble(
            index = integer(0),
            where = numeric(0)
          )),
          weight = numeric(0),
          .block = integer(0)
        )
      } else {
        ips_edge[[k]] <- tibble::tibble(
          x = fm_as_MGG_bary(loc = cbind(
            rep(interedge$start$index, m_ips),
            loc_simpson[(weight_simpson > 0)]
          )),
          weight = weight_simpson[(weight_simpson > 0)] * theweight,
          .block = the.block
        )
      }

      colnames(ips_edge[[k]])[1] <- name
    }
    ips_edge <- do.call(dplyr::bind_rows, ips_edge)
    ips[[j]] <- ips_edge
  }
  ips <- do.call(dplyr::bind_rows, ips)
  if (NROW(ips) == 0) {
    ips <- tibble::tibble(
      x = fm_as_MGG_bary(tibble::tibble(
        index = integer(0),
        where = numeric(0)
      )),
      weight = numeric(0),
      .block = integer(0)
    )
    colnames(ips)[1] <- name
  }
  ips
}




# MetricGraph specific functions----
#' @title Make a (`mesh`, `fm_bary`) object from MGG coordinates
#' @description
#' Create a (`mesh`, `fm_bary`) object from MGG coordinates.
#'
#' @param graph metric_graph that the location should be mapped to.
#' @param coord MGG coordinates
#' @author Karina Lilleborge \email{karina.lilleborge@@gmail.com}
#' @returns A (`mesh`, `fm_bary`) object
#' @export
#' @family object creation and conversion
#' @examples
#' if (requireNamespace("MetricGraph", quietly = TRUE)) {
#'   edge1 <- rbind(c(0, 0), c(1, 0))
#'   edge2 <- rbind(c(0, 0), c(0, 1))
#'   edge3 <- rbind(c(0, 1), c(-1, 1))
#'   theta <- seq(from = pi, to = 3 * pi / 2, length.out = 20)
#'   edge4 <- cbind(sin(theta), 1 + cos(theta))
#'   edges <- list(edge1, edge2, edge3, edge4)
#'   graph <- MetricGraph::metric_graph$new(edges = edges)
#'   graph$build_mesh(h = 0.005)
#'   mgm <- MGG_to_MGM(
#'     fm_as_MGG_bary(cbind(1, 0.5)),
#'     graph
#'   )
#'   mgm
#' }
#'
#' @keywords internal
MGG_to_MGM <- function(coord, graph) {
  stopifnot(inherits(coord, "fm_MGG_bary"))
  if (is.null(fm_MG_graph(graph)[["mesh"]])) {
    stop("There is no mesh")
  }

  res <- fm_MG_graph(graph)$.__enclos_env__$private$PtE_to_mesh(
    cbind(coord$index, coord$where[, 2])
  )
  res <- fm_as_MGM_bary(tibble::tibble(
    index = res[, 1],
    where = cbind(
      1 - as.numeric(res[, 2]),
      as.numeric(res[, 2])
    )
  ))
  return(res)

  mesh_MGG <- fm_as_MGG_bary(fm_MG_graph(graph)$mesh$VtE)
  mesh_edge_len <- fm_MG_graph(graph)$mesh$h_e
  # storage for new coordinates
  new_coord <- fm_as_MGM_bary(tibble::tibble(
    index = integer(NROW(coord)),
    where = numeric(NROW(coord))
  ))
  for (i in seq_len(NROW(coord))) {
    # which mesh vertices are on the given edge
    ids <- (mesh_MGG$index == as.integer(coord$index[i]))
    # MGG coordinates for those vertices
    edge_MGG <- mesh_MGG[ids, ]
    if (sum(ids) == 0) {
      # we have no mesh locations here, so we need to determine the vertices on
      # either side
      vertices_MGG <- fm_MG_graph(graph)$E[coord$index[i], ]

      index_MGM <- which.max(
        (fm_MG_graph(graph)$mesh$E[, 1] == vertices_MGG[1]) &
          (fm_MG_graph(graph)$mesh$E[, 2] == vertices_MGG[2])
      )
      where_MGM <- coord$where[i, 2]
    } else if (sum(ids) == 1) {
      # there is only one mesh vertex on the edge
      index_on_edge <- which(ids)
      # coord[i, ] is before or past the mesh node
      further <- ((as.numeric(coord$where[i, 2]) - edge_MGG$where[2]) >= 0)
      # find the vertex on the other side of coord[i, ]
      graph_vertex <- fm_MG_graph(graph)$E[edge_MGG$index, further * 1 + 1]
      if (further) {
        # find the edge index that connects (mesh_vertex, end_vertex)
        index_MGM <- which.max(
          (fm_MG_graph(graph)$mesh$E[, 1] == index_on_edge) &
            (fm_MG_graph(graph)$mesh$E[, 2] == graph_vertex)
        )
        mesh_h_e <- mesh_edge_len[index_MGM]
        where_MGM <- as.numeric(
          (as.numeric(coord$where[i, 2]) - edge_MGG$where[, 2]) / mesh_h_e
        )
      } else {
        # find the edge index that connects (start_vertex, mesh_vertex)
        index_MGM <- which.max(
          (fm_MG_graph(graph)$mesh$E[, 1] == graph_vertex) &
            (fm_MG_graph(graph)$mesh$E[, 2] == index_on_edge)
        )
        mesh_h_e <- mesh_edge_len[index_MGM]
        where_MGM <- 1 - (as.numeric((edge_MGG$where[, 2] -
          as.numeric(coord$where[i, 2])) /
          mesh_h_e))
      }
    } else {
      # order the mesh_MGG locations:
      ordering <- order(edge_MGG$where[, 2])
      edge_MGG_o <- edge_MGG[ordering, ]
      # find the mesh point index where coord[i,] is next to
      index_on_edge <- which.max((edge_MGG_o$where[, 2] -
        as.numeric(coord$where[i, 2])) >= 0)
      # coord[i, ] is between these two mesh locs
      mesh_indices <-
        which(ids)[(ordering[c(index_on_edge - 1, index_on_edge)])]
      index_MGM <- which.max(
        (fm_MG_graph(graph)$mesh$E[, 1] == mesh_indices[1]) &
          (fm_MG_graph(graph)$mesh$E[, 2] == mesh_indices[2])
      )
      mesh_h_e <- mesh_edge_len[index_MGM]
      where_MGM <- 1 - as.numeric((edge_MGG_o$where[index_on_edge, 2] -
        as.numeric(coord$where[i, 2])) /
        mesh_h_e) # normalized
    }
    if (length(c(index_MGM, where_MGM)) != 2) {
      stop(paste0(
        "We found ", sum(ids),
        " mesh locations and we got length(index_MGM)=",
        length(index_MGM), " and length(where_MGM)=",
        length(where_MGM), ". Coordinate: (",
        coord$index, ",", coord$where[, 2], ")"
      ))
    }
    new_coord[i, ] <- tibble::tibble(
      index = index_MGM,
      where = cbind(1 - where_MGM, where_MGM)
    )
  }
  new_coord
}

#' @title Convert (`mesh`, `fm_bary`) coordinates to a (`graph`, `fm_bary`)
#'   coordinates
#' @description
#' Create a (`graph`, `fm_bary`) object from MGM coordinates.
#'
#' @param graph metric_graph that the location should be mapped to.
#' @param coord MGM coordinates
#' @author Karina Lilleborge \email{karina.lilleborge@@gmail.com}
#' @returns An (`graph`, `fm_bary`) object
#' @export
#' @family object creation and conversion
#' @examples
#' if (requireNamespace("MetricGraph", quietly = TRUE)) {
#'   edge1 <- rbind(c(0, 0), c(1, 0))
#'   edge2 <- rbind(c(0, 0), c(0, 1))
#'   edge3 <- rbind(c(0, 1), c(-1, 1))
#'   theta <- seq(from = pi, to = 3 * pi / 2, length.out = 20)
#'   edge4 <- cbind(sin(theta), 1 + cos(theta))
#'   edges <- list(edge1, edge2, edge3, edge4)
#'   graph <- MetricGraph::metric_graph$new(edges = edges)
#'   graph$build_mesh(h = 0.01)
#'   mgg <- MGM_to_MGG(
#'     fm_as_MGM_bary(cbind(5, 1)),
#'     graph
#'   )
#'   mgg
#' }
#'
#' @keywords internal
MGM_to_MGG <- function(coord, graph) {
  stopifnot(inherits(coord, "fm_MGM_bary"))
  if (is.null(fm_MG_graph(graph)[["mesh"]])) {
    stop("There is no mesh")
  }
  mesh_MGG <- fm_as_MGG_bary(fm_MG_graph(graph)$mesh$VtE)
  eps <- min(fm_MG_graph(graph)$mesh$h_e)
  new_coord <- fm_as_MGG_bary(tibble::tibble(
    index = integer(NROW(coord)),
    where = numeric(NROW(coord))
  ))
  for (i in seq_len(NROW(coord))) {
    # mesh vertex indices for neighboring loc
    mesh_vertices <- fm_MG_graph(graph)$mesh$E[coord$index[i], ]
    # which graph coordinates for the neighboring vertices
    graph_edge_r <- mesh_MGG[mesh_vertices[2], ]
    graph_edge_l <- mesh_MGG[mesh_vertices[1], ]
    # as long as they are on the same edge
    if (graph_edge_l$index == graph_edge_r$index) {
      # neighbouring mesh nodes are on the same edge
      new_coord[i, ] <- fm_as_MGG_bary(
        tibble::tibble(
          index = graph_edge_l$index,
          where = (1 - coord$where[i, 2]) * graph_edge_l$where[1, 2] +
            coord$where[i, 2] * graph_edge_r$where[1, 2]
        )
      )
    } else {
      # neighbouring mesh nodes are not on same edge (one is a vertex)
      on_vertex_r <- (abs(graph_edge_r$where[1, 2] - c(0, 1)) < eps)
      on_vertex_l <- (abs(graph_edge_l$where[1, 2] - c(0, 1)) < eps)
      if (any(on_vertex_r)) {
        new_coord[i, ] <- fm_as_MGG_bary(
          tibble::tibble(
            index = graph_edge_l$index,
            where = (1 - coord$where[i, 2]) * graph_edge_l$where[1, 2] +
              coord$where[i, 2] * c(0, 1)[on_vertex_r]
          )
        )
      } else {
        new_coord[i, ] <- fm_as_MGG_bary(
          tibble::tibble(
            index = graph_edge_r$index,
            where = (1 - coord$where[i, 2]) * c(0, 1)[on_vertex_l] +
              coord$where[i, 2] * graph_edge_r$where[1, 2]
          )
        )
      }
    }
  }
  fm_as_MGG_bary(new_coord)
}


#' @title Make a (`mesh`, `fm_bary`) object
#' @description
#' Create a (`mesh`, `fm_bary`) object
#'
#' @param loc MGM coordinates
#' @param graph metric_graph that the location should be mapped to (must be
#'   provided if loc should be converted)
#' @author Karina Lilleborge \email{karina.lilleborge@@gmail.com}
#' @returns An (`graph`, `fm_bary`) object from `matrix`, `data.frame`, `list`,
#'   `tibble` etc.
#' @export
#' @family object creation and conversion
#' @examples
#' if (requireNamespace("MetricGraph", quietly = TRUE)) {
#'   edge1 <- rbind(c(0, 0), c(1, 0))
#'   edge2 <- rbind(c(0, 0), c(0, 1))
#'   edge3 <- rbind(c(0, 1), c(-1, 1))
#'   theta <- seq(from = pi, to = 3 * pi / 2, length.out = 20)
#'   edge4 <- cbind(sin(theta), 1 + cos(theta))
#'   edges <- list(edge1, edge2, edge3, edge4)
#'   graph <- MetricGraph::metric_graph$new(edges = edges)
#'   graph$build_mesh(h = 0.01)
#'   m <- fm_as_MGM_bary(
#'     cbind(1, 1),
#'     graph
#'   )
#'   class(m) # "fm_MGM_bary", "fm_bary", "tbl_df", "tbl", "data.frame"
#' }
#'
#' @export
fm_as_MGM_bary <- function(loc, graph = NULL) {
  if (inherits(loc, "fm_MGM_bary")) {
    return(loc)
  }
  if (inherits(loc, "fm_MGG_bary")) {
    if (is.null(graph)) {
      stop("Graph must be provided to convert from MGG to MGM.")
    }
    return(MGG_to_MGM(coord = loc, graph = graph))
  }
  if (is.matrix(loc)) {
    res <- tibble::tibble(
      index = as.integer(loc[, 1]),
      where = cbind(
        1 - as.numeric(loc[, 2]),
        as.numeric(loc[, 2])
      )
    )
  } else if (!tibble::is_tibble(loc)) {
    res <- tibble::tibble(
      index = as.integer(loc[[1]]),
      where = cbind(
        1 - as.numeric(loc[[2]]),
        as.numeric(loc[[2]])
      )
    )
  } else {
    stopifnot(tibble::is_tibble(loc))
    if (is.matrix(loc$where)) {
      res <- loc
    } else {
      res <- tibble::tibble(
        index = as.integer(loc$index),
        where = cbind(
          1 - as.numeric(loc$where),
          as.numeric(loc$where)
        )
      )
    }
  }
  coord <-
    structure(
      res,
      class = c("fm_MGM_bary", "fm_bary", "tbl_df", "tbl", "data.frame")
    )
}


#' @title Make a (`graph`, `fm_bary`) object
#' @description
#' Create a (`graph`, `fm_bary`) object from `matrix`, `data.frame`, `list`,
#' `tibble` etc.
#'
#' @param loc MGG coordinates
#' @param graph metric_graph that the location should be mapped to (must be
#'   provided if loc should be converted)
#' @author Karina Lilleborge \email{karina.lilleborge@@gmail.com}
#' @returns A (`mesh`, `fm_bary`) object
#' @export
#' @family object creation and conversion
#' @examples
#' if (requireNamespace("MetricGraph", quietly = TRUE)) {
#'   edge1 <- rbind(c(0, 0), c(1, 0))
#'   edge2 <- rbind(c(0, 0), c(0, 1))
#'   edge3 <- rbind(c(0, 1), c(-1, 1))
#'   theta <- seq(from = pi, to = 3 * pi / 2, length.out = 20)
#'   edge4 <- cbind(sin(theta), 1 + cos(theta))
#'   edges <- list(edge1, edge2, edge3, edge4)
#'   graph <- MetricGraph::metric_graph$new(edges = edges)
#'   m <- fm_as_MGG_bary(cbind(1, 0.5))
#'   class(m) # "fm_MGG_bary", "fm_bary", "tbl_df", "tbl", "data.frame"
#' }
#'
#' @export
fm_as_MGG_bary <- function(loc, graph = NULL) {
  if (inherits(loc, "fm_MGG_bary")) {
    return(loc)
  }
  if (inherits(loc, "fm_MGM_bary")) {
    if (is.null(graph)) {
      stop("Graph must be provded to convert from MGM to MGG.")
    }
    return(MGM_to_MGG(coord = loc, graph = graph))
  }
  # TO DO: add check here for sf (and then convert to correct crs)

  if (is.matrix(loc)) {
    res <- tibble::tibble(
      index = as.integer(loc[, 1]),
      where = cbind(
        1 - as.numeric(loc[, 2]),
        as.numeric(loc[, 2])
      )
    )
  } else if (!tibble::is_tibble(loc)) {
    res <- tibble::tibble(
      index = as.integer(loc[[1]]),
      where = cbind(
        1 - as.numeric(loc[[2]]),
        as.numeric(loc[[2]])
      )
    )
  } else {
    stopifnot(tibble::is_tibble(loc))
    if (is.matrix(loc$where)) {
      res <- loc
    } else {
      res <- tibble::tibble(
        index = as.integer(loc$index),
        where = cbind(
          1 - as.numeric(loc$where),
          as.numeric(loc$where)
        )
      )
    }
  }
  coord <-
    structure(
      res,
      class = c("fm_MGG_bary", "fm_bary", "tbl_df", "tbl", "data.frame")
    )
}


#' @title Make inter edge intervals on graph object
#' @description
#' Create an `fm_MGG_interval` object.
#'
#' @param start Start locations for inter-graph-edge intervals.
#' @param end End locations for inter-graph-edge intervals.
#' @param graph `metric_graph` that the intervals should be mapped to. Must be
#'   provided if input should be converted
#' @author Karina Lilleborge \email{karina.lilleborge@@gmail.com}
#' @returns An `fm_MGG_interval` object
#' @export
#' @family object creation and conversion
#' @examples
#' if (requireNamespace("MetricGraph", quietly = TRUE)) {
#'   edge1 <- rbind(c(0, 0), c(1, 0))
#'   edge2 <- rbind(c(0, 0), c(0, 1))
#'   edge3 <- rbind(c(0, 1), c(-1, 1))
#'   theta <- seq(from = pi, to = 3 * pi / 2, length.out = 20)
#'   edge4 <- cbind(sin(theta), 1 + cos(theta))
#'   edges <- list(edge1, edge2, edge3, edge4)
#'   graph <- MetricGraph::metric_graph$new(edges = edges)
#'   int <- fm_as_MGG_intervals(
#'     cbind(1, 0.8),
#'     cbind(1, 0.5),
#'     graph
#'   )
#'   int
#' }
#'
#' @export
fm_as_MGG_intervals <- function(start, end, graph = NULL) {
  start <- fm_as_MGG_bary(start, graph = graph)
  end <- fm_as_MGG_bary(end, graph = graph)
  if (any(start$index != end$index)) {
    stop("Not all start- and end points are inter edge intervals")
  }
  inter_edge_interval <- structure(
    tibble::tibble(
      start = start,
      end = end
    ),
    class = c("fm_MGG_intervals", "tbl_df", "tbl", "data.frame")
  )
}



#' @title Make an interval on graph object from path
#' @description
#' Create a `fm_MGG_interval` object from specified path (either geometric line
#' or a simple path)
#'
#' @param graph metric_graph that the interval should be mapped to.
#' @param path `sf::st_geometry` (`LINESTRING`) which is a subset of the graph
#'   or a list with three elements; start MGG coordinate, an ordered list of
#'   edges and end MGG coordinate
#' @author Karina Lilleborge \email{karina.lilleborge@@gmail.com}
#' @returns A `fm_MGG_interval` object
#' @export
#' @family object creation and conversion
#' @examples
#' if (requireNamespace("MetricGraph")) {
#'   edge1 <- rbind(c(0, 0), c(1, 0))
#'   edge2 <- rbind(c(0, 0), c(0, 1))
#'   edge3 <- rbind(c(0, 1), c(-1, 1))
#'   theta <- seq(from = pi, to = 3 * pi / 2, length.out = 20)
#'   edge4 <- cbind(sin(theta), 1 + cos(theta))
#'   edges <- list(edge1, edge2, edge3, edge4)
#'   graph <- MetricGraph::metric_graph$new(edges = edges)
#'   path <- list(
#'     c(1, 0.5),
#'     c(2,3),
#'     c(4, 0.2)
#'   )
#'   intervals <- fm_MGG_intervals(
#'     graph,
#'     path
#'   )
#'   path
#' }
#'
#' @export
fm_MGG_intervals <- function(graph, path){
  if(is.list(path)){
    if(length(path)==3){
      return(simple_path_MGG(graph,
                             path[[1]],
                             path[[2]],
                             path[[3]]))
    } else{
      stop("The list provided as path is not of the correct length 3.")
    }
  } else if(inherits(path, "sfc_LINESTRING")){
    return(geom_path_to_path_MGG(path, graph))
  } else{
    stop("The path provided is not of any of the supported object types.")
  }
}

#' @title Make an interval on graph object from simple path
#' @description
#' Create a `fm_MGG_interval` object from known start (MGG), end (MGG) and
#' visiting edges (MGG).
#'
#' @param graph metric_graph that the interval should be mapped to.
#' @param start `fm_MGG_bary` or [fm_as_MGG_bary()] compatible coordinates for
#'   start
#' @param edges Ordered list of edge indices related to MGG
#' @param end `fm_MGG_bary` or [fm_as_MGG_bary()] compatible coordinates for end
#' @author Karina Lilleborge \email{karina.lilleborge@@gmail.com}
#' @returns A `fm_MGG_interval` object
#' @export
#' @family object creation and conversion
#' @examples
#' if (requireNamespace("MetricGraph")) {
#'   edge1 <- rbind(c(0, 0), c(1, 0))
#'   edge2 <- rbind(c(0, 0), c(0, 1))
#'   edge3 <- rbind(c(0, 1), c(-1, 1))
#'   theta <- seq(from = pi, to = 3 * pi / 2, length.out = 20)
#'   edge4 <- cbind(sin(theta), 1 + cos(theta))
#'   edges <- list(edge1, edge2, edge3, edge4)
#'   graph <- MetricGraph::metric_graph$new(edges = edges)
#'   path <- simple_path_MGG(
#'     graph,
#'     start = cbind(1, 0.5),
#'     edges = c(2),
#'     end = cbind(3, 0.6)
#'   )
#'   path
#' }
#'
#' @keywords internal
simple_path_MGG <- function(graph,
                            start,
                            edges,
                            end) {
  # check if the graph does have circles
  start_MGG <- fm_as_MGG_bary(start, graph = graph)
  end_MGG <- fm_as_MGG_bary(end, graph = graph)
  if (length(edges) > 0) {
    # check direction from start to edges[1]
    v1 <- fm_MG_graph(graph)$E[as.integer(start_MGG$index), ]
    v2 <- fm_MG_graph(graph)$E[as.integer(edges[1L]), ]
    if (sum(v2 %in% v1[1]) > 0) {
      # if the (e,0) vertex is in v2
      end_vertex <- 0
    }
    if (sum(v2 %in% v1[2]) > 0) {
      # if the (e,1) vertex is in v2
      end_vertex <- 1
    }
    # make storage for the inter edge intervals for each of the edge
    # (start, end)
    n <- length(edges) + 2
    inter_edge_intervals <-
      fm_as_MGG_intervals(
        fm_as_MGG_bary(tibble::tibble(
          index = integer(n),
          where = numeric(n)
        )),
        fm_as_MGG_bary(tibble::tibble(
          index = integer(n),
          where = numeric(n)
        ))
      )
    inter_edge_intervals[1, ] <- tibble::tibble(
      start = start_MGG,
      end = fm_as_MGG_bary(
        cbind(
          as.numeric(start_MGG$index),
          as.numeric(end_vertex)
        )
      )
    )
    # start and end must be determined
    # check direction for edges
    for (i in seq_len(length(edges))) {
      end_vertex <- c(0, 1)[!(v2 %in% v1[end_vertex + 1])]
      if (end_vertex == 0) start_vertex <- 1
      if (end_vertex == 1) start_vertex <- 0
      inter_edge_intervals[i + 1, ] <- tibble::tibble(
        start = fm_as_MGG_bary(cbind(
          as.integer(edges[i]),
          as.numeric(start_vertex)
        )),
        end = fm_as_MGG_bary(cbind(
          as.integer(edges[i]),
          as.numeric(end_vertex)
        ))
      )
      v1 <- fm_MG_graph(graph)$E[as.integer(edges[i]), ]
      v2 <- fm_MG_graph(graph)$E[as.integer(edges[i + 1]), ]
    }
    v2 <- fm_MG_graph(graph)$E[as.integer(end_MGG$index), ]
    # check direction
    start_vertex <- c(0:1)[(v2 %in% v1[end_vertex + 1])]
    inter_edge_intervals[length(edges) + 2, ] <- tibble::tibble(
      start = fm_as_MGG_bary(cbind(
        as.integer(end_MGG$index),
        as.numeric(start_vertex)
      )),
      end = fm_as_MGG_bary(cbind(
        as.integer(end_MGG$index),
        as.numeric(end_MGG$where[, 2L])
      ))
    )
  } else { # there are no whole edges visited (edges=c())
    if (as.integer(start_MGG$index) == as.integer(end_MGG$index)) { # same edge
      inter_edge_intervals <- fm_as_MGG_intervals(
        start = fm_as_MGG_bary(tibble::tibble(
          index = integer(1),
          where = numeric(1)
        )),
        end = fm_as_MGG_bary(tibble::tibble(
          index = integer(1),
          where = numeric(1)
        ))
      )
      inter_edge_intervals[1, ] <- tibble::tibble(
        start = fm_as_MGG_bary(cbind(
          as.integer(start_MGG$index),
          as.numeric(start_MGG$where[, 2])
        )),
        end = fm_as_MGG_bary(cbind(
          as.integer(start_MGG$index),
          as.numeric(end_MGG$where[, 2])
        ))
      )
    } else { # neighboring edges
      v1 <- fm_MG_graph(graph)$E[as.integer(start_MGG$index), ]
      v2 <- fm_MG_graph(graph)$E[as.integer(end_MGG$index), ]
      if (sum(v2 %in% v1[1]) > 0) {
        # if the (e,0) vertex is in v2
        end_vertex <- 0
      }
      if (sum(v2 %in% v1[2]) > 0) {
        # if the (e,1) vertex is in v2
        end_vertex <- 1
      }
      # make storage for the inter edge intervals for each of the edge
      # index, start and end (fm_MGG_interval)
      inter_edge_intervals <- fm_as_MGG_intervals(
        start = fm_as_MGG_bary(tibble::tibble(
          index = integer(2),
          where = numeric(2)
        )),
        end = fm_as_MGG_bary(tibble::tibble(
          index = integer(2),
          where = numeric(2)
        ))
      )
      inter_edge_intervals[1, ] <- tibble::tibble(
        start = fm_as_MGG_bary(cbind(
          as.integer(start_MGG$index),
          as.numeric(start_MGG$where[, 2])
        )),
        end = fm_as_MGG_bary(cbind(
          as.integer(start_MGG$index),
          as.numeric(end_vertex)
        ))
      )
      v2 <- fm_MG_graph(graph)$E[as.integer(end_MGG$index), ]
      # check direction
      start_vertex <- c(0:1)[(v2 %in% v1[end_vertex + 1L])]
      inter_edge_intervals[2, ] <- tibble::tibble(
        start = fm_as_MGG_bary(cbind(
          as.integer(end_MGG$index),
          as.numeric(start_vertex)
        )),
        end = fm_as_MGG_bary(cbind(
          as.integer(end_MGG$index),
          as.numeric(end_MGG$where[, 2])
        ))
      )
    }
  }
  # construct object
  path <- structure(
    inter_edge_intervals,
    class = c("fm_MGG_intervals", "tbl_df", "tbl", "data.frame")
  )
  path
}

#' @title Make an interval on graph object from sf object
#' @description
#' Create a tibble of `fm_MGG_interval` objects from `sf::st_geometry`
#' (`LINESTRING`) and a column with "id".
#'
#' @param graph metric_graph that the interval should be mapped to.
#' @param geom_path `sf::st_geometry` (`LINESTRING`) on a graph
#' @author Karina Lilleborge \email{karina.lilleborge@@gmail.com}
#' @returns A tibble containing a set of `fm_MGG_interval` objects and id
#'   referring to what `LINESTRING` it was constructed from.
#' @export
#' @family object creation and conversion
#' @examples
#' if (requireNamespace("MetricGraph") &&
#'   requireNamespace("sf")) {
#'   edge1 <- rbind(c(0, 0), c(1, 0))
#'   edge2 <- rbind(c(0, 0), c(0, 1))
#'   edge3 <- rbind(c(0, 1), c(-1, 1))
#'   theta <- seq(from = pi, to = 3 * pi / 2, length.out = 20)
#'   edge4 <- cbind(sin(theta), 1 + cos(theta))
#'   edges <- list(edge1, edge2, edge3, edge4)
#'   graph <- MetricGraph::metric_graph$new(edges = edges)
#'   graph$build_mesh(h = 0.01)
#'   geom_path <- sf::st_linestring(cbind(
#'     c(0, 0, 0.7),
#'     c(0.5, 0, 0)
#'   ))
#'   path <- geom_path_to_path_MGG(
#'     sf::st_geometry(geom_path),
#'     graph
#'   )
#'   path
#' }
#'
#' @keywords internal
geom_path_to_path_MGG <- function(geom_path, graph) {
  if (!inherits(geom_path, "sfc_LINESTRING")) {
    stop("Method not implemented. Input must be sfc_LINESTRING")
  }
  if (is.null(fm_MG_graph(graph)$mesh)) {
    stop("The graph has no mesh!")
  }
  # transform geom_path to the same coordinate reference system as the mesh (if
  # any)
  if(!is.null(fm_MG_graph(graph)$.__enclose_env__$private$crs)){
    geom_path <- sf::st_transform(geom_path, fm_MG_graph(graph)$.__enclose_env__$private$crs)
  }
  # create eps (tolerance)
  eps <- 0.01 * min(fm_MG_graph(graph)$mesh$h_e)
  # matrix with colnames X Y and L1:
  internal_XY <- sf::st_coordinates(geom_path)
  paths <- list()
  ids <- list()
  l <- 0
  for (k in unique(internal_XY[, "L1"])) {
    # a line should give us one path
    line <- internal_XY[internal_XY[, "L1"] == k, ]
    line_MGG <- fm_bary(
      fm_as_MG(graph, MGG = TRUE),
      as.matrix(line[, c("X", "Y")])
    )
    # index for number of segments added
    j <- 0
    # storing the segments (start and end separately)
    start_seg <- fm_as_MGG_bary(cbind(
      rep(NA, 2 * (NROW(line_MGG) - 1)),
      rep(NA, 2 * (NROW(line_MGG) - 1))
    ))
    end_seg <- start_seg
    for (i in seq_len(NROW(line_MGG) - 1)) {
      if (line_MGG$index[i] == line_MGG$index[i + 1]) {
        # points are on the same edge, so we can add them as they are
        # before adding, check previous point
        if (j > 0) {
          if (end_seg$index[j] == line_MGG$index[i] &&
            abs(end_seg$where[j, 2] - line_MGG$where[i, 2]) < eps) {
            # We just extend the segment
            end_seg$where[j, 2] <- line_MGG$where[i + 1, 2]
            end_seg$where[j, 1] <- 1 - end_seg$where[j, 2]
          } else {
            j <- j + 1
            start_seg[j, ] <- line_MGG[i, ]
            end_seg[j, ] <- line_MGG[i + 1, ]
          }
        } else {
          j <- j + 1
          start_seg[j, ] <- line_MGG[i, ]
          end_seg[j, ] <- line_MGG[i + 1, ]
        }
      } else {
        # need to find out if the points are giving a valid path
        start_on_vertex <- (abs(line_MGG$where[i, 2] - c(0, 1)) < eps)
        end_on_vertex <- (abs(line_MGG$where[i + 1, 2] - c(0, 1)) < eps)
        if (any(start_on_vertex)) {
          # we know that the start is a vertex -> find which index
          start_vertex <-
            fm_MG_graph(graph)$E[line_MGG$index[i], start_on_vertex]

          if (any(end_on_vertex)) {
            # the end point is also a vertex
            # the end is a vertex -> find which index
            end_vertex <-
              fm_MG_graph(graph)$E[line_MGG$index[i + 1], end_on_vertex]
            if (start_vertex != end_vertex) {
              # they are not the same vertex
              edges_start <- fm_MG_graph(graph)$E == start_vertex
              edges_end <- fm_MG_graph(graph)$E == end_vertex
              # find what edge both start vertex and end vertex are present
              edges <- rowSums(edges_start) > 0 & rowSums(edges_end) > 0
              if (sum(edges) > 1) {
                stop("Ambiguous geom_path: Multiple edge candidates.")
              }
              if (sum(edges) == 0) {
                stop(paste0(
                  "Unclear geom_path: ",
                  "Not directly connected subsequent end points."
                ))
              }
              # there is only one edge match
              edge_index <- which(edges)
              if (isTRUE(edges_start[edge_index, 1])) {
                # the start vertex is the start of the edge v = (edge_index, 0)
                j <- j + 1
                start_seg[j, ] <- fm_as_MGG_bary(cbind(edge_index, 0))
                end_seg[j, ] <- fm_as_MGG_bary(cbind(edge_index, 1))
              } else {
                # the start is the end of the edge v = (edge_index, 1)
                j <- j + 1
                start_seg[j, ] <- fm_as_MGG_bary(cbind(edge_index, 1))
                end_seg[j, ] <- fm_as_MGG_bary(cbind(edge_index, 0))
              }
            }
          } else {
            # !any(end_on_vertex) (the end is not on a vertex, but the start is)
            # we check that the start is not ambiguous wrt the end point
            edge_index <- line_MGG$index[i + 1]
            # candidates for the start
            edge_nodes <- fm_MG_graph(graph)$E[edge_index, ]
            if (!any(edge_nodes == start_vertex)) {
              stop(paste0(
                "Unclear geom_path: Not directly connected points for line",
                k, " and points ", i, " and ", i + 1, "."
              ))
            }
            # we have a valid start
            if (edge_nodes[1L] == start_vertex) {
              j <- j + 1
              start_seg[j, ] <- fm_as_MGG_bary(cbind(edge_index, 0))
              end_seg[j, ] <- line_MGG[i + 1, ]
            } else {
              j <- j + 1
              start_seg[j, ] <- fm_as_MGG_bary(cbind(edge_index, 1))
              end_seg[j, ] <- line_MGG[i + 1, ]
            }
          }
        } else if (any(end_on_vertex)) {
          # start is not a vertex, end is a vertex
          end_vertex <-
            fm_MG_graph(graph)$E[line_MGG$index[i + 1], end_on_vertex]
          # not any start_on_vertex
          edge_index <- line_MGG$index[i]
          edge_nodes <- fm_MG_graph(graph)$E[edge_index, ]
          if (!any(edge_nodes == end_vertex)) {
            stop(paste0(
              "Unclear geom_path: Not directly connected points for line",
              k, " and points ", i, " and ", i + 1, "."
            ))
          }
          if (edge_nodes[1L] == end_vertex) {
            j <- j + 1
            start_seg[j, ] <- line_MGG[i, ]
            end_seg[j, ] <- fm_as_MGG_bary(cbind(edge_index, 0))
          } else {
            j <- j + 1
            start_seg[j, ] <- line_MGG[i, ]
            end_seg[j, ] <- fm_as_MGG_bary(cbind(edge_index, 1))
          }
        } else { # !any(start_on_vertex) & !any(end_on_vertex)
          # are the points on connected edges?
          v1 <- fm_MG_graph(graph)$E[line_MGG$index[i], ]
          v2 <- fm_MG_graph(graph)$E[line_MGG$index[i + 1], ]
          if (sum(v1[1] %in% v2) == 1) {
            j <- j + 1
            start_seg[j, ] <- line_MGG[i, ]
            end_seg[j, ] <- fm_as_MGG_bary(cbind(line_MGG$index[i], 0))
            # and next interedge
            j <- j + 1
            start_seg[j, ] <- fm_as_MGG_bary(cbind(
              line_MGG$index[i + 1],
              c(0, 1)[v2 %in% v1[1]]
            ))
            end_seg[j, ] <- line_MGG[i + 1, ]
          } else if (sum(v1[2] %in% v2) == 1) {
            j <- j + 1
            start_seg[j, ] <- line_MGG[i, ]
            end_seg[j, ] <- fm_as_MGG_bary(cbind(line_MGG$index[i], 1))
            # and next interedge
            j <- j + 1
            start_seg[j, ] <- fm_as_MGG_bary(cbind(
              line_MGG$index[i + 1],
              c(0, 1)[v2 %in% v1[2]]
            ))
            end_seg[j, ] <- line_MGG[i + 1, ]
          } else {
            stop(paste0(
              "Unclear geom_path: Subsequent points on not ",
              "directly connected edges on path ",
              l + 1, " for points indexed by ", i, " and ", i + 1, "."
            ))
          }
        }
      }
    }
    if (j < NROW(start_seg)) {
      start_seg <- start_seg[seq_len(j), , drop = FALSE]
      end_seg <- end_seg[seq_len(j), , drop = FALSE]
    }

    # storage for the start_seq & end_seq as a well-defined path
    # they should all be on the same edge:
    path_MGG <- fm_as_MGG_intervals(
      start = start_seg,
      end = end_seg,
      graph = graph
    )


    l <- l + 1
    paths[[l]] <- path_MGG
    ids[[l]] <- rep(k, NROW(path_MGG))
  }
  paths <- do.call(dplyr::bind_rows, paths)
  ids <- unlist(ids)
  paths <- tibble::tibble(paths = paths, ID = ids)
  paths
}
