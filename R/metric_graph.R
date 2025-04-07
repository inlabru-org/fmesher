# inlabru MAPPER ----

# required packages/suggest packages
#' @rawNamespace S3method(inlabru::bru_get_mapper,rspde_metric_graph)
#' @rawNamespace S3method(inlabru::ibm_n, bm_metric_graph)
#' @rawNamespace S3method(inlabru::ibm_values, bm_metric_graph)
#' @rawNamespace S3method(inlabru::ibm_jacobian, bm_metric_graph)
#' @rawNamespace S3method(inlabru::bru_mapper, metric_graph)

#' @title Wrapper that calls bru_mapper with correct input
#' @param model Model class (contains a metric graph object)
#' @param \dots arguments passed to sub-methods
#' @rdname bru_get_mapper_rspde_metric_graph
bru_get_mapper.rspde_metric_graph <- function(model, ...) {
  if ((model[["f"]]$n) %% (fm_dof(model[["mesh"]])) != 0) {
    stop(paste0(
      "Incompatible degrees of freedom. SPDE: ",
      model[["f"]]$n, " and mesh: ", fm_dof(model["mesh"])
    ))
  }
  inlabru::bru_mapper(model[["mesh"]],
    n_rep = (model[["f"]]$n) / (fm_dof(model[["mesh"]]))
  )
}

#' @title bru_mapper for the metric_graph class
#' @param mesh a metric_graph object
#' @param n_rep number of components in linear predictor
#' @param \dots arguments passed to sub-methods
#' @rdname bm_metric_graph
bru_mapper.metric_graph <- function(mesh, n_rep = 1, ...) {
  mapper <- inlabru::bru_mapper_fmesher(mesh)
  if (n_rep > 1) {
    mapper <- inlabru::bru_mapper_repeat(mapper, n_rep = n_rep)
  }
  mapper
}

#' @describeIn bm_metric_graph Returns the degrees of freedom (number of
#'   vertices in the mesh)
#' @param mapper A `bm_metric_graph` object
ibm_n.bm_metric_graph <- function(mapper, ...) {
  mesh <- mapper[["mesh"]]
  n_rep <- mapper[["n_rep"]]
  return(n_rep * fmesher::fm_dof(mesh))
}
#' @describeIn bm_metric_graph Returns a vector with indices for the
#'   degrees of freedom
ibm_values.bm_metric_graph <- function(mapper, ...) {
  seq_len(inlabru::ibm_n(mapper))
}
#' @describeIn bm_metric_graph Returns the mapping matrix between
#' @param input Data input for the mapper
ibm_jacobian.bm_metric_graph <- function(mapper, input, ...) {
  mesh <- mapper[["mesh"]] # metric graph object
  n_rep <- mapper[["n_rep"]]
  if (is.null(input)) {
    return(Matrix::Matrix(0, 0, inlabru::ibm_n(mapper)))
  }
  # pte_tmp <- mesh$mesh$VtE
  # input_list <- lapply(seq_len(nrow(input)), function(i){input[i,]})
  # pte_tmp_list <- lapply(seq_len(nrow(pte_tmp)), function(i){pte_tmp[i,]})
  # idx_tmp <- match(input_list, pte_tmp_list)
  A_tmp <- fm_basis(mesh, input) # idx_tmp
  return(fm_row_kron(Matrix::Matrix(1, NROW(A_tmp), n_rep), A_tmp))
}

# fmesher functions ----

#' @title Internal helper functions for metric graph evaluation
#'
#' @description Methods called internally by [fm_basis()] methods.
#' @param x metric_graph object
#' @param loc Observation locations, can be either MGG coordinates, MGM
#'   coordinates or Euclidean coordinates (passed to fm_bary())
#' @param weights Optional weight vector, one weight for each location
#' @inheritParams fm_basis
#' @export
#' @keywords internal
#' @returns A `fm_basis` object; a list of evaluator information objects,
#' at least a matrix `A` and logical vector `ok`.
fm_basis.metric_graph <- function(x,
                                  loc,
                                  weights = NULL,
                                  ...,
                                  full = FALSE) {
  if (is.null(weights)) {
    weights <- rep(1.0, NROW(loc))
  } else if (length(weights) == 1L) {
    weights <- rep(weights, NROW(loc))
  }

  # derivatives <- !is.null(derivatives) && derivatives
  info <- list()
  # use metric graph function to get basis functions
  # obtain bary wrt to MGM
  barys <- fm_bary(x, loc, MGG = FALSE)
  n <- NROW(barys)
  info$A <- Matrix::sparseMatrix(
    i = c(seq_len(n), seq_len(n)),
    j = c(x$mesh$E[barys$index, 1], x$mesh$E[barys$index, 2]),
    x = c(weights * barys$where[, 1], weights * barys$where[, 2]),
    dims = c(n, fm_dof(x))
  )
  info[["ok"]] <- rep(TRUE, n)

  fm_basis(
    structure(
      info,
      class = "fm_basis"
    ),
    full = full
  )
}


#' @describeIn fm_bary Return a tibble with elements
#' `fm_bary`
#'
#' @param MGG indicator for the barycentric coordinates related to the graph
#'   (MGG) or mesh (MGM). Default is MGG coords
#' @export
fm_bary.metric_graph <- function(mesh,
                                 loc,
                                 MGG = TRUE,
                                 ...) {
  if (!MGG) {
    if (is.null(mesh$mesh)) {
      stop("There is no mesh.")
    }
  }
  if (inherits(loc, "fm_bary")) {
    if (inherits(loc, "graph")) { # TO DO: name this class
      if (MGG) {
        bary_coord <- loc
      } else {
        # bary_coord <- MGG_to_MGM(loc, mesh)
        bary_coord <- as_MGM(mesh$.__enclos_env__$private$PtE_to_mesh(
          cbind(loc$index, loc$where[, 2])
        ))
      }
    } else if (inherits(loc, "mesh")) {
      if (MGG) {
        bary_coord <- MGM_to_MGG(loc, mesh)
      } else {
        bary_coord <- loc
      }
    }
  } else if (inherits(loc, "sfg") || inherits(loc, "sf") ||
    inherits(loc, "sfc")) {
    # check the crs of point and convert to the same crs as graph (or
    # coordinates handles this)
    res <- Euclidean_to_graph(sf::st_coordinates(loc), mesh)
    bary_coord <- res # res$bary[res$ok, ]
    if (!MGG) {
      # bary_coord <- MGG_to_MGM(bary_coord, mesh)
      bary_coord <- as_MGM(mesh$.__enclos_env__$private$PtE_to_mesh(
        cbind(bary_coord$index, bary_coord$where[, 2])
      ))
    }
  } else {
    # Or should it only call as_MGG/as_MGM depending on "MGG"?
    # cat("loc is interpreted as Euclidean coordinates")
    res <- Euclidean_to_graph(loc, mesh)
    if (MGG) {
      bary_coord <- res # res$bary
    } else {
      # bary_coord <- MGG_to_MGM(res, mesh) # MGG_to_MGM(res$bary, mesh)
      bary_coord <- as_MGG(mesh$.__enclos_env__$private$PtE_to_mesh(
        cbind(res$index, res$where[, 2])
      ))
    }
  }
  return(bary_coord)
}


#' @rdname fm_manifold
#' @export
fm_manifold_get.metric_graph <- function() {
  return("G1")
}


#' @rdname fm_dof
#' @export
fm_dof.metric_graph <- function(x) {
  NROW(x[["mesh"]][["VtE"]])
}


#' @export
#' @describeIn fm_int `metric_graph` integration. Supported samplers: * `NULL`
#'   for integration over the entire domain; * A tibble with a named column
#'   containing a matrix with single edge intervals (ordered), and optionally a
#'   `weight` column.
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
#'     start_MGG = cbind(1, 0.2),
#'     edges = c(2),
#'     end_MGG = cbind(3, 0.8)
#'   )
#'   samplers <- tibble::tibble(x = list(p1), weight = c(1))
#'   ips <- fm_int(
#'     graph,
#'     samplers
#'   )
#' }
fm_int.metric_graph <- function(domain,
                                samplers = NULL,
                                name = "x",
                                ...) {
  ips <- list()
  if (is.null(domain$mesh)) {
    stop("There is no mesh")
  }

  mesh_MGG <- as_MGG(domain$mesh$VtE) # All mesh locations in MGG format
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
      if (!inherits(interedge, "graph_interval")) {
        interedge <- as_graph_interval(
          graph = domain,
          start_MGG = as_MGG(list(
            interedge$start$index,
            interedge$start$where[1, 2]
          )),
          end_MGG = as_MGG(list(
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
        inside <- (loc_mid <= interedge$start$where[1, 2]) & (loc_mid >= interedge$end$where[1, 2])
      } else {
        inside <- (loc_mid >= interedge$start$where[1, 2]) & (loc_mid <= interedge$end$where[1, 2])
      }
      # convert to MGM (call outside of for-loop and only get the desired rows in this step)
      # loc_mid_MGM <- MGG_to_MGM(
      #   coord = as_MGG(cbind(interedge$start$index, loc_mid)),
      #   graph = domain
      # )
      loc_mid_MGM <- domain$.__enclos_env__$private$PtE_to_mesh(
        cbind(interedge$start$index, loc_mid)
      )
      loc_mid_MGM <- as_MGM(loc_mid_MGM)
      # get the edge lengths for each mesh
      weight_mid <- domain$mesh$h_e[loc_mid_MGM$index]
      weight_mid[!inside] <- 0.0

      weight_trap <- c(weight_mid / 2, 0) + c(0, weight_mid / 2)
      loc_simpson <- c(loc_trap, loc_mid)
      weight_simpson <- c(weight_trap / 3, weight_mid * 2 / 3)

      m_ips <- sum(weight_simpson > 0)
      if (m_ips == 0) {
        ips_edge[[k]] <- tibble::tibble(
          x = as_MGG(tibble::tibble(
            index = integer(0),
            where = numeric(0)
          )),
          weight = numeric(0),
          .block = integer(0)
        )
      } else {
        ips_edge[[k]] <- tibble::tibble(
          x = as_MGG(loc = cbind(
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
      x = numeric(0),
      weight = numeric(0),
      .block = integer(0)
    )
    colnames(ips)[1] <- name
  }
  ips
}




# MetricGraph specific functions----
#' @title Make a ("graph", "fm_bary") object from Euclidean coordinates
#' @description
#' Create a (`graph`, `fm_bary`) object from Euclidean coordinates.
#'
#' @param loc Euclidean coords (if not on graph, they are mapped to the closest
#'   point on graph)
#' @param graph metric_graph that the location should be mapped to.
#' @author Karina Lilleborge \email{karina.lilleborge@@gmail.com}
#' @returns A (`mesh`, `fm_bary`) object
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
#'   m <- Euclidean_to_graph(
#'     cbind(0, 1),
#'     graph
#'   )
#'   # c(2,1)
#'   m
#' }
#'
Euclidean_to_graph <- function(loc, graph) {
  res <- graph$coordinates(XY = loc)
  # check distance from original points:
  # tolerance <- min(graph$edge_lengths) / 2
  # tmp_loc <- graph$coordinates(PtE = res, normalized = TRUE)
  # norm_XY <- sf::st_distance(loc, sf::st_point(tmp_loc, crs= sf::st_crs(loc)))
  # ok_ <- (norm_XY < tolerance)
  # graph_coords <- as_MGG(loc = res)
  # return(list(bary = graph_coords, ok = ok_))
  return(as_MGG(loc = res))
}

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
#'     as_MGG(cbind(1, 0.5)),
#'     graph
#'   )
#'   mgm
#' }
#'
MGG_to_MGM <- function(coord, graph) {
  if (is.null(graph$mesh)) {
    stop("There is no mesh")
  }
  stopifnot(inherits(coord, "graph") && inherits(coord, "fm_bary"))
  mesh_MGG <- as_MGG(graph$mesh$VtE)
  mesh_edge_len <- graph$mesh$h_e
  # storage for new coordinates
  new_coord <- as_MGM(tibble::tibble(
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
      vertices_MGG <- graph$E[coord$index[i], ]

      index_MGM <- which.max((graph$mesh$E[, 1] == vertices_MGG[1]) &
        (graph$mesh$E[, 2] == vertices_MGG[2]))
      where_MGM <- coord$where[i, 2]
    } else if (sum(ids) == 1) {
      # there is only one mesh vertex on the edge
      index_on_edge <- which(ids)
      # coord[i, ] is before or past the mesh node
      further <- ((as.numeric(coord$where[i, 2]) - edge_MGG$where[2]) >= 0)
      # find the vertex on the other side of coord[i, ]
      graph_vertex <- graph$E[edge_MGG$index, further * 1 + 1]
      if (further) {
        # find the edge index that connects (mesh_vertex, end_vertex)
        index_MGM <- which.max((graph$mesh$E[, 1] == index_on_edge) &
          (graph$mesh$E[, 2] == graph_vertex))
        mesh_h_e <- mesh_edge_len[index_MGM]
        where_MGM <- as.numeric((as.numeric(coord$where[i, 2]) - edge_MGG$where[, 2]) / mesh_h_e)
      } else {
        # find the edge index that connects (start_vertex, mesh_vertex)
        index_MGM <- which.max((graph$mesh$E[, 1] == graph_vertex) &
          (graph$mesh$E[, 2] == index_on_edge))
        mesh_h_e <- mesh_edge_len[index_MGM]
        where_MGM <- 1 - (as.numeric((edge_MGG$where[, 2] - as.numeric(coord$where[i, 2])) / mesh_h_e))
      }
    } else {
      # order the mesh_MGG locations:
      ordering <- order(edge_MGG$where[, 2])
      edge_MGG_o <- edge_MGG[ordering, ] # unique(c(0, edge_MGG[ordering, ], 1)) # EDIT 21.01 from: edge_MGG[ordering, ]
      # edge_MGG_o <- rbind(tibble::tibble(index=coord$index[i], where=cbind(1,0)),
      #                     edge_MGG[ordering, ],
      #                     tibble::tibble(index=coord$index[i], where=cbind(0,1))) #these are ordered #edge_MGG[ordering, ]
      # #remove duplicates:
      # edge_MGG_o <- edge_MGG_o[!dplyr::duplicated(edge_MGG_o$where[,2]), ]
      # find the mesh point index where coord[i,] is next to
      index_on_edge <- which.max((edge_MGG_o$where[, 2] - as.numeric(coord$where[i, 2])) >= 0)
      if (index_on_edge) {
        # coord[i, ] is between these two mesh locs
        mesh_indices <- which(ids)[(ordering[c(index_on_edge - 1, index_on_edge)])]
      }
      index_MGM <- which.max((graph$mesh$E[, 1] == mesh_indices[1]) &
        (graph$mesh$E[, 2] == mesh_indices[2]))
      mesh_h_e <- mesh_edge_len[index_MGM]
      where_MGM <- 1 - as.numeric((edge_MGG_o$where[index_on_edge, 2] - as.numeric(coord$where[i, 2])) / mesh_h_e) # normalized
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
#'     as_MGM(cbind(5, 1)),
#'     graph
#'   )
#'   mgg
#' }
#'
MGM_to_MGG <- function(coord, graph) {
  stopifnot(inherits(coord, "mesh") && inherits(coord, "fm_bary"))
  mesh_loc <- graph$mesh$VtE
  eps <- min(graph$mesh$h_e)
  new_coord <- as_MGG(tibble::tibble(
    index = integer(NROW(coord)),
    where = numeric(NROW(coord))
  ))
  for (i in seq_len(NROW(coord))) {
    # mesh vertex indices for neighboring loc
    mesh_vertices <- graph$mesh$E[coord$index[i], ]
    # which graph coordinates for the neighboring vertices
    graph_edge_r <- mesh_loc[mesh_vertices[2], ]
    graph_edge_l <- mesh_loc[mesh_vertices[1], ]
    # as long as they are on the same edge
    if (graph_edge_l[1] == graph_edge_r[1]) {
      # neighbouring mesh nodes are on the same edge
      new_coord[i, ] <- as_MGG(
        tibble::tibble(
          index = graph_edge_l[1],
          where = (1 - coord$where[i, 2]) * graph_edge_l[2] + coord$where[i, 2] * graph_edge_r[2]
        )
      )
    } else {
      # neighbouring mesh nodes are not on same edge (one is a vertex)
      on_vertex_r <- (abs(graph_edge_r[2] - c(0, 1)) < eps)
      on_vertex_l <- (abs(graph_edge_l[2] - c(0, 1)) < eps)
      if (any(on_vertex_r)) {
        new_coord[i, ] <- as_MGG(
          tibble::tibble(
            index = graph_edge_l[1],
            where = (1 - coord$where[i, 2]) * graph_edge_l[2] + coord$where[i, 2] * c(0, 1)[on_vertex_r]
          )
        )
      } else {
        new_coord[i, ] <- as_MGG(
          tibble::tibble(
            index = graph_edge_r[1],
            where = (1 - coord$where[i, 2]) * c(0, 1)[on_vertex_l] + coord$where[i, 2] * graph_edge_r[2]
          )
        )
      }
    }
  }
  as_MGG(new_coord)
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
#'   m <- as_MGM(
#'     cbind(1, 1),
#'     graph
#'   )
#'   class(m) # "mesh", "fm_bary", "tbl_df", "tbl", "data.frame"
#' }
#'
as_MGM <- function(loc, graph = NULL) {
  if (inherits(loc, "mesh") && inherits(loc, "fm_bary")) {
    return(loc)
  }
  if (inherits(loc, "graph") && inherits(loc, "fm_bary")) {
    if (is.null(graph)) {
      stop("Graph must be provided to convert from MGG to MGM.")
    }
    res <- (graph$.__enclos_env__$private$PtE_to_mesh(cbind(loc$index, loc$where[, 2]))) # MGG_to_MGM(coord = loc, graph = graph)
    res <- tibble::tibble(
      index = res[, 1],
      where = cbind(
        1 - as.numeric(res[, 2]),
        as.numeric(res[, 2])
      )
    )
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
      class = c("mesh", "fm_bary", "tbl_df", "tbl", "data.frame")
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
#'   m <- as_MGG(cbind(1, 0.5))
#'   class(m) # "graph", "fm_bary", "tbl_df", "tbl", "data.frame"
#' }
#'
as_MGG <- function(loc, graph = NULL) {
  if (inherits(loc, "graph") && inherits(loc, "fm_bary")) {
    return(loc)
  }
  if (inherits(loc, "mesh") && inherits(loc, "fm_bary")) {
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
      class = c("graph", "fm_bary", "tbl_df", "tbl", "data.frame")
    )
}


#' @title Make an inter edge interval on graph object
#' @description
#' Create an `graph_interval` object.
#'
#' @param start_MGG Start location for inter edge interval.
#' @param end_MGG End location for inter edge interval.
#' @param graph `metric_graph` that the interval should be mapped to. Must be
#'   provided if input should be converted
#' @author Karina Lilleborge \email{karina.lilleborge@@gmail.com}
#' @returns An `graph_interval` object
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
#'   int <- as_graph_interval(
#'     cbind(1, 0.8),
#'     cbind(1, 0.5),
#'     graph
#'   )
#'   int
#' }
#'
as_graph_interval <- function(start_MGG,
                              end_MGG,
                              graph = NULL) {
  if (!(inherits(start_MGG, "graph") && inherits(start_MGG, "fm_bary"))) {
    start_MGG <- as_MGG(start_MGG, graph = graph)
  }
  if (!(inherits(end_MGG, "graph") && inherits(end_MGG, "fm_bary"))) {
    end_MGG <- as_MGG(end_MGG, graph = graph)
  }
  if (!(sum(start_MGG$index == end_MGG$index) == NROW(start_MGG))) {
    stop("Not all start- and end points are inter edge intervals")
  }
  inter_edge_interval <- structure(
    tibble::tibble(
      start = start_MGG,
      end = end_MGG
    ),
    class = c("graph_interval", "tbl_df", "tbl", "data.frame")
  )
}

#' @title Make an interval on graph object
#' @description
#' Create a `graph_interval` object from known start (MGG), end (MGG) and
#' visiting edges (MGG).
#'
#' @param graph metric_graph that the interval should be mapped to.
#' @param start_MGG MGG coordinates for start
#' @param edges Ordered list of edge indices related to MGG
#' @param end_MGG MGG coordinates for end
#' @author Karina Lilleborge \email{karina.lilleborge@@gmail.com}
#' @returns A `graph_interval` object
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
#'     start_MGG = cbind(1, 0.5),
#'     edges = c(2),
#'     end_MGG = cbind(3, 0.6)
#'   )
#'   path
#' }
#'
simple_path_MGG <- function(graph,
                            start_MGG,
                            edges,
                            end_MGG) {
  # check if the graph does have circles
  start_MGG <- as_MGG(start_MGG)
  end_MGG <- as_MGG(end_MGG)
  if (length(edges) > 0) {
    # check direction from start to edges[1]
    v1 <- graph$E[as.integer(start_MGG$index), ]
    v2 <- graph$E[as.integer(edges[1L]), ]
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
      as_graph_interval(
        as_MGG(tibble::tibble(
          index = integer(n),
          where = numeric(n)
        )),
        as_MGG(tibble::tibble(
          index = integer(n),
          where = numeric(n)
        ))
      )
    inter_edge_intervals[1, ] <- tibble::tibble(
      start = start_MGG,
      end = as_MGG(
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
        start = as_MGG(cbind(as.integer(edges[i]), as.numeric(start_vertex))),
        end = as_MGG(cbind(as.integer(edges[i]), as.numeric(end_vertex)))
      )
      v1 <- graph$E[as.integer(edges[i]), ]
      v2 <- graph$E[as.integer(edges[i + 1]), ]
    }
    v2 <- graph$E[as.integer(end_MGG$index), ]
    # check direction
    start_vertex <- c(0:1)[(v2 %in% v1[end_vertex + 1])]
    inter_edge_intervals[length(edges) + 2, ] <- tibble::tibble(
      start = as_MGG(cbind(
        as.integer(end_MGG$index),
        as.numeric(start_vertex)
      )),
      end = as_MGG(cbind(
        as.integer(end_MGG$index),
        as.numeric(end_MGG$where[, 2L])
      ))
    )
  } else { # there are no whole edges visited (edges=c())
    if (as.integer(start_MGG$index) == as.integer(end_MGG$index)) { # same edge
      inter_edge_intervals <- as_graph_interval(
        start_MGG = as_MGG(tibble::tibble(
          index = integer(1),
          where = numeric(1)
        )),
        end_MGG = as_MGG(tibble::tibble(
          index = integer(1),
          where = numeric(1)
        ))
      )
      inter_edge_intervals[1, ] <- tibble::tibble(
        start = as_MGG(cbind(
          as.integer(start_MGG$index),
          as.numeric(start_MGG$where[, 2])
        )),
        end = as_MGG(cbind(
          as.integer(start_MGG$index),
          as.numeric(end_MGG$where[, 2])
        ))
      )
    } else { # neighboring edges
      v1 <- graph$E[as.integer(start_MGG$index), ]
      v2 <- graph$E[as.integer(end_MGG$index), ]
      if (sum(v2 %in% v1[1]) > 0) {
        # if the (e,0) vertex is in v2
        end_vertex <- 0
      }
      if (sum(v2 %in% v1[2]) > 0) {
        # if the (e,1) vertex is in v2
        end_vertex <- 1
      }
      # make storage for the inter edge intervals for each of the edge
      # index, start and end (graph_interval)
      inter_edge_intervals <- as_graph_interval(
        start_MGG = as_MGG(tibble::tibble(
          index = integer(2),
          where = numeric(2)
        )),
        end_MGG = as_MGG(tibble::tibble(
          index = integer(2),
          where = numeric(2)
        ))
      )
      inter_edge_intervals[1, ] <- tibble::tibble(
        start = as_MGG(cbind(
          as.integer(start_MGG$index),
          as.numeric(start_MGG$where[, 2])
        )),
        end = as_MGG(cbind(
          as.integer(start_MGG$index),
          as.numeric(end_vertex)
        ))
      )
      v2 <- graph$E[as.integer(end_MGG$index), ]
      # check direction
      start_vertex <- c(0:1)[(v2 %in% v1[end_vertex + 1L])]
      inter_edge_intervals[2, ] <- tibble::tibble(
        start = as_MGG(cbind(
          as.integer(end_MGG$index),
          as.numeric(start_vertex)
        )),
        end = as_MGG(cbind(
          as.integer(end_MGG$index),
          as.numeric(end_MGG$where[, 2])
        ))
      )
    }
  }
  # construct object
  path <- structure(
    inter_edge_intervals,
    class = c("graph_interval", "tbl_df", "tbl", "data.frame")
  )
  return(path)
}

#' @title Make an interval on graph object from sf object
#' @description
#' Create a tibble of `graph_interval` objects from `sf::st_geometry`
#' (`LINESTRING`) and a column with "id".
#'
#' @param graph metric_graph that the interval should be mapped to.
#' @param geom_path `sf::st_geometry` (`LINESTRING`) on a graph
#' @author Karina Lilleborge \email{karina.lilleborge@@gmail.com}
#' @returns A tibble containing a set of `graph_interval` objects and id
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
geom_path_to_path_MGG <- function(geom_path, graph) {
  # new function name for this (as_graph_interval(input) check what input is)
  if (!inherits(geom_path, "sfc_LINESTRING")) {
    stop("Method not implemented. Input must be sfc_LINESTRING")
  }
  # TO DO: convert to correct crs: Handled by MetricGraph$coordinates.
  # create eps
  if (is.null(graph$mesh)) {
    stop("The graph has no mesh!")
  }
  eps <- 0.01 * min(graph$mesh$h_e)
  # matrix with colnames X Y and L1:
  internal_XY <- sf::st_coordinates(geom_path)
  paths <- list()
  ids <- list()
  l <- 0
  for (k in unique(internal_XY[, "L1"])) {
    # a line should give us one path
    line <- internal_XY[internal_XY[, "L1"] == k, ]
    line_MGG <- fm_bary(graph, as.matrix(line[, c("X", "Y")]), MGG = TRUE)
    # line_MGG <- graph$coordinates(XY = line[, c("X", "Y")])
    # convert to ("graph", "fm_bary")
    # line_MGG <- as_MGG(line_MGG)
    # index for number of segments added
    j <- 0
    # storing the segments (start and end separately)
    start_seg <- as_MGG(cbind(
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
          start_vertex <- graph$E[line_MGG$index[i], start_on_vertex]

          if (any(end_on_vertex)) {
            # the end point is also a vertex
            # the end is a vertex -> find which index
            end_vertex <- graph$E[line_MGG$index[i + 1], end_on_vertex]
            if (start_vertex != end_vertex) {
              # they are not the same vertex
              edges_start <- graph$E == start_vertex
              edges_end <- graph$E == end_vertex
              # find what edge both start vertex and end vertex are present
              edges <- rowSums(edges_start) > 0 & rowSums(edges_end) > 0
              if (sum(edges) > 1) {
                stop("Ambiguous geom_path: Multiple edge candidates.")
              }
              if (sum(edges) == 0) {
                stop("Unclear geom_path: Not directly connected subsequent end points.")
              }
              # there is only one edge match
              edge_index <- which(edges)
              if (isTRUE(edges_start[edge_index, 1])) {
                # the start vertex is the start of the edge v = (edge_index, 0)
                j <- j + 1
                start_seg[j, ] <- as_MGG(cbind(edge_index, 0))
                end_seg[j, ] <- as_MGG(cbind(edge_index, 1))
              } else {
                # the start is the end of the edge v = (edge_index, 1)
                j <- j + 1
                start_seg[j, ] <- as_MGG(cbind(edge_index, 1))
                end_seg[j, ] <- as_MGG(cbind(edge_index, 0))
              }
            }
          } else {
            # !any(end_on_vertex) (the end is not on a vertex, but the start is)
            # we check that the start is not ambiguous wrt the end point
            edge_index <- line_MGG$index[i + 1]
            # candidates for the start
            edge_nodes <- graph$E[edge_index, ]
            if (!any(edge_nodes == start_vertex)) {
              stop(paste0(
                "Unclear geom_path: Not directly connected points for line",
                k, " and points ", i, " and ", i + 1, "."
              ))
            }
            # we have a valid start
            if (edge_nodes[1L] == start_vertex) {
              j <- j + 1
              start_seg[j, ] <- as_MGG(cbind(edge_index, 0))
              end_seg[j, ] <- line_MGG[i + 1, ]
            } else {
              j <- j + 1
              start_seg[j, ] <- as_MGG(cbind(edge_index, 1))
              end_seg[j, ] <- line_MGG[i + 1, ]
            }
          }
        } else if (any(end_on_vertex)) { # start is not a vertex, end is a vertex

          end_vertex <- graph$E[line_MGG$index[i + 1], end_on_vertex]
          # not any start_on_vertex
          edge_index <- line_MGG$index[i]
          edge_nodes <- graph$E[edge_index, ]
          if (!any(edge_nodes == end_vertex)) {
            stop(paste0(
              "Unclear geom_path: Not directly connected points for line",
              k, " and points ", i, " and ", i + 1, "."
            ))
          }
          if (edge_nodes[1L] == end_vertex) {
            j <- j + 1
            start_seg[j, ] <- line_MGG[i, ]
            end_seg[j, ] <- as_MGG(cbind(edge_index, 0))
          } else {
            j <- j + 1
            start_seg[j, ] <- line_MGG[i, ]
            end_seg[j, ] <- as_MGG(cbind(edge_index, 1))
          }
        } else { # !any(start_on_vertex) & !any(end_on_vertex)
          # are the points on connected edges?
          v1 <- graph$E[line_MGG$index[i], ]
          v2 <- graph$E[line_MGG$index[i + 1], ]
          if (sum(v1[1] %in% v2) == 1) {
            j <- j + 1
            start_seg[j, ] <- line_MGG[i, ]
            end_seg[j, ] <- as_MGG(cbind(line_MGG$index[i], 0))
            # and next interedge
            j <- j + 1
            start_seg[j, ] <- as_MGG(cbind(line_MGG$index[i + 1], c(0, 1)[v2 %in% v1[1]]))
            end_seg[j, ] <- line_MGG[i + 1, ]
          } else if (sum(v1[2] %in% v2) == 1) {
            j <- j + 1
            start_seg[j, ] <- line_MGG[i, ]
            end_seg[j, ] <- as_MGG(cbind(line_MGG$index[i], 1))
            # and next interedge
            j <- j + 1
            start_seg[j, ] <- as_MGG(cbind(line_MGG$index[i + 1], c(0, 1)[v2 %in% v1[2]]))
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
    path_MGG <- as_graph_interval(
      start_MGG = start_seg,
      end_MGG = end_seg,
      graph = graph
    )


    l <- l + 1
    paths[[l]] <- path_MGG
    ids[[l]] <- rep(k, NROW(path_MGG))
  }
  paths <- do.call(dplyr::bind_rows, paths)
  ids <- unlist(ids)
  paths <- tibble::tibble(paths = paths, ID = ids)

  return(paths)
}
