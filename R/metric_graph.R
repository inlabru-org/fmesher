# MAPPER ----

# inlabru functions ----
## INFO FROM inlabru: ----
### Constructor ----
# @param mapper For `bru_mapper_define`, a prototype mapper object, see
#   Details. For `bru_mapper_scale`, a mapper to be scaled.
# @param new_class If non-`NULL`, this is added at the front of the class
#   definition
#
# @describeIn bru_mapper Adds the `new_class` and "bru_mapper" class names to
#   the inheritance list for the input `mapper` object, unless the object
#   already inherits from these.
#
# To register mapper classes and methods in scripts, use `.S3method()`
# to register the methods, e.g.
# `.S3method("ibm_jacobian", "my_mapper_class", ibm_jacobian.my_mapper_class)`.
#
# In packages with `Suggests: inlabru`, add method information for delayed
# registration, e.g.:
# ```
# # @rawNamespace S3method(inlabru::bru_get_mapper, inla_rspde)
# # @rawNamespace S3method(inlabru::ibm_n, bru_mapper_inla_rspde)
# # @rawNamespace S3method(inlabru::ibm_values, bru_mapper_inla_rspde)
# # @rawNamespace S3method(inlabru::ibm_jacobian, bru_mapper_inla_rspde)
# ```
# or before each method, use `@exportS3Method`:
# ```
# # @exportS3Method inlabru::bru_get_mapper
# ```
# etc., which semi-automates it.


# styler::style_pkg()

# see functions:
# metric_graph$public_methods$function_name

# required packages/suggest packages


# @title Wrapper that calls bru_mapper with correct input
# @param model Model class (contains a metric graph object)
# bru_get_mapper.inla_rspde <- function(model,...){
#   bru_mapper(model[["mesh"]])
# }

#' @rawNamespace S3method(inlabru::ibm_n, bru_mapper_metric_graph)
#' @rawNamespace S3method(inlabru::ibm_values, bru_mapper_metric_graph)
#' @rawNamespace S3method(inlabru::ibm_jacobian, bru_mapper_metric_graph)
#' @rawNamespace S3method(inlabru::bru_mapper, metric_graph)



#' @title bru_mapper for the metric_graph class
#' @param mesh a metric_graph object
#' @param n_eta number of components in linear predictor
#' @param \dots arguments passed to sub-methods
#' @rdname bru_mapper_metric_graph
bru_mapper.metric_graph <- function(mesh, n_eta = 1, ...) {
  mapper <- list(mesh = mesh, n_eta = n_eta)
  inlabru::bru_mapper_define(mapper, new_class = "bru_mapper_metric_graph")
}


###### ------

#' @describeIn bru_mapper_metric_graph Returns the degrees of freedom (number of vertices in the mesh)
#' @param mapper A `bru_mapper_metric_graph` object
ibm_n.bru_mapper_metric_graph <- function(mapper, ...) {
  mesh <- mapper[["mesh"]]
  n_eta <- mapper[["n_eta"]]
  # should return the degrees of freedom
  return(n_eta * fmesher::fm_dof(mesh))
}
#' @describeIn bru_mapper_metric_graph Returns a vector with indices for the degrees of freedom
ibm_values.bru_mapper_metric_graph <- function(mapper, ...) {
  seq_len(inlabru::ibm_n(mapper))
}
#' @describeIn bru_mapper_metric_graph Returns the mapping matrix between
#' @param input Data input for the mapper
ibm_jacobian.bru_mapper_metric_graph <- function(mapper, input, ...) {
  mesh <- mapper[["mesh"]] # metric graph object
  n_eta <- mapper[["n_eta"]]
  if (is.null(input)) {
    return(Matrix::Matrix(0, 0, inlabru::ibm_n(mapper)))
  }
  # pte_tmp <- mesh$mesh$VtE
  # input_list <- lapply(seq_len(nrow(input)), function(i){input[i,]})
  # pte_tmp_list <- lapply(seq_len(nrow(pte_tmp)), function(i){pte_tmp[i,]})
  # idx_tmp <- match(input_list, pte_tmp_list)
  A_tmp <- fm_basis(mesh, input, n_eta = n_eta) # idx_tmp
  return(A_tmp)
}

# fmesher functions ----

#' @title Internal helper functions for metric graph evaluation
#'
#' @description Methods called internally by [fm_basis()] methods.
#' @param x metric_graph object
#' @param loc Observation locations, can be either
#'            MGG coordinates, MGM coordinates or Euclidean coordinates (passed to fm_bary())
#' @param weights Optional weight vector, one weight for each location
#' @inheritParams fm_basis
#' @export
#' @keywords internal
#' @returns A `fm_basis` object; a list of evaluator information objects,
#' at least a matrix `A` and logical vector `ok`.
fm_basis.metric_graph <- function(x,
                                  loc,
                                  weights = NULL,
                                  n_eta = 1,
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
    x = c(weights * (1 - barys$where), weights * barys$where),
    dims = c(n, n_eta * fm_dof(x))
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
  if (is.null(mesh$mesh)) {
    if (!MGG) {
      stop("There is no mesh")
    }
  }
  if (inherits(loc, "fm_bary_MGG")) {
    if (MGG) {
      bary_coord <- loc
    } else {
      bary_coord <- graph_to_mesh_coord(mesh, loc)
    }
  } else if (inherits(loc, "fm_bary_MGM")) {
    if (MGG) {
      bary_coord <- mesh_to_graph_coord(mesh, loc)
    } else {
      bary_coord <- loc
    }
  } else {
    cat("loc is interpreted as Euclidean coordinates")
    res <- Euclidean_to_MGG(mesh, loc)
    if (MGG) {
      bary_coord <- res
    } else {
      bary_coord <- graph_to_mesh_coord(mesh, res)
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
  NROW(x$mesh$V)
}


#' @export
#' @describeIn fm_int `metric_graph` integration. Supported samplers:
#' * `NULL` for integration over the entire domain;
#' * A tibble with a named column containing a matrix with single edge intervals (ordered), and optionally a
#'  `weight` column.
#' @examples
#' if (requireNamespace("MetricGraph")) {
#'   edge1 <- rbind(c(0, 0), c(1, 0))
#'   edge2 <- rbind(c(0, 0), c(0, 1))
#'   edge3 <- rbind(c(0, 1), c(-1, 1))
#'   theta <- seq(from = pi, to = 3 * pi / 2, length.out = 20)
#'   edge4 <- cbind(sin(theta), 1 + cos(theta))
#'   edges <- list(edge1, edge2, edge3, edge4)
#'   graph <- MetricGraph::metric_graph$new(edges = edges)
#'   graph$build_mesh(h=0.01)
#'   p1 <- path_MGG(
#'     graph = graph,
#'     start_MGG = matrix(c(1, 0.2), nrow = 1),
#'     edges = c(2),
#'     end_MGG = matrix(c(3, 0.8), nrow = 1)
#'   )
#'   samplers <- tibble::tibble(x = list(p1), weight = c(1))
#'   #ips <- fm_int(
#'   #  graph,
#'   #  samplers
#'   #)
#' }
fm_int.metric_graph <- function(domain, samplers = NULL, name = "x", int.args = NULL, ...) {
  int.args.default <- list(method = "stable", nsub1 = 30, nsub2 = 9)
  if (is.null(int.args)) {
    int.args <- list()
  }
  missing.args <- setdiff(names(int.args.default), names(int.args))
  int.args[missing.args] <- int.args.default[missing.args]
  if (!is.null(int.args[["nsub"]])) {
    int.args[["nsub1"]] <- int.args[["nsub"]]
  }
  ips <- list()
  if (is.null(domain$mesh)) {
    stop("There is no mesh")
  }

  mesh_locs <- domain$mesh$VtE # All mesh locations in MGG format
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

  for (j in seq_len(length(samplers))) {
    # for a graph interval
    subsampler <- samplers[[name]][[j]]
    theweight <- samplers[["weight"]][[j]]
    ips_edge <- list()
    for (k in seq_len(nrow(subsampler))) {
      interedge <- subsampler[k, , drop = TRUE] # class=MGG_interval
      if (!inherits(interedge, "MGG_interval")) {
        interedge <- MGG_interval(
          domain,
          matrix(c(interedge$index, interedge$where[1, 1]), nrow = 1),
          matrix(c(interedge$index, interedge$where[1, 2]), nrow = 1)
        )
      }
      the.block <- .block[j]
      # Simpson's rule integration
      # mesh vertices on edge of interest (+ end points)
      loc_trap <- sort(unique(c(0, mesh_locs[mesh_locs[, 1] == interedge$index, 2], 1)))
      loc_mid <- (loc_trap[-1] + loc_trap[-length(loc_trap)]) / 2
      # Detect mid-points inside the interval
      if (interedge$where[1, 1] > interedge$where[1, 2]) {
        inside <- (loc_mid <= interedge$where[1, 1]) & (loc_mid >= interedge$where[1, 2])
      } else {
        inside <- (loc_mid >= interedge$where[1, 1]) & (loc_mid <= interedge$where[1, 2])
      }
      # convert to MGM
      loc_mid_MGM <- graph_to_mesh_coord(domain, matrix(c(
        rep(interedge$index, NROW(loc_mid)),
        loc_mid
      ), nrow = NROW(loc_mid)))
      # get the edge lengths for each mesh
      weight_mid <- domain$mesh$h_e[loc_mid_MGM$index]
      weight_mid[!inside] <- 0.0

      weight_trap <- c(weight_mid / 2, 0) + c(0, weight_mid / 2)
      loc_simpson <- c(loc_trap, loc_mid)
      weight_simpson <- c(weight_trap / 3, weight_mid * 2 / 3)

      m_ips <- sum(weight_simpson > 0)
      ips_edge[[k]] <- tibble::tibble(
        x = fm_bary_MGG(domain,
          loc = matrix(c(
            rep(interedge$index, m_ips),
            loc_simpson[(weight_simpson > 0)]
          ), ncol = 2)
        ),
        weight = weight_simpson[(weight_simpson > 0)] * theweight,
        .block = the.block
      )
      colnames(ips_edge[[k]])[1] <- name
    }
    ips_edge <- do.call(dplyr::bind_rows, ips_edge)
    ips[[j]] <- ips_edge
  }
  ips <- do.call(dplyr::bind_rows, ips)
  if (NROW(ips) == 0) {
    ips <- tibble::tibble(x = numeric(0), weight = numeric(0), .block = integer(0))
    colnames(ips)[1] <- name
  }
  ips
}




# object creation and conversion----
#' @title Make a fm_bary_MGG object from Euclidean coordinates
#' @description
#' Create a `fm_bary_MGG` object from Euclidean coordinates.
#'
#' @param graph metric_graph that the location should be mapped to.
#' @param loc Euclidean coords (if not on graph, they are mapped to the closest point on graph)
#' @author Karina Lilleborge \email{karina.lilleborge@@gmail.com}
#' @returns An `fm_bary_MGM` object
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
#'   m <- Euclidean_to_MGG(graph,
#'     matrix(c(0, 1),nrow=1)
#'   )
#'   # c(2,1)
#'   m
#' }
#'
Euclidean_to_MGG <- function(graph, loc) {
  res <- graph$coordinates(XY = loc)
  graph_coords <- fm_bary_MGG(graph = graph, loc = res)
  return(graph_coords)
}

#' @title Make a fm_bary_MGM object from MGG coordinates
#' @description
#' Create a `fm_bary_MGM` object from MGG coordinates.
#'
#' @param graph metric_graph that the location should be mapped to.
#' @param coord MGG coordinates
#' @author Karina Lilleborge \email{karina.lilleborge@@gmail.com}
#' @returns An `fm_bary_MGM` object
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
#'   graph$build_mesh(h=0.005)
#'   mgm <- graph_to_mesh_coord(
#'     graph,
#'     matrix(c(1, 0.5),nrow=1)
#'   )
#'   mgm
#' }
#'
graph_to_mesh_coord <- function(graph,
                                coord) {
  if (is.null(graph$mesh)) {
    # error
    stop("There is no mesh")
  }
  mesh_loc <- graph$mesh$VtE
  mesh_edge_len <- graph$mesh$h_e
  new_coord <- matrix(nrow = NROW(coord), ncol = NCOL(coord))
  for (i in seq_len(NROW(coord))) {
    # which mesh edges are in the given edge
    ids <- (mesh_loc[, 1] == as.integer(coord[i, 1]))
    mesh_same_edge <- mesh_loc[ids, ]
    mesh_same_edge_len <- mesh_edge_len[ids]
    # order the mesh_same_edge locations:
    ordering <- order(mesh_same_edge[, 2])
    mesh_same_edge_o <- mesh_same_edge[ordering, ]
    mesh_same_edge_len_o <- mesh_same_edge_len[ordering]
    # find the mesh edge index where the point is in
    index_on_edge <- which.max((mesh_same_edge_o[, 2] - as.numeric(coord[i, 2])) >= 0)

    mesh_index <- which(ids)[(ordering[index_on_edge])]
    mesh_where <- as.numeric((mesh_same_edge_o[index_on_edge, 2] - as.numeric(coord[i, 2])) / mesh_same_edge_len_o[index_on_edge]) # normalized
    new_coord[i, ] <- c(mesh_index, mesh_where)
  }
  fm_bary_MGM(graph, new_coord)
}

#' @title Convert fm_bary_MGM coordinates to a fm_bary_MGG coordinates
#' @description
#' Create a `fm_bary_MGG` object from MGM coordinates.
#'
#' @param graph metric_graph that the location should be mapped to.
#' @param coord MGM coordinates
#' @author Karina Lilleborge \email{karina.lilleborge@@gmail.com}
#' @returns An `fm_bary_MGG` object
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
#'   graph$build_mesh(h=0.01)
#'   mgg <- mesh_to_graph_coord(
#'     graph,
#'     matrix(c(5, 1),nrow=1)
#'   )
#'   mgg
#' }
#'
mesh_to_graph_coord <- function(graph,
                                coord) {
  mesh_loc <- graph$mesh$VtE
  new_coord <- matrix(nrow = nrow(coord), ncol = ncol(coord))
  for (i in seq_len(NROW(coord))) {
    # mesh vertex indices for neighboring loc
    mesh_vertices <- graph$mesh$E[coord[i, 1], ]
    # which graph coordinates for the neighboring vertices
    graph_edge_r <- mesh_loc[mesh_vertices[2], ]
    graph_edge_l <- mesh_loc[mesh_vertices[1], ]
    # as long as they are on the same edge
    if (graph_edge_l[1] == graph_edge_r[1]) {
      new_coord[i, ] <- c(graph_edge_l[1], (1 - coord[i, 2]) * graph_edge_l[2] + coord[i, 2] * graph_edge_r[2])
    } else {

    }
  }
  fm_bary_MGG(graph, new_coord)
}


#' @title Make a MGM object
#' @description
#' Create a `fm_bary_MGM` object
#'
#' @param graph metric_graph that the location should be mapped to.
#' @param loc MGM coordinates
#' @author Karina Lilleborge \email{karina.lilleborge@@gmail.com}
#' @returns An `fm_bary_MGG` object
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
#'   graph$build_mesh(h=0.01)
#'   m <- fm_bary_MGM(graph,
#'     matrix(c(1, 1),nrow=1)
#'   )
#'   class(m) # "fm_bary_MGM", "fm_bary", "tbl_df", "tbl", "data.frame"
#' }
#'
fm_bary_MGM <- function(graph,
                        loc) {
  res <- tibble::tibble(
    index = as.integer(loc[, 1]),
    where = as.numeric(loc[, 2])
  )
  coord <-
    structure(
      res,
      class = c("fm_bary_MGM", "fm_bary", "tbl_df", "tbl", "data.frame")
    )
  return(coord)
}


#' @title Make a graph_coord object from normalized PtE
#' @description
#' Create a `graph_coord` object from normalized PtE.
#'
#' @param graph metric_graph that the location should be mapped to.
#' @param loc PtE format (normalized)
#' @author Karina Lilleborge \email{karina.lilleborge@@gmail.com}
#' @returns An `fm_bary_MGM` object
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
#'   m <- fm_bary_MGG(graph, matrix(c(1, 0.5),nrow=1))
#'   class(m) # "fm_bary_MGG", "fm_bary", "tbl_df", "tbl", "data.frame"
#' }
#'
fm_bary_MGG <- function(graph,
                        loc) {
  res <- tibble::tibble(
    index = as.integer(loc[, 1]),
    where = as.numeric(loc[, 2])
  )
  coord <-
    structure(
      res,
      class = c("fm_bary_MGG", "fm_bary", "tbl_df", "tbl", "data.frame")
    )
  return(coord)
}


#' @title Make an inter edge interval on graph object
#' @description
#' Create an `MGG_interval` object.
#'
#' @param graph `metric_graph` that the interval should be mapped to.
#' @param start_MGG Start location for inter edge interval
#' @param end_MGG End location for inter edge interval
#' @author Karina Lilleborge \email{karina.lilleborge@@gmail.com}
#' @returns An `MGG_interval` object
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
#'   int <- MGG_interval(
#'     graph,
#'     matrix(c(1, 0.8),nrow=1),
#'     matrix(c(1, 0.5),nrow=1)
#'   )
#'   int
#' }
#'
MGG_interval <- function(graph,
                         start_MGG,
                         end_MGG) {
  inter_edge_interval <- structure(
    tibble::tibble(
      index = as.integer(start_MGG[, 1]),
      where = matrix(c(start_MGG[, 2], end_MGG[, 2]), ncol = 2)
    ),
    class = c("MGG_interval", "tbl_df", "tbl", "data.frame")
  )
  return(inter_edge_interval)
}

#' @title Make an interval on graph object
#' @description
#' Create a `path_MGG` object from known start, end and visiting edges.
#'
#' @param graph metric_graph that the interval should be mapped to.
#' @param start_MGG MGG coordinates for start
#' @param edges Ordered list of edge indices related to MGG
#' @param end_MGG MGG coordinates for end
#' @author Karina Lilleborge \email{karina.lilleborge@@gmail.com}
#' @returns A `path_MGG` object
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
#'   path <- path_MGG(
#'     graph,
#'     start_MGG = matrix(c(1, 0.5), nrow = 1),
#'     edges = c(2),
#'     end_MGG = matrix(c(3, 0.6), nrow = 1)
#'   )
#'   # a tibble with three interedge intervals
#'   # 1  0.5  0
#'   # 2  0    1
#'   # 3  0    0.6
#'   path
#' }
#'
path_MGG <- function(graph,
                     start_MGG,
                     edges,
                     end_MGG) {
  # check the graph does have circles
  if (length(edges) > 0) {
    # check direction from start to edges[1]
    v1 <- graph$E[as.integer(start_MGG[1L, 1L]), ]
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
    # index, start and end (MGG_interval)
    inter_edge_intervals <- matrix(nrow = length(edges) + 2, ncol = 3)
    inter_edge_intervals[1, ] <- c(
      as.integer(start_MGG[1, 1L]),
      as.numeric(start_MGG[1, 2L]),
      as.numeric(end_vertex)
    )
    # start and end must be determined
    # check direction for edges
    for (i in seq_len(length(edges))) {
      end_vertex <- c(0, 1)[!(v2 %in% v1[end_vertex + 1])]
      if (end_vertex == 0) start_vertex <- 1
      if (end_vertex == 1) start_vertex <- 0
      inter_edge_intervals[i + 1, ] <- c(
        as.integer(edges[i]),
        as.numeric(start_vertex),
        as.numeric(end_vertex)
      )
      v1 <- graph$E[as.integer(edges[i]), ]
      v2 <- graph$E[as.integer(edges[i + 1]), ]
    }
    v2 <- graph$E[as.integer(end_MGG[1L, 1L]), ]
    # check direction
    start_vertex <- c(0:1)[(v2 %in% v1[end_vertex + 1])]
    inter_edge_intervals[length(edges) + 2, ] <- c(
      as.integer(end_MGG[1L, 1L]),
      as.numeric(start_vertex),
      as.numeric(end_MGG[1L, 2L])
    )
  } else {
    # there are no whole edges visited
    v1 <- graph$E[as.integer(start_MGG[1L, 1L]), ]
    v2 <- graph$E[as.integer(end_MGG[1L, 1L]), ]
    if (sum(v2 %in% v1[1]) > 0) {
      # if the (e,0) vertex is in v2
      end_vertex <- 0
    }
    if (sum(v2 %in% v1[2]) > 0) {
      # if the (e,1) vertex is in v2
      end_vertex <- 1
    }
    # make storage for the inter edge intervals for each of the edge
    # index, start and end (MGG_interval)
    inter_edge_intervals <- matrix(nrow = 2, ncol = 3)
    inter_edge_intervals[1, ] <- c(
      as.integer(start_MGG[, 1L]),
      as.numeric(start_MGG[, 2L]),
      as.numeric(end_vertex)
    )

    v2 <- graph$E[as.integer(end_MGG[1L, 1L]), ]
    # check direction
    start_vertex <- c(0:1)[(v2 %in% v1[end_vertex + 1L])]
    inter_edge_intervals[2, ] <- c(
      as.integer(end_MGG[1L, 1L]),
      as.numeric(start_vertex),
      as.numeric(end_MGG[1L, 2L])
    )
  }
  path <- structure(
    tibble::tibble(
      index = as.integer(inter_edge_intervals[, 1L]),
      where = inter_edge_intervals[, -1L]
    ),
    class = c("path_MGG", "tbl_df", "tbl", "data.frame")
  )
  return(path)
}

#' @title Make an interval on graph object from sf object (NOT FINISHED)
#' @description
#' Create an `graph_interval` object from `sf::st_geometry` (`LINESTRING`)
#'
#' @param graph metric_graph that the interval should be mapped to.
#' @param geom_path `sf::st_geometry` (`LINESTRING`) on a graph
#' @author Karina Lilleborge \email{karina.lilleborge@@gmail.com}
#' @returns A `path_MGG` object
#' @export
#' @family object creation and conversion
#' @examples
#' if (requireNamespace("MetricGraph") && requireNamespace("lwgeom") &&
#'   requireNamespace("sf")) {
#'   edge1 <- rbind(c(0, 0), c(1, 0))
#'   edge2 <- rbind(c(0, 0), c(0, 1))
#'   edge3 <- rbind(c(0, 1), c(-1, 1))
#'   theta <- seq(from = pi, to = 3 * pi / 2, length.out = 20)
#'   edge4 <- cbind(sin(theta), 1 + cos(theta))
#'   edges <- list(edge1, edge2, edge3, edge4)
#'   graph <- MetricGraph::metric_graph$new(edges = edges)
#'   geom_path <- sf::st_sfc(sf::st_linestring(matrix(
#'     c(-1., 0., 0., 1., 1., 0.2),
#'     nrow = 3
#'   )))
#'   path <- geom_path_to_path_MGG(graph, geom_path)
#'   path
#' }
#'
geom_path_to_path_MGG <- function(graph, geom_path) {
  # get the start coordinates
  start_XY <- sf::st_coordinates(lwgeom::st_startpoint((geom_path)))
  start_MGG <- graph$coordinates(XY = start_XY)
  # get the end coordinates
  end_XY <- sf::st_coordinates(lwgeom::st_endpoint((geom_path)))
  end_MGG <- graph$coordinates(XY = end_XY)
  # matrix with colnames X Y and L1 (if geom_path is multiple linestrings)
  internal_XY <- sf::st_coordinates(geom_path)
  # determine the set of edges that connect start_MGG and end_MGG
  list_edges <- list()
  i <- 1
  for (l in unique(internal_XY[, "L1"])) {
    # a line
    line <- internal_XY[internal_XY[, "L1"] == l, ]
    line_MGG <- graph$coordinates(XY = line[, c("X", "Y")])
    # unique edges (consecutive repeats are removed)
    edges <- rle(line_MGG[, 1L])$values
    # remove the first entry if it is the same as the start
    if (edges[1] == start_MGG[l, 1L]) {
      edges <- edges[-1]
    }
    # remove the last entry if it is the same as the end
    if (edges[length(edges)] == end_MGG[l, 1L]) {
      edges <- edges[-(length(edges))]
    }
    # what can happen, is that edges are included now because an endpoint
    # of a certain edge is in the path, but not the full edge...
    # need to deal with this..

    list_edges[[i]] <- edges
    i <- i + 1L
  }

  # should call path MGG
  paths <- list()
  for (i in seq_len(length(list_edges))) {
    path_i <- path_MGG(
      graph = graph,
      start_MGG = fm_bary_MGG(
        graph,
        matrix(start_MGG[i, ], nrow = 1)
      ),
      edges = list_edges[[i]],
      end_MGG = fm_bary_MGG(
        graph,
        matrix(end_MGG[i, ], nrow = 1)
      )
    )
    paths[[i]] <- path_i
  }

  return(paths)
}
