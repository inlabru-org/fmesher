# Testing of the metric graph extensions to inlabru/fmesher
# devtools::test(filter="metric_graph")
# devtools::document()
local_bru_test_graph <- function() {
  edge1 <- rbind(c(0, 0), c(1, 0))
  edge2 <- rbind(c(0, 0), c(0, 1))
  edge3 <- rbind(c(0, 1), c(-1, 1))
  theta <- seq(from = pi, to = 3 * pi / 2, length.out = 20)
  edge4 <- cbind(sin(theta), 1 + cos(theta))
  edge5 <- rbind(c(0, 1), c(1, 1))
  edge6 <- rbind(c(1, 1), c(1, 0))
  edge7 <- rbind(c(1, 1), c(2, 1))
  edges <- list(edge1, edge2, edge3, edge4, edge5, edge6, edge7)
  graph <- MetricGraph::metric_graph$new(edges = edges)
  return(graph)
}

test_that("MGG bary", {
  skip_if_not_installed("MetricGraph")
  graph0 <- local_bru_test_graph()
  graph0$build_mesh(h = 0.005)
  # Euclidean
  locs <- rbind(c(0, 0.6), c(1, 0.20))
  b <-
    fm_bary(
      mesh = fm_as_MG(graph0, MGG = TRUE),
      loc = locs
    )
  expect_equal(
    c(
      b$index[1, drop = FALSE],
      b$index[2, drop = FALSE],
      b$where[1, 2, drop = FALSE],
      b$where[2, 2, drop = FALSE]
    ),
    c(2, 6, 0.6, 0.8)
  )
  # MGG
  locs <- as_MGG(rbind(c(2, 0.6), c(6, 0.8)))
  b <-
    fm_bary(
      mesh = fm_as_MG(graph0, MGG = TRUE),
      loc = locs
    )
  expect_equal(
    c(
      b$index[1, drop = FALSE],
      b$index[2, drop = FALSE],
      b$where[1, 2, drop = FALSE],
      b$where[2, 2, drop = FALSE]
    ),
    c(2, 6, 0.6, 0.8)
  )
  # MGM
  locs <- as_MGM(rbind(c(260, 1), c(400, 0.8)))
  b <-
    fm_bary(
      mesh = fm_as_MG(graph0, MGG = FALSE),
      loc = locs
    )
  expect_equal(
    c(
      b$index[1, drop = FALSE],
      b$index[2, drop = FALSE],
      b$where[1, 2, drop = FALSE],
      b$where[2, 2, drop = FALSE]
    ),
    c(260, 400, 1, 0.8)
  )

  # Euclidean coordinates
  locs <- sf::st_multipoint(rbind(c(0, 1), c(1, 0.8)))
  b <-
    fm_bary(
      mesh = fm_as_MG(graph0, MGG = TRUE),
      loc = locs
    )
  expect_equal(
    c(
      b$index[1, drop = FALSE],
      b$index[2, drop = FALSE],
      b$where[1, 2, drop = FALSE],
      b$where[2, 2, drop = FALSE]
    ),
    c(2, 6, 1, 0.2)
  )
})

test_that("MGG to MGM", {
  skip_if_not_installed("MetricGraph")
  graph0 <- local_bru_test_graph()
  graph0$build_mesh(h = 0.005)
  locs <- as_MGG(rbind(
    c(1, 0.6),
    c(3, 0.20)
  ))
  expect_error(MGM_to_MGG(
    graph = graph0,
    coord = locs
  ))
  mgm <-
    MGG_to_MGM(
      graph = graph0,
      coord = locs
    )
  mgm2 <- fm_bary(
    mesh = fm_as_MG(graph0, MGG = FALSE),
    loc = locs
  )
  expect_equal(
    c(
      mgm$index[1, drop = FALSE],
      mgm$index[2, drop = FALSE]
    ),
    c(120, 440)
  )
  expect_equal(
    c(
      mgm$where[1, 2, drop = FALSE],
      mgm$where[2, 2, drop = FALSE]
    ),
    c(1, 1)
  )
  expect_equal(
    c(
      mgm$index,
      mgm$where
    ),
    c(
      mgm2$index,
      mgm2$where
    )
  )
})


test_that("MGM to MGG", {
  skip_if_not_installed("MetricGraph")
  graph0 <- local_bru_test_graph()
  graph0$build_mesh(h = 0.005)
  locs <- rbind(c(300, 0.5), c(1250, 1))
  expect_error(MGM_to_MGG(
    graph = graph0,
    coord = locs
  ))
  mgg <-
    MGM_to_MGG(
      graph = graph0,
      coord = as_MGM(locs)
    )

  expect_equal(
    c(
      mgg$index[1, drop = FALSE],
      mgg$index[2, drop = FALSE]
    ),
    c(
      2,
      6
    )
  )
  expect_equal(
    c(
      mgg$where[1, 2, drop = FALSE],
      mgg$where[2, 2, drop = FALSE]
    ),
    c(
      0.4975,
      0.6750
    )
  )
})

test_that("bary MGG to MGG", {
  skip_if_not_installed("MetricGraph")
  graph0 <- local_bru_test_graph()
  locs <- as_MGG(cbind(c(2, 5), c(0.8, 0.2)))
  b <-
    fm_bary(
      mesh = fm_as_MG(graph0, MGG = TRUE),
      loc = locs
    )
  expect_equal(
    c(
      b$where[1, 2, drop = FALSE],
      b$where[2, 2, drop = FALSE]
    ),
    c(
      0.8,
      0.2
    )
  )
})

test_that("bary MGM to MGM", {
  skip_if_not_installed("MetricGraph")
  graph0 <- local_bru_test_graph()
  graph0$build_mesh(h = 0.005)
  locs <- as_MGM(cbind(c(300, 1250), c(0.5, 1.0)))
  b <-
    fm_bary(
      mesh = fm_as_MG(graph0, MGG = FALSE),
      loc = locs
    )
  expect_equal(
    c(
      b$index[1, drop = FALSE],
      b$index[2, drop = FALSE]
    ),
    c(
      300,
      1250
    )
  )
})


test_that("path construction", {
  skip_if_not_installed("MetricGraph")
  graph0 <- local_bru_test_graph()

  # same edge interval
  start <- cbind(2, 0.5)
  edges <- c()
  end <- cbind(2, 0.8)
  p <-
    simple_path_MGG(
      graph = graph0,
      start_MGG = start,
      edges = edges,
      end_MGG = end
    )

  expect_equal(
    c(
      p$start$index,
      p$start$where[, 2, drop = FALSE],
      p$end$where[, 2, drop = FALSE]
    ),
    c(
      2,
      0.5,
      0.8
    )
  )

  # neighboring edges
  start <- cbind(2, 0.5)
  edges <- c()
  end <- cbind(4, 0.8)
  p <-
    simple_path_MGG(
      graph = graph0,
      start_MGG = start,
      edges = edges,
      end_MGG = end
    )

  expect_equal(
    c(
      p$start$index,
      p$start$where[, 2, drop = FALSE],
      p$end$where[, 2, drop = FALSE]
    ),
    c(
      2, 4,
      0.5, 0,
      0, 0.8
    )
  )

  # with edges
  start <- cbind(2, 0.5)
  edges <- c(1, 6, 5)
  end <- cbind(3, 0.8)
  p <-
    simple_path_MGG(
      graph = graph0,
      start_MGG = start,
      edges = edges,
      end_MGG = end
    )

  expect_equal(
    c(
      p$start$index,
      p$start$where[, 2, drop = FALSE],
      p$end$where[, 2, drop = FALSE]
    ),
    c(
      2, 1, 6, 5, 3,
      0.5, 0, 1, 1, 0,
      0, 1, 0, 0, 0.8
    )
  )
})


# making single path, multiple paths etc

# fm_int test
test_that("integration one path", {
  skip_if_not_installed("MetricGraph")
  edge1 <- rbind(c(0, 0), c(1, 0))
  edge2 <- rbind(c(0, 0), c(0, 1))
  graph0 <- MetricGraph::metric_graph$new(edges = list(
    edge1,
    edge2
  ))
  start1 <- cbind(1, 0.5)
  end1 <- cbind(1, 0.8)
  p1 <-
    simple_path_MGG(
      graph = graph0,
      start_MGG = start1,
      edges = c(),
      end_MGG = end1
    )
  test_sampler <- tibble::tibble(x = list(p1), weight = 1)

  # there is no mesh in the graph yet, test that fm_int checks for the mesh
  expect_error(
    fm_int(graph0, samplers = test_sampler),
    "There is no mesh"
  )

  # build mesh and check output is correct
  graph0$build_mesh(h = 0.005)
  # expect no error with NA
  expect_error(
    fm_int(graph0, samplers = test_sampler),
    NA
  )
  ips <- fm_int(graph0, samplers = test_sampler)
  expect_equal(
    c(
      unique(ips$x[["index"]])
    ),
    1
  )
  expect_equal(
    c(
      sum(ips$weight)
    ),
    0.3
  )
})

test_that("integration two paths", {
  skip_if_not_installed("MetricGraph")
  graph0 <- local_bru_test_graph()
  start1 <- cbind(2, 0.5)
  edges1 <- c(1, 6, 5)
  end1 <- cbind(3, 0.8)
  p1 <-
    simple_path_MGG(
      graph = graph0,
      start_MGG = start1,
      edges = edges1,
      end_MGG = end1
    )
  start2 <- cbind(7, 0.2)
  edges2 <- c(5, 3, 4, 1)
  end2 <- cbind(6, 0.8)
  p2 <-
    simple_path_MGG(
      graph = graph0,
      start_MGG = start2,
      edges = edges2,
      end_MGG = end2
    )

  test_sampler <- tibble::tibble(x = list(p1, p2), weight = c(2, 1))

  # there is no mesh in the graph yet, test that fm_int checks for the mesh
  expect_error(
    fm_int(graph0, samplers = test_sampler),
    "There is no mesh"
  )


  # build mesh and check output is correct
  graph0$build_mesh(h = 0.005)
  # expect no error with NA
  expect_error(
    fm_int(graph0, samplers = test_sampler),
    NA
  )
  ips <- fm_int(graph0, samplers = test_sampler)
  expect_equal(
    c(
      unique(ips$x[["index"]])
    ),
    unique(c(2, 1, 6, 5, 3, 7, 5, 3, 4, 1, 6))
  )
})

test_that("fm_basis paths", {
  skip_if_not_installed("MetricGraph")
  graph0 <- local_bru_test_graph()
  start1 <- cbind(2, 0.5)
  edges1 <- c(1, 6, 5)
  end1 <- cbind(3, 0.8)
  p1 <-
    simple_path_MGG(
      graph = graph0,
      start_MGG = start1,
      edges = edges1,
      end_MGG = end1
    )
  start2 <- cbind(7, 0.2)
  edges2 <- c(5, 3, 4, 1)
  end2 <- cbind(6, 0.8)
  p2 <-
    simple_path_MGG(
      graph = graph0,
      start_MGG = start2,
      edges = edges2,
      end_MGG = end2
    )
  test_sampler <- tibble::tibble(x = list(p1, p2), weight = c(1, 1))
  graph0$build_mesh(h = 0.005)
  ips <- fm_int(fm_as_MG(graph0, MGG = FALSE), samplers = test_sampler)
  basis <- fm_basis(x = fm_as_MG(graph0, MGG = FALSE),
                    loc = ips$x,
                    weights = ips$weight)
  n <- NROW(ips)
  MGM_locs <- as_MGM(ips$x, graph = graph0)
  true_A <- Matrix::sparseMatrix(
    i = c(seq_len(n), seq_len(n)),
    j = c(graph0$mesh$E[MGM_locs$index, 1], graph0$mesh$E[MGM_locs$index, 2]),
    x = c(ips$weight * MGM_locs$where[, 1], ips$weight * MGM_locs$where[, 2]),
    dims = c(n, NROW(graph0$mesh$V))
  )
  expect_equal(
    basis,
    true_A
  )
  expect_equal(
    sum(apply(basis[ips$.block == 1, ], 2, sum)),
    4.3
  )
  expect_equal(
    sum(apply(basis[ips$.block == 2, ], 2, sum)),
    3.4 + pi / 2,
    tolerance = midtol
  )
})


# detect a error with the error message
# expect_error(code, "message to match")


test_that("ibm values", {
  skip_if_not_installed("MetricGraph")
  graph0 <- local_bru_test_graph()
  graph0$build_mesh(h = 0.005)
  mapper <- bru_mapper_metric_graph(graph0, n_rep = 2)
  values <- inlabru::ibm_values(mapper)
  expect_equal(
    values,
    seq_len(2 * NROW(graph0$mesh$V))
  )
})

test_that("sf to MGG", {
  skip_if_not_installed("MetricGraph")
  skip_if_not_installed("sf")
  skip_if_not_installed("lwgeom")
  graph0 <- local_bru_test_graph()
  graph0$build_mesh(h = 0.005)
  line1 <- sf::st_linestring(cbind(c(0, 0, 1), c(0.5, 0, 0)))
  line1_g <- sf::st_geometry(line1)
  path_MGG1 <-
    geom_path_to_path_MGG(
      graph = graph0,
      geom_path = line1_g
    )
  expect_equal(
    path_MGG1$paths$start$index,
    c(2, 1)
  )
  expect_equal(
    cbind(path_MGG1$paths$start$where[, 2], path_MGG1$paths$end$where[, 2]),
    cbind(c(0.5, 0), c(0, 1))
  )
  line2 <- sf::st_linestring(cbind(c(-1, 0, 1), c(1, 1, 1)))
  lines <- sf::st_sfc(list(line1, line2))
  lines <- sf::st_geometry(lines)
  path_MGGs <-
    geom_path_to_path_MGG(
      graph = graph0,
      geom_path = lines
    )

  expect_equal(
    path_MGGs$paths$start$index,
    c(2, 1, 3, 5)
  )
  expect_equal(
    cbind(path_MGGs$paths$start$where[, 2], path_MGGs$paths$end$where[, 2]),
    cbind(c(0.5, 0, 1, 0), c(0, 1, 0, 1))
  )
})
