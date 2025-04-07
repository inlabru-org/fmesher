test_that("Flat CDT works", {
  max.edge0 <- 0.25
  min.angle0 <- 21
  mesh <- fm_rcdt_2d_inla(cbind(0, 0),
    extend = list(offset = 1, n = 16),
    refine = list(max.edge = max.edge0, min.angle = min.angle0)
  )

  #  expect_equal(fm_dof(mesh), 135) # 138 on M1?

  expect_equal(fm_diameter(mesh), 2.039182, tolerance = midtol)

  expect_equal(fm_manifold(mesh), "R2")
  expect_equal(fm_manifold_type(mesh), "R")
  expect_equal(fm_manifold_dim(mesh), 2)
  expect_equal(fm_manifold(mesh, "R"), TRUE)
  expect_equal(fm_manifold(mesh, "2"), TRUE)
  expect_equal(fm_manifold(mesh, "R2"), TRUE)
  expect_equal(fm_manifold(mesh, "S"), FALSE)
  expect_equal(fm_manifold(mesh, "1"), FALSE)
  expect_equal(fm_manifold(mesh, "S1"), FALSE)
  # Check issue #16, where it ignored all but the first type option:
  expect_equal(fm_manifold(mesh, c("R", "S")), TRUE)
  expect_equal(fm_manifold(mesh, c("S", "R")), TRUE)

  edges <- list(
    mesh$loc[mesh$graph$tv[, 2], ] - mesh$loc[mesh$graph$tv[, 1], ],
    mesh$loc[mesh$graph$tv[, 3], ] - mesh$loc[mesh$graph$tv[, 2], ],
    mesh$loc[mesh$graph$tv[, 1], ] - mesh$loc[mesh$graph$tv[, 3], ]
  )
  min.angle <-
    180 / pi * min(acos(pmin(1, pmax(
      -1, c(
        -rowSums(edges[[1]] * edges[[2]]) /
          (rowSums(edges[[1]]^2)^0.5 * rowSums(edges[[2]]^2)^0.5),
        -rowSums(edges[[2]] * edges[[3]]) /
          (rowSums(edges[[2]]^2)^0.5 * rowSums(edges[[3]]^2)^0.5),
        -rowSums(edges[[3]] * edges[[1]]) /
          (rowSums(edges[[3]]^2)^0.5 * rowSums(edges[[1]]^2)^0.5)
      )
    ))))
  max.edge <- max(rowSums(do.call(rbind, edges)^2)^0.5)

  expect_lte(max.edge, max.edge0 + lowtol)
  expect_gte(min.angle, min.angle0 - lowtol)
})

test_that("Spherical CDT works", {
  max.edge0 <- 0.5
  mesh <- fm_rcdt_2d_inla(globe = 1, refine = list(max.edge = max.edge0))
  #  expect_equal(fm_dof(mesh), 108) # 106 or 107 on M1?

  expect_equal(fm_diameter(mesh), pi)

  expect_equal(fm_manifold(mesh), "S2")
  expect_equal(fm_manifold_type(mesh), "S")
  expect_equal(fm_manifold_dim(mesh), 2)
  expect_equal(fm_manifold(mesh, "S"), TRUE)
  expect_equal(fm_manifold(mesh, "2"), TRUE)
  expect_equal(fm_manifold(mesh, "S2"), TRUE)
  expect_equal(fm_manifold(mesh, "R"), FALSE)
  expect_equal(fm_manifold(mesh, "1"), FALSE)
  expect_equal(fm_manifold(mesh, "R1"), FALSE)

  edges <- list(
    mesh$loc[mesh$graph$tv[, 2], ] - mesh$loc[mesh$graph$tv[, 1], ],
    mesh$loc[mesh$graph$tv[, 3], ] - mesh$loc[mesh$graph$tv[, 2], ],
    mesh$loc[mesh$graph$tv[, 1], ] - mesh$loc[mesh$graph$tv[, 3], ]
  )
  sums <- list(
    mesh$loc[mesh$graph$tv[, 2], ] + mesh$loc[mesh$graph$tv[, 1], ],
    mesh$loc[mesh$graph$tv[, 3], ] + mesh$loc[mesh$graph$tv[, 2], ],
    mesh$loc[mesh$graph$tv[, 1], ] + mesh$loc[mesh$graph$tv[, 3], ]
  )
  euc_len <- rowSums(do.call(rbind, edges)^2)^0.5
  sum_len <- rowSums(do.call(rbind, sums)^2)^0.5

  len <- 2 * atan2(euc_len, sum_len)
  max.edge <- max(len)

  expect_lte(max.edge, max.edge0 + lowtol)
})


test_that("fm_lattice_2d ordering", {
  latt1 <- fm_lattice_2d(1:4, 1:3)
  latt2 <- fm_lattice_2d(rev(1:4), 1:3)
  expect_equal(latt1$x, latt2$x)
  expect_equal(latt1$y, latt2$y)
  expect_equal(latt1$loc, latt2$loc)
})


test_that("interior should be single object", {
  # create boundary polygon
  bnd <- sf::st_sfc(sf::st_polygon(list(as.matrix(
    data.frame(x = c(0, 0, 10, 10, 0), y = c(0, 10, 10, 0, 0))
  ))))
  # create interior polygon
  interior <- sf::st_sfc(sf::st_polygon(list(as.matrix(
    data.frame(x = c(3, 3, 7, 7, 3), y = c(0, 5, 5, 0, 0))
  ))))

  # simple mesh
  expect_error(fm_mesh_2d_inla(boundary = bnd, max.edge = 1), NA)
  # with interior
  expect_error(
    fm_mesh_2d_inla(
      boundary = bnd,
      interior = interior,
      max.edge = 1
    ),
    NA
  )
  # extended boundary
  expect_error(fm_mesh_2d_inla(boundary = bnd, max.edge = c(1, 2)), NA)
  # or extended by running non-convex for the boundary
  expect_error(fm_mesh_2d_inla(
    boundary = list(bnd, fm_nonconvex_hull(bnd)),
    max.edge = 1
  ), NA)

  # combining extended and interior segment should not fail
  expect_error(
    fm_mesh_2d_inla(boundary = bnd, interior = interior, max.edge = c(1, 2)),
    NA
  )
  expect_error(
    fm_mesh_2d_inla(
      boundary = list(bnd, fm_nonconvex_hull(bnd)),
      interior = interior, max.edge = 1
    ),
    NA
  )

  # List interior supported from 0.2.0.9016
  expect_error(
    fm_mesh_2d_inla(
      boundary = bnd,
      interior = list(interior),
      max.edge = c(1, 2)
    ),
    NA
  )
  expect_error(
    fm_mesh_2d_inla(
      boundary = list(bnd, fm_nonconvex_hull(bnd)),
      interior = list(interior), max.edge = 1
    ),
    NA
  )
})
