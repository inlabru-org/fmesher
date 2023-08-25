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

  edges <- list(
    mesh$loc[mesh$graph$tv[, 2], ] - mesh$loc[mesh$graph$tv[, 1], ],
    mesh$loc[mesh$graph$tv[, 3], ] - mesh$loc[mesh$graph$tv[, 2], ],
    mesh$loc[mesh$graph$tv[, 1], ] - mesh$loc[mesh$graph$tv[, 3], ]
  )
  min.angle <-
    180 / pi * min(acos(pmin(1, pmax(
      -1, c(
        -rowSums(edges[[1]] * edges[[2]]) / rowSums(edges[[1]]^2)^0.5 / rowSums(edges[[2]]^2)^0.5,
        -rowSums(edges[[2]] * edges[[3]]) / rowSums(edges[[2]]^2)^0.5 / rowSums(edges[[3]]^2)^0.5,
        -rowSums(edges[[3]] * edges[[1]]) / rowSums(edges[[3]]^2)^0.5 / rowSums(edges[[1]]^2)^0.5
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
