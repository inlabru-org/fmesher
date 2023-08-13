test_that("Flat CDT works", {
  mesh <- fm_rcdt_2d_inla(cbind(0, 0),
    extend = list(offset = 1, n = 16),
    refine = list(max.edge = 0.25)
  )

  expect_equal(fm_dof(mesh), 135)

  expect_equal(fm_diameter(mesh), 2, tolerance = 0.1)

  expect_equal(fm_manifold(mesh), "R2")
  expect_equal(fm_manifold_type(mesh), "R")
  expect_equal(fm_manifold_dim(mesh), 2)
  expect_equal(fm_manifold(mesh, "R"), TRUE)
  expect_equal(fm_manifold(mesh, "2"), TRUE)
  expect_equal(fm_manifold(mesh, "R2"), TRUE)
  expect_equal(fm_manifold(mesh, "S"), FALSE)
  expect_equal(fm_manifold(mesh, "1"), FALSE)
  expect_equal(fm_manifold(mesh, "S1"), FALSE)
})

test_that("Spherical CDT works", {
  mesh <- fm_rcdt_2d_inla(globe = 1, refine = list(max.edge = 0.5))
  #  expect_equal(fm_dof(mesh), 108) # 107 on macos!

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
})
