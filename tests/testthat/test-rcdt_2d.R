test_that("Spherical CDT works", {
  mesh <- fm_rcdt_2d_inla(globe = 1, refine = list(max.edge = 0.25))

  expect_equal(fm_dof(mesh), 484)

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
