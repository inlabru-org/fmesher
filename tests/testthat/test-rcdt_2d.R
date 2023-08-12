test_that("Spherical CDT works", {
  mesh <- fm_rcdt_2d_inla(globe = 1, refine = list(max.edge = 0.5))

  expect_equal(fm_dof(mesh), 108)

  expect_equal(fm_diameter(mesh), pi)
})
