mesh <- fm_rcdt_2d_inla(globe = 1)
plot(mesh)

fem1 <- fm_fem(mesh, order = 2)
fem2 <- fm_fem(mesh, order = 2, aniso = list(
  gamma = rep(1, mesh$n),
  v = matrix(0, mesh$n, 3)
))
fem3 <- fm_fem_aniso(mesh, aniso = list(
  kappa = rep(1, mesh$n),
  v = matrix(0, mesh$n, 3)
))
# Test for C0
test_that("C0 matrices are equal", {
  expect_equal(fem1$c0, fem2$c0, label = "fem1$c0 is not equal to fem2$c0")
  expect_equal(fem1$c0, fem3$c0, label = "fem1$c0 is not equal to fem3$c0")
  expect_equal(fem2$c0, fem3$c0, label = "fem2$c0 is not equal to fem3$c0")
})

# Test for C1
test_that("C1 matrices are equal", {
  expect_equal(fem1$c1, fem2$c1, label = "fem1$c1 is not equal to fem2$c1")
  expect_equal(fem1$c1, fem3$c1, label = "fem1$c1 is not equal to fem3$c1")
  expect_equal(fem2$c1, fem3$c1, label = "fem2$c1 is not equal to fem3$c1")
})

# Test for G1
test_that("G1 matrices are equal", {
  expect_equal(fem1$g1, fem2$g1, label = "fem1$g1 is not equal to fem2$g1")
  expect_equal(fem1$g1, fem3$g1, label = "fem1$g1 is not equal to fem3$g1")
  expect_equal(fem2$g1, fem3$g1, label = "fem2$g1 is not equal to fem3$g1")
})

# Test for G2
test_that("G2 matrix is equal", {
  expect_equal(fem1$g2, fem2$g2, label = "fem1$g2 is not equal to fem2$g2")
  expect_equal(fem1$g2, fem3$g2, label = "fem1$g2 is not equal to fem3$g2")
  expect_equal(fem2$g2, fem3$g2, label = "fem2$g2 is not equal to fem3$g2")
})
