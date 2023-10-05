test_that("GMRF methods, 1d", {
  mesh <- fm_mesh_1d(1:5, boundary = "neumann", degree = 2)

  Q <- fm_matern_precision(mesh, alpha = 2, rho = 1, sigma = 1)
  expect_equal(dim(Q), c(1, 1) * fm_dof(mesh))

  set.seed(12345L)
  samp <- fm_sample(n = 2, Q = Q)
  expect_equal(dim(samp), c(fm_dof(mesh), 2))

  set.seed(12345L)
  samp_constr <- fm_sample(n = 2, Q = Q, constr = list(A = matrix(1, 1, fm_dof(mesh)), e = 2))
  expect_equal(dim(samp), c(fm_dof(mesh), 2))
  expect_equal(colSums(samp_constr), c(2, 2))
})

test_that("GMRF methods, 2d", {
  mesh <- fm_mesh_2d_inla(
    boundary = fm_extensions(cbind(2, 1), 3),
    max.edge = 0.5
  )

  Q <- fm_matern_precision(mesh, alpha = 2, rho = 1, sigma = 1)
  expect_equal(dim(Q), c(1, 1) * fm_dof(mesh))

  set.seed(12345L)
  samp <- fm_sample(n = 2, Q = Q)
  expect_equal(dim(samp), c(fm_dof(mesh), 2))

  set.seed(12345L)
  samp_constr <- fm_sample(n = 2, Q = Q, constr = list(A = matrix(1, 1, fm_dof(mesh)), e = 2))
  expect_equal(dim(samp), c(fm_dof(mesh), 2))
  expect_equal(colSums(samp_constr), c(2, 2))
})

test_that("GMRF methods, anisotropic", {
  mesh <- fm_mesh_2d_inla(
    boundary = fm_extensions(cbind(2, 1), 3),
    max.edge = 0.5
  )

  Q <- fm_aniso_precision(mesh, list(rep(1, mesh$n), matrix(1, mesh$n, 2)))
  expect_equal(dim(Q), c(1, 1) * fm_dof(mesh))

  set.seed(12345L)
  samp <- fm_sample(n = 2, Q = Q)
  expect_equal(dim(samp), c(fm_dof(mesh), 2))

  set.seed(12345L)
  samp_constr <- fm_sample(n = 2, Q = Q, constr = list(A = matrix(1, 1, fm_dof(mesh)), e = 2))
  expect_equal(dim(samp), c(fm_dof(mesh), 2))
  expect_equal(colSums(samp_constr), c(2, 2))
})
