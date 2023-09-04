test_that("Basic 2d FEM works", {
  mesh <- fm_rcdt_2d_inla(globe = 1)

  fem1 <- fm_fem(mesh, order = 2)
  names_fem <- c("b1", "c0", "c1", "g1", "g2", "k1", "k2", "ta", "va")
  expect_setequal(names(fem1), names_fem)

  expect_error(
    fm_fem(mesh, order = 2, aniso = list()),
    regexp = "must be NULL or a list of length 2"
  )

  fem2 <- fm_fem(mesh, order = 2, aniso = list(
    gamma = rep(1, mesh$n),
    v = matrix(1, mesh$n, 3)
  ))
  names_aniso <- union(names_fem, c("g1aniso", "g2aniso"))
  expect_setequal(names(fem2), names_aniso)

  fem3 <- fm_fem_aniso(mesh, aniso = list(
    kappa = rep(1, mesh$n),
    v = matrix(1, mesh$n, 3)
  ))
  names_aniso2 <- c("c0", "c1", "g1", "g2", "ta", "va")
  expect_setequal(names(fem3), names_aniso2)
})

test_that("Basic 1d FEM works", {
  names_fem0 <- c("c0", "c1", "g1", "g2")
  names_fem2 <- union(names_fem0, c("g01", "g02", "g12"))
  names_fem <- list(names_fem0, names_fem0, names_fem2)

  configs <-
    expand.grid(
      degree = 0:2,
      bnd1 = c("neumann", "dirichlet", "free", "cyclic"),
      bnd2 = c("neumann", "dirichlet", "free", "cyclic"),
      free.clamped1 = c(FALSE, TRUE),
      free.clamped2 = c(FALSE, TRUE)
    )
  nok <- ((configs$bnd1 != "free") & configs$free.clamped1) |
    ((configs$bnd2 != "free") & configs$free.clamped2) |
    ((configs$bnd1 == "cyclic") & (configs$bnd2 != "cyclic")) |
    ((configs$bnd1 != "cyclic") & (configs$bnd2 == "cyclic"))
  configs <- configs[!nok, ]

  for (k in seq_len(nrow(configs))) {
    mesh <- fm_mesh_1d(
      c(1, 2, 4, 8),
      degree = configs$degree[k],
      boundary = as.character(c(configs$bnd1[k], configs$bnd2[k])),
      free.clamped = c(configs$free.clamped1[k], configs$free.clamped2[k])
    )

    fem <- fm_fem(mesh)

    expect_setequal(names(fem), names_fem[[configs$degree[k] + 1]])
  }
})
