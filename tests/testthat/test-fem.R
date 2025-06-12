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


test_that("Basic fm_collect FEM works", {
  mesh1 <- fmexample$mesh
  fm_crs(mesh1) <- fm_crs("longlat_globe")
  # Two transformations that both have kilometres as units:
  mesh1 <- fm_transform(mesh1, fm_crs("mollweide_globe"))
  mesh2 <- fm_transform(mesh1, fm_crs("globe"))
  mesh3 <- fm_rcdt_2d(globe = 1, crs = fm_crs("globe"))

  expect_error(
    {
      fm_collect(list(mesh1, mesh2, mesh3))
    },
    "All function spaces in a collection need to be of the same manifold type"
  )
  expect_no_error({
    mesh <- fm_collect(list(mesh2, mesh3))
  })

  loc_sf <- fmexample$loc_sf
  sf::st_crs(loc_sf) <- fm_crs("longlat_globe")
  basis <- fm_basis(
    mesh,
    list(
      c(loc_sf, loc_sf),
      rep(1:2, each = NROW(loc_sf))
    )
  )

  fem <- fm_fem(mesh, order = 4)
  names_fem <- c("cc", "c0", "c1", "g1", "g2", "g3", "g4", "ta", "va")
  expect_setequal(names(fem), names_fem)
})
