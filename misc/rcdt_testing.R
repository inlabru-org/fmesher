rcdt_testing <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
  initial_sink <- sink.number()
  on.exit(if (sink.number() > initial_sink) sink(NULL))

  sink(file.path(path, "test_rcdt_R2_plain.txt"))
  meshR2 <- fm_rcdt_2d_inla(cbind(0, 0),
                            extend = list(offset = 1, n = 16))
  sink(NULL)

  sink(file.path(path, "test_rcdt_S2_plain.txt"))
  meshS2 <- fm_rcdt_2d_inla(globe = 1)
  sink(NULL)

  sink(file.path(path, "test_rcdt_R2_refine.txt"))
  meshR2_ref <- fm_rcdt_2d_inla(cbind(0, 0),
                                extend = list(offset = 1, n = 16),
                                refine = list(max.edge = 0.25))
  sink(NULL)

  sink(file.path(path, "test_rcdt_S2_refine.txt"))
  meshS2_ref <- fm_rcdt_2d_inla(globe = 1, refine = list(max.edge = 0.5))
  sink(NULL)

  invisible()
}

  #  expect_equal(fm_dof(meshR2_ref), 135) # 138 on M1?
  #  expect_equal(fm_dof(meshS2_ref), 108) # 106 or 107 on M1?
