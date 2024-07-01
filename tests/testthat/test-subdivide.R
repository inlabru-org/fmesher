test_that("fmesher_subdivide works", {
  mesh <- fm_rcdt_2d_inla(
    rbind(c(0, 0), c(1, 0), c(0, 1)),
    rbind(c(0, 1, 2)) + 1L
  ) # ,
  #                          refine=list(max.edge = 0.1))

  # 2
  # 0, 1
  sub0.tv <- rbind(c(0, 1, 2))
  sub0 <- fmesher_subdivide(
    mesh$loc,
    mesh$graph$tv - 1L,
    mesh$segm$bnd$idx - 1L,
    mesh$segm$int$idx - 1L,
    0L,
    list()
  )
  expect_equal(sub0$loc, mesh$loc)
  expect_equal(sub0$tv, sub0.tv)

  # 2
  # 5, 4
  # 0, 3, 1
  sub1.tv <- rbind(
    c(0, 3, 5),
    c(3, 1, 4),
    c(4, 5, 3),
    c(5, 4, 2)
  )
  sub1.loc <- rbind(
    mesh$loc,
    cbind(rbind(c(0.5, 0), c(0.5, 0.5), c(0, 0.5)), 0)
  )
  sub1 <- fmesher_subdivide(
    mesh$loc,
    mesh$graph$tv - 1L,
    mesh$segm$bnd$idx - 1L,
    mesh$segm$int$idx - 1L,
    1L,
    list()
  )
  expect_equal(sub1$loc, sub1.loc)
  expect_equal(sub1$tv, sub1.tv)

  # 2
  # 7, 6
  # 8, 9, 5
  # 0, 3, 4, 1
  sub2.tv <- rbind(
    c(0, 3, 8),
    c(3, 4, 9),
    c(4, 1, 5),
    c(9, 8, 3),
    c(5, 9, 4),
    c(8, 9, 7),
    c(9, 5, 6),
    c(6, 7, 9),
    c(7, 6, 2)
  )
  sub2.loc <- rbind(
    mesh$loc,
    cbind(rbind(
      c(1 / 3, 0), c(2 / 3, 0),
      c(2 / 3, 1 / 3), c(1 / 3, 2 / 3),
      c(0, 2 / 3), c(0, 1 / 3),
      c(1 / 3, 1 / 3)
    ), 0)
  )
  sub2 <- fmesher_subdivide(
    mesh$loc,
    mesh$graph$tv - 1L,
    mesh$segm$bnd$idx - 1L,
    mesh$segm$int$idx - 1L,
    2L,
    list()
  )
  expect_equal(sub2$loc, sub2.loc)
  expect_equal(sub2$tv, sub2.tv)

  m0 <- fm_rcdt_2d_inla(sub0$loc, sub0$tv + 1L)
  m1 <- fm_rcdt_2d_inla(sub1$loc, sub1$tv + 1L)
  m2 <- fm_rcdt_2d_inla(sub2$loc, sub2$tv + 1L)

  m0_ <- fm_subdivide(mesh, 0)
  m1_ <- fm_subdivide(mesh, 1)
  m2_ <- fm_subdivide(mesh, 2)

  expect_identical(m0, m0_)
  expect_identical(m1, m1_)
  expect_identical(m2, m2_)
})

test_that("fm_subdivide works for S2 meshes", {
  mesh_s2 <- fm_rcdt_2d_inla(globe = 0)
  # 4-division twice:
  mesh_s2_1 <- fm_subdivide(mesh_s2, 1)
  mesh_s2_1_1 <- fm_subdivide(mesh_s2_1, 1)
  # 3-division once:
  mesh_s2_3 <- fm_subdivide(mesh_s2, 3)

  # Should have points on the same sphere:
  expect_equal(
    mean(rowSums(mesh_s2_1$loc^2)),
    mean(rowSums(mesh_s2$loc^2))
  )
  expect_equal(
    mean(rowSums(mesh_s2_1_1$loc^2)),
    mean(rowSums(mesh_s2$loc^2))
  )
  expect_equal(
    mean(rowSums(mesh_s2_3$loc^2)),
    mean(rowSums(mesh_s2$loc^2))
  )

  expect_equal(fm_dof(mesh_s2_1), 12 + 20 * 3 * 1 / 2)
  expect_equal(fm_dof(mesh_s2_1_1), 12 + 20 * (3 * 3 / 2 + 3))
  expect_equal(fm_dof(mesh_s2_3), 12 + 20 * (3 * 3 / 2 + 3))

  # Should have the same number of nodes, but not the same locations:
  expect_equal(fm_dof(mesh_s2_3), fm_dof(mesh_s2_1_1))
})
