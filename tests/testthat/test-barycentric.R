test_that("point location", {
  mesh_loc <- rbind(
    c(0, 0),
    c(1, 0),
    c(1, 1),
    c(0, 1)
  )
  mesh_tv <- rbind(
    c(1L, 2L, 3L),
    c(3L, 4L, 1L)
  )
  loc <- cbind(
    c(0.8, 0.2),
    c(0.1, 0.8)
  )
  b <-
    fmesher_bary(
      mesh_loc = mesh_loc,
      mesh_tv = mesh_tv - 1L,
      loc = loc,
      options = list()
    )

  expect_equal(
    rbind(
      b$bary[1, , drop = FALSE] %*%
        mesh_loc[mesh_tv[b$t[1] + 1L, ], , drop = FALSE],
      b$bary[2, , drop = FALSE] %*%
        mesh_loc[mesh_tv[b$t[2] + 1L, ], , drop = FALSE]
    ),
    loc
  )
})

test_that("point location", {
  local_fm_safe_inla()

  mesh_loc <- rbind(
    c(0, 0),
    c(1, 0),
    c(1, 1),
    c(0, 1)
  )
  mesh_tv <- rbind(
    c(1L, 2L, 3L),
    c(3L, 4L, 1L)
  )
  mesh <- INLA::inla.mesh.create(loc = mesh_loc, tv = mesh_tv)
  loc <- cbind(
    c(0.8, 0.2),
    c(0.1, 0.8)
  )
  b <- fm_bary(mesh, loc)

  expect_equal(
    rbind(
      b$bary[1, , drop = FALSE] %*%
        mesh_loc[mesh_tv[b$t[1], ], , drop = FALSE],
      b$bary[2, , drop = FALSE] %*%
        mesh_loc[mesh_tv[b$t[2], ], , drop = FALSE]
    ),
    loc
  )
})
