test_that("fm_nonconvex_hull works", {
  mat <- matrix(c(1:6, 6:1), 6, 2, byrow = TRUE)
  expect_error(
    {
      hull <- fm_nonconvex_hull(mat)
    },
    NA
  )

  expect_s3_class(hull, "sfc_POLYGON")

  expect_error(
    {
      hull <- fm_extensions(mat)
    },
    NA
  )

  expect_s3_class(hull[[1]], "sfc_POLYGON")
})

test_that("fm_nonconvex_hull_inla works", {
  mat <- matrix(c(1:6, 6:1), 6, 2, byrow = TRUE)
  expect_error(
    {
      hull <- fm_nonconvex_hull_inla(mat)
    },
    NA
  )

  expect_s3_class(hull, "fm_segm")
})
