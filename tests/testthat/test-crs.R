test_that("CRS/WKT predefinitions", {
  crs <- fm_crs("longlat_globe")
  expect_s3_class(crs, "crs")

  skip_if_not(fm_safe_sp())

  expect_warning(
    {
      crs <- fm_CRS("longlat")
    },
    "Use of old predefined projection"
  )
  expect_s4_class(crs, "CRS")

  crs <- fm_CRS("longlat_globe")
  expect_s4_class(crs, "CRS")
})


test_that("crs object creation", {
  crs <- fm_crs("+proj=longlat +R=1 +no_defs")
  expect_s3_class(crs, "crs")

  crs <- fm_crs(NA_character_)
  expect_s3_class(crs, "crs")
})

test_that("fm_crs object creation", {
  crs <- fm_crs("+proj=longlat +R=1 +no_defs", oblique = 1:4)
  expect_s3_class(crs, "fm_crs")

  crs <- fm_crs(NA_character_, oblique = 1:2)
  expect_s3_class(crs, "fm_crs")

  crs <- fm_crs("globe", oblique = 1:4)
  expect_s3_class(crs, "fm_crs")
})
