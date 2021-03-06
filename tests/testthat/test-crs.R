local_fm_testthat_setup()

test_that("Basic CRS predefinitions", {
  expect_warning({crs <- fm_CRS("longlat")},
                 "Use of old predefined projection")
  expect_s4_class(crs, "CRS")

  crs <- fm_CRS("longlat_globe")
  expect_s4_class(crs, "CRS")
})
