test_that("Call stack utilities", {
  fun_ok <- function() {
    1 + 1
  }
  fun_not_ok <- function() {
    stop("An error test")
  }

  expect_equal(fm_try_callstack(fun_ok()), 2)
  expect_s3_class(fm_try_callstack(fun_not_ok()), "try-error")
})

test_that("Package requirements handling", {
  expect_message(
    {
      # jarl-ignore internal_function: Internal test
      require_fmesher <- fmesher:::fm_require_message("fmesher")
    },
    NA
  )
  expect_true(require_fmesher)

  expect_message(
    {
      # jarl-ignore internal_function: Internal test
      require_unknown <-
        fmesher:::fm_require_message("this_package_really_should_not_exist")
    },
    "Please install 'this_package_really_should_not_exist'"
  )
  expect_false(require_unknown)

  expect_error(
    {
      # jarl-ignore internal_function: Internal test
      fmesher:::fm_require_stop("this_package_really_should_not_exist")
    },
    "Please install 'this_package_really_should_not_exist'"
  )
})
