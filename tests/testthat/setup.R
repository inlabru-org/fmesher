tryCatch(
  {
    local_fm_testthat_setup(envir = testthat::teardown_env())
  },
  error = function(e) {
    local_fm_testthat_setup(envir = .GlobalEnv)
  }
)
