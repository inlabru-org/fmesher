test_that("base graphics plotting", {
  skip_on_cran()
  skip_if_not(interactive())

  plot(fmexample$mesh)

  plot(fmexample$boundary_fm)
})

test_that("rgl plotting", {
  skip_on_cran()
  skip_if_not(interactive() && FALSE)
  skip_if_not_installed("rgl")

  plot_rgl(fmexample$mesh)

  plot_rgl(fmexample$boundary_fm)
})
