# fmesher 0.1.1

# fmesher 0.1.1

* Simplify LICENSE information in the R package DESCRIPTION

* Fix example for legacy `plot_PolySet` method

* Fix C++ unused-variable warning detected by clang-tidy version 16.0.6

* Add workaround for `sf::st_buffer` not supporting negative distance
  for `s2` with longlat coordinates. Fixes #5

# fmesher 0.1.0

* Full conversion of the fmesher interface from the standalone binary in R-INLA
  (https://www.r-inla.org) to an interface powered by `Rcpp`
  (https://cran.r-project.org/package=Rcpp)

* Tools for 2D and 1D function spaces, see `fm_mesh_2d()` and `fm_mesh_1d()`

* Tools for CRS handling see `fm_crs()` and `fm_transform()`

* Plotting support for base graphics, `rgl`, and `ggplot2`

* Added basic GMRF tools, see `fm_matern_precision()`

# fmesher 0.0.9

* Basic fmesher library I/O interface

* Added a `NEWS.md` file to track changes to the package.
