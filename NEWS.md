# fmesher (development version)

* Fix bug in `fm_rcdt_2d_inla()` that gave different defaults for the `extend` and
  `refine` arguments when equal to `TRUE`, that should be treated the same as `list()`
  (version `0.1.2.9001`)
* Fix bug in `fm_rcdt_2d_inla()` to properly update `n` basis counter element
  when removing unused vertices. (version `0.1.2.9002`)

# fmesher 0.1.2

* Various documentation improvements, in particular for INLA compatibility
* Modify mesh refinement tests to directly check the refinement criteria
  instead of the specific mesh result, to account for differences in
  floating point behaviour on M1 processors.
* Modify tests of non-mesh-generation features to use precomputed meshes
  or meshes with stable properties
* Protect against invalid `tv` inputs
* Revert from `\text{}` to `\textrm{}`, as AMS extensions are only supported
  from R 4.2.2 (https://www.stats.bris.ac.uk/R/doc/manuals/r-devel/R-exts.pdf
  2023-08-24, page 90), and CRAN oldrel for macOS is 4.2.0, not 4.2.3

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
