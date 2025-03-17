# fmesher 0.3.0

## New features

* Convert `fm_bary()` output to a dedicated `fm_bary` class, with elements 'index'
  (index of the containing simplex) and 'where' (matrix of barycentric weights).
  (version `0.2.0.9001`)
* Add `fm_bary_simplex()` generic method to extract the simplex vertex indices for an
  `fm_bary` object.
  (version `0.2.0.9001`)
* Add `fm_bary_loc()` generic method for converting `fm_bary` information to
  Euclidean coordinates (version `0.2.0.9005`)
* Add support for `fm_mesh_3d` and new `fm_lattice_Nd` class
  (version `0.2.0.9008` and `0.2.0.9011`)
* Add `fm_assess()` method, replicating the old `INLA::inla.mesh.assessment()`
  method (version `0.2.0.9010`)
* Add R implementation of `fm_qinv()` for computing sparse matrix partial
  inverses (version `0.2.0.9010`)
  
## Improved features

* Handle `NA` location inputs to `fm_basis.fm_mesh_1d()` (version `0.2.0.9002`)
* Simplify `fm_basis` object creation, and add `custom_classes` developers
  vignette (version `0.2.0.9004`)
* Add `fm_basis(mesh, loc = fm_bary(mesh, ...))` support for `fm_mesh_1d`
  and `fm_mesh_2d` objects (version `0.2.0.9006`)
* Add `list()` input support for `fm_int.numeric`, to allow multiple integration
  blocks for discrete domains (version `0.2.0.9012`)
* Add `mappings` and `defs` support to the `geom_fm.fm_mesh_1d()` method,
  allowing separate control of the `ggplot2` aesthetics for basis functions,
  knots, and function evaluations (version `0.2.0.9013`)
* Add support for `character` block input to `fm_block` methods, to automate
  multi-domain integration support from `fm_int()` (version `0.2.0.9017`)
  
## Bug fixes

* Make `fm_try_callstack()` more robust against large callstack sizes; solves
  "C stack" crash issue for `inla()` error reporting (version `0.2.0.9007`)
* Fix `fm_bbox.fm_mesh_2d()` bug for `"S2"` and `"M2"` manifold meshes (version `0.2.0.9008`)
* Fix bug in `fm_rcdt_2d_inla()` that lead to ignoring the `quality.spec` argument
  (version `0.2.0.9009`)
* Fix bug in `fm_mesh_2d_inla()` that improperly ignored negative `offset`
  values when creating meshes based on only boundary information (version `0.2.0.9014`)
* Minor bugfix for `print.fm_segm()` for empty `fm_segm` objects (version `0.2.0.9015`)
* Bugfix for `fm_mesh_2d_inla(interior = ...)` where part of the code incorrectly
  assumed `interior` would be a list. Now converts a list into a single `fm_segm`
  object (version `0.2.0.9016`)

## Deprecation updates

* Remove long deprecated `inla.mesh` etc legacy methods; need to explicitly
  convert old objects. Retaining the `inla.mesh` etc class suffixes for now.
  (version `0.2.0.9001`)

# fmesher 0.2.0

## New methods

* Add print methods for `fm_basis` and `fm_evaluator` objects (version `0.1.7.9003`)
* Add `fm_manifold_get()` generic method to extract manifold information from general
  objects, so that external objects can implement their own manifold information
  storage (version `0.1.7.9005`)
* Add `fm_crs_plot()` method for plotting `fm_crs` objects with optional
  graticules and Tissot indicatrices (version `0.1.7.9010`)

## New method options

* Add `full` argument to `fm_basis()` to toggle between matrix and full `fm_basis`
  object output (version `0.1.7.9002`)
* Add `loc` plotting option to `geom_fm.fm_mesh_2d` and modify the `ggplot`
  mapping interface for interior and boundary segments in the same method
  (version `0.1.7.9009`)
* Add `format="loc"` argument to `fm_as_sfc.fm_mesh_2d()` for converting mesh
  node coordinates to `sfc_POINT` format (version `0.1.7.9009`)

## Minor updates

* Update documentation and vector coordinate inputs to `fm_lattice_2d()`
  to clarify input interpretation and ensure correct boundary orientation (version `0.1.7.9001`)
* Add some length unit handling to `fm_crs_bounds()` (version `0.1.7.9008`)
* Add control argument `max_batch_size` to `fm_bary.fm_mesh_2d()`, that can be
  supplied via `fm_basis()`, for optional override of the default maximal batch
  calculation size, see #14 (version `0.1.7.9011`)

## Deprecation updates

* Remove `sp` objects from `fmexample` data. Use `fmexample_sp()` to access them if needed
  (version `0.1.7.9004`)
* Move `sp` dependency to Suggests, and remove `inlabru` dependency (version `0.1.7.9006`)
* Further `sp` use protection (version `0.1.7.9007`)
* Increased deprecation warning and error messages for old unsupported methods

# fmesher 0.1.7

* Fix regression bug in `fm_manifold()` that made it ignore all but the first given type options. See #16 (version `0.1.6.9001)
* Fix `plot.fm_mesh_2d` vectorisation bug (version `0.1.6.9002`)
* Add new `fm_subdivide()` method for `fm_mesh_2d` meshes (version `0.1.6.9003`)

# fmesher 0.1.6

* Fix for hiding away-facing triangle edges in `plot.fm_mesh_2d()` and `lines.fm_segm()`.
* Fix duplicated `fm_unify_coords.sf()` method. Thanks to Pablo Paccioretti for debugging and reporting, see #13 (version `0.1.5.9001`)
* Use batched `fm_bary.fm_mesh_2d()` computations for data sizes of `2e5` and upwards.
  This improves performance of `fm_basis()`/`inla.spde.make.A()` for large data sets, see #14 (version `0.1.5.9002`)
* Automatically convert raw `sfc` samplers objects to `sf` objects in `fm_int.list()` (version `0.1.5.9003`)
* Detect and warn about unsupported use of `fm_segm` objects in `fm_int.list()` (version `0.1.5.9004`)
* Add `fm_basis()` and `fm_fem()` support for `fm_tensor` function spaces (version `0.1.5.9005`)
* Add `fm_CRS()` support for `terra` objects (version `0.1.5.9006`)

# fmesher 0.1.5

* Fix bug in `fm_matern_sample()` that prevented the use of a non-NULL `loc` argument.
  For earlier versions, the workaround is to make a separate call to `fm_evaluate()`
  after calling `fm_matern_sample()`.
* Improve speed of `fm_block_log_shift()` by an order of magnitude for multi-block
  cases.
* Fix bug in `plot.fm_mesh_2d` for meshes with 2D coordinate storage (version `0.1.4.9002`)

# fmesher 0.1.4

* Work around `std::get<variant>()` lack of support for MacOS `< 10.14`.

# fmesher 0.1.3

* Fix bug in `fm_rcdt_2d_inla()` that gave different defaults for the `extend` and
  `refine` arguments when equal to `TRUE`, that should be treated the same as `list()`
  (version `0.1.2.9001`)
* Fix bug in `fm_rcdt_2d_inla()` to properly update `n` basis counter element
  when removing unused vertices. (version `0.1.2.9002`)
* Fix bug giving a spurious warning about inconsistent `is.bnd` for empty
  `fm_segm` objects, and inconsistent `grp` vector lengths. (version `0.1.2.9003`)
* Convert some of the old potentially unsafe C++ pointer methods to type safe
  C++17 features

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
