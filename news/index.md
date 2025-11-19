# Changelog

## fmesher (development version)

### New features

- Add
  [`fm_subset()`](https://inlabru-org.github.io/fmesher/reference/fm_subset.md)
  method for constructing a subset of a mesh based on a set of triangle
  (for `fm_mesh_2d`) or tetrahedron (for `fm_mesh_3d`) indices. (version
  `0.5.0.9003`)
- Add
  [`fm_zm()`](https://inlabru-org.github.io/fmesher/reference/fm_zm.md)
  method for adding/removing/unifying the Z/M dimensions of coordinate
  matrices and `sf` objects (version `0.5.0.9010`)
- New method
  [`fm_int_object()`](https://inlabru-org.github.io/fmesher/reference/fm_int_object.md)
  to construct tibbles with the same output format as the
  [`fm_int()`](https://inlabru-org.github.io/fmesher/reference/fm_int.md)
  method, for user-defined integration schemes (version `0.5.0.9011`)

### Improved features

- Add `bary=fm_bary()` information to
  [`fm_subdivide()`](https://inlabru-org.github.io/fmesher/reference/fm_subdivide.md)
  output, mapping the new mesh locations to the original mesh locations,
  e.g. for interpolating functions from the original mesh to the new
  mesh (version `0.5.0.9002`)
- Speed up
  [`fm_int()`](https://inlabru-org.github.io/fmesher/reference/fm_int.md)
  for polygons by bulk pre-computing
  [`fm_bary()`](https://inlabru-org.github.io/fmesher/reference/fm_bary.md)
  information instead of separate calls in
  [`fm_vertex_projection()`](https://inlabru-org.github.io/fmesher/reference/fm_vertex_projection.md)
  (version `0.5.0.9004`)
- Handle heterogeneous `sf` geometry XY/XYZ dimensions in
  [`fm_bary()`](https://inlabru-org.github.io/fmesher/reference/fm_bary.md)/[`fm_basis()`](https://inlabru-org.github.io/fmesher/reference/fm_basis.md)
  via
  [`fm_zm()`](https://inlabru-org.github.io/fmesher/reference/fm_zm.md)
  method that is called by `fm_onto_mesh()` and
  [`fm_unify_coords()`](https://inlabru-org.github.io/fmesher/reference/fm_unify_coords.md)
  to promote XY to XYZ when needed, before calling
  [`sf::st_coordinates()`](https://r-spatial.github.io/sf/reference/st_coordinates.html),
  as
  [`sf::st_coordinates()`](https://r-spatial.github.io/sf/reference/st_coordinates.html)
  otherwise fails. (version `0.5.0.9005`)
- Allow
  [`fm_int.fm_mesh_1d()`](https://inlabru-org.github.io/fmesher/reference/fm_int.md)
  to handle lists of matrices (for interval integration) and vectors
  (sums over point sets), for more flexible blockwise integration and
  summation schemes. (version `0.5.0.9007`)
- Allow
  [`fm_int.list()`](https://inlabru-org.github.io/fmesher/reference/fm_int.md)
  to include non-domain variables in the output object, e.g. for
  including per-transect covariates in the integration scheme. (version
  `0.5.0.9008`)
- Allow
  [`fm_int()`](https://inlabru-org.github.io/fmesher/reference/fm_int.md)
  for numeric/character/factor/`fm_mesh_1d` to handle nested list
  samplers (version `0.5.0.9009`)

### Bug fixes

- Make
  [`fm_subdivide()`](https://inlabru-org.github.io/fmesher/reference/fm_subdivide.md)
  store the indexing information for the original mesh locations in
  `$idx$loc` (version `0.5.0.9001`)
- Correct off-by-one indexing error in `$graph$vt` triangle indices
  (version `0.5.0.9002`)
- Propagate correct `crs` information in
  [`fm_bary_loc()`](https://inlabru-org.github.io/fmesher/reference/fm_bary_loc.md)
  for 2D spaces (version `0.5.0.9005`)
- Fix bugs in
  [`fm_mesh_intersection()`](https://inlabru-org.github.io/fmesher/reference/fm_mesh_intersection.md)
  and
  [`fm_bary()`](https://inlabru-org.github.io/fmesher/reference/fm_bary.md)
  C++ code that caused incorrect behaviour for locating points on meshes
  on subsets of the sphere. Also allow
  [`fm_mesh_intersection()`](https://inlabru-org.github.io/fmesher/reference/fm_mesh_intersection.md)
  to generate non-Delaunay triangles, allowing the generated meshes to
  be used for stable integration schemes (version `0.5.0.9006`)

## fmesher 0.5.0

CRAN release: 2025-07-07

### New features

- Add
  [`fm_area()`](https://inlabru-org.github.io/fmesher/reference/fm_area.md)
  method for `fm_segm` area calculations (version `0.4.0.9002`)
- Make
  [`fm_nonconvex_hull()`](https://inlabru-org.github.io/fmesher/reference/fm_nonconvex_hull.md)
  a unified function for “fm” and “sf” construction methods and
  input/output formats, with default method “fm” and output format “sf”,
  making
  [`fm_nonconvex_hull_inla()`](https://inlabru-org.github.io/fmesher/reference/fm_nonconvex_hull_inla.md)
  deprecated. (version `0.4.0.9004`)

### Improved features

- Rename
  [`fm_mesh_components()`](https://inlabru-org.github.io/fmesher/reference/fm_components.md)
  to
  [`fm_components()`](https://inlabru-org.github.io/fmesher/reference/fm_components.md)
  and add support for `fm_segm` objects (version `0.4.0.9001`)
- Add support for polygon output in
  [`fm_as_sfc()`](https://inlabru-org.github.io/fmesher/reference/fm_as_sfc.md)
  for closed boundary `fm_segm` objects (version `0.4.0.9002`)
- Add support for `fm_segm` as boundary input to
  [`fm_hexagon_lattice()`](https://inlabru-org.github.io/fmesher/reference/fm_hexagon_lattice.md)
  (version `0.4.0.9005`)
- Add `fm_segm` integration support for `fm_mesh_2d` objects in
  [`fm_int()`](https://inlabru-org.github.io/fmesher/reference/fm_int.md)
  (version `0.4.0.9005`)

### Bug fixes

- Fix indexing bug in
  [`fm_basis.fm_mesh_1d()`](https://inlabru-org.github.io/fmesher/reference/fm_basis.md)
  for `degree = 2` and `NA` locations (version `0.4.0.9003`)
- Detect unnamed `sfc` objects in
  [`fm_int()`](https://inlabru-org.github.io/fmesher/reference/fm_int.md)
  for multi-domain integration and give an error, as the user must
  provide a geometry name (version `0.4.0.9005`)
- Check that every
  [`fm_int()`](https://inlabru-org.github.io/fmesher/reference/fm_int.md)
  sampler for domain lists has at least one corresponding domain
  (version `0.4.0.9007`)
- Remove support for `character` `.block` information in
  [`fm_int()`](https://inlabru-org.github.io/fmesher/reference/fm_int.md)/[`fm_cprod()`](https://inlabru-org.github.io/fmesher/reference/fm_cprod.md),
  forcing the use of `integer`, as `character` could lead to incorrect
  block aggregation output ordering. The equivalent information is now
  available in `.block_origin`. Also clean up of
  [`fm_cprod()`](https://inlabru-org.github.io/fmesher/reference/fm_cprod.md)
  and `sf` handling (version `0.4.0.9006`)

## fmesher 0.4.0

CRAN release: 2025-06-12

### New features

- Add
  [`fm_hexagon_lattice()`](https://inlabru-org.github.io/fmesher/reference/fm_hexagon_lattice.md)
  for creating regular hexagonal lattice points to use with
  [`fm_mesh_2d()`](https://inlabru-org.github.io/fmesher/reference/fm_mesh_2d.md),
  from Man Ho Suen (version `0.3.0.9001` and `0.3.0.9006`)
- Add
  [`fm_mesh_components()`](https://inlabru-org.github.io/fmesher/reference/fm_components.md)
  method for extracting (dis)connected components from 2D and 3D meshes
  (version `0.3.0.9005`)
- Add
  [`fm_collect()`](https://inlabru-org.github.io/fmesher/reference/fm_collect.md)
  method for creating a collection of meshes of the same type (version
  `0.3.0.9009`)

### Improved features

- Construct better representative boundary points for `mid` data of
  `fm_mesh_1d` for `degree = 2, boundary = "free"` (version
  `0.3.0.9002`)
- Add argument `delaunay` to allow non-Delaunay mesh construction in
  [`fm_subdivide()`](https://inlabru-org.github.io/fmesher/reference/fm_subdivide.md)
  (version `0.3.0.9003`)
- Better handling of line integration when triangle edges and line
  transects are co-linear (version `0.3.0.9004`)
- Add support for `order > 2` for
  [`fm_fem.fm_mesh_1d()`](https://inlabru-org.github.io/fmesher/reference/fm_fem.md)
  (version `0.3.0.9011`) and generally for
  [`fm_matern_precision()`](https://inlabru-org.github.io/fmesher/reference/fm_gmrf.md)
  (version `0.3.0.9012`)
- Add argument `units` to
  [`fm_crs()`](https://inlabru-org.github.io/fmesher/reference/fm_crs.md)
  to allow setting the length unit for the CRS on creation/extraction
  (version `0.3.0.9013`)

### Bug fixes

- Improved bug fix in
  [`fm_mesh_2d_inla()`](https://inlabru-org.github.io/fmesher/reference/fm_mesh_2d.md)
  from `0.2.0.9014` to allow negative offset in the second layer when
  the inner layer is specified only through a boundary polygon (version
  `0.3.0.9007`)
- Bug fix for
  [`fm_mesh_2d_inla()`](https://inlabru-org.github.io/fmesher/reference/fm_mesh_2d.md)
  to ensure S2 manifold meshes covering the entire sphere aren’t
  accidentally limited to a subset of the globe (version `0.3.0.9008`)
- Handle colour palettes with fewer than requested colours in
  [`fm_generate_colors()`](https://inlabru-org.github.io/fmesher/reference/fm_generate_colors.md)
  (version `0.3.0.9008`)
- Fix bug in
  [`fm_detect_manifold()`](https://inlabru-org.github.io/fmesher/reference/fm_detect_manifold.md)
  that caused it to return `"R2"` for `"S2"` and `"M2"` manifolds
  (version `0.3.0.9010`)
- Fix bug in
  [`fm_as_mesh_3d_list()`](https://inlabru-org.github.io/fmesher/reference/fm_as_mesh_3d.md)
  that caused it to return the mesh surface as an `fm_mesh_2d` object
  instead of the full `fm_mesh_3d` object (version `0.3.0.9014`)

## fmesher 0.3.0

CRAN release: 2025-03-18

### New features

- Convert
  [`fm_bary()`](https://inlabru-org.github.io/fmesher/reference/fm_bary.md)
  output to a dedicated `fm_bary` class, with elements ‘index’ (index of
  the containing simplex) and ‘where’ (matrix of barycentric weights).
  (version `0.2.0.9001`)
- Add
  [`fm_bary_simplex()`](https://inlabru-org.github.io/fmesher/reference/fm_bary_simplex.md)
  generic method to extract the simplex vertex indices for an `fm_bary`
  object. (version `0.2.0.9001`)
- Add
  [`fm_bary_loc()`](https://inlabru-org.github.io/fmesher/reference/fm_bary_loc.md)
  generic method for converting `fm_bary` information to Euclidean
  coordinates (version `0.2.0.9005`)
- Add support for `fm_mesh_3d` and new `fm_lattice_Nd` class (version
  `0.2.0.9008` and `0.2.0.9011`)
- Add
  [`fm_assess()`](https://inlabru-org.github.io/fmesher/reference/fm_assess.md)
  method, replicating the old `INLA::inla.mesh.assessment()` method
  (version `0.2.0.9010`)
- Add R implementation of
  [`fm_qinv()`](https://inlabru-org.github.io/fmesher/reference/fm_qinv.md)
  for computing sparse matrix partial inverses (version `0.2.0.9010`)

### Improved features

- Handle `NA` location inputs to
  [`fm_basis.fm_mesh_1d()`](https://inlabru-org.github.io/fmesher/reference/fm_basis.md)
  (version `0.2.0.9002`)
- Simplify `fm_basis` object creation, and add `custom_classes`
  developers vignette (version `0.2.0.9004`)
- Add `fm_basis(mesh, loc = fm_bary(mesh, ...))` support for
  `fm_mesh_1d` and `fm_mesh_2d` objects (version `0.2.0.9006`)
- Add [`list()`](https://rdrr.io/r/base/list.html) input support for
  `fm_int.numeric`, to allow multiple integration blocks for discrete
  domains (version `0.2.0.9012`)
- Add `mappings` and `defs` support to the
  [`geom_fm.fm_mesh_1d()`](https://inlabru-org.github.io/fmesher/reference/geom_fm.md)
  method, allowing separate control of the `ggplot2` aesthetics for
  basis functions, knots, and function evaluations (version
  `0.2.0.9013`)
- Add support for `character` block input to `fm_block` methods, to
  automate multi-domain integration support from
  [`fm_int()`](https://inlabru-org.github.io/fmesher/reference/fm_int.md)
  (version `0.2.0.9017`)

### Bug fixes

- Make
  [`fm_try_callstack()`](https://inlabru-org.github.io/fmesher/reference/call-stack.md)
  more robust against large callstack sizes; solves “C stack” crash
  issue for `inla()` error reporting (version `0.2.0.9007`)
- Fix
  [`fm_bbox.fm_mesh_2d()`](https://inlabru-org.github.io/fmesher/reference/fm_bbox.md)
  bug for `"S2"` and `"M2"` manifold meshes (version `0.2.0.9008`)
- Fix bug in
  [`fm_rcdt_2d_inla()`](https://inlabru-org.github.io/fmesher/reference/fm_rcdt_2d.md)
  that lead to ignoring the `quality.spec` argument (version
  `0.2.0.9009`)
- Fix bug in
  [`fm_mesh_2d_inla()`](https://inlabru-org.github.io/fmesher/reference/fm_mesh_2d.md)
  that improperly ignored negative `offset` values when creating meshes
  based on only boundary information (version `0.2.0.9014`)
- Minor bugfix for
  [`print.fm_segm()`](https://inlabru-org.github.io/fmesher/reference/fmesher-print.md)
  for empty `fm_segm` objects (version `0.2.0.9015`)
- Bugfix for `fm_mesh_2d_inla(interior = ...)` where part of the code
  incorrectly assumed `interior` would be a list. Now converts a list
  into a single `fm_segm` object (version `0.2.0.9016`)

### Deprecation updates

- Remove long deprecated `inla.mesh` etc legacy methods; need to
  explicitly convert old objects. Retaining the `inla.mesh` etc class
  suffixes for now. (version `0.2.0.9001`)

## fmesher 0.2.0

CRAN release: 2024-11-06

### New methods

- Add print methods for `fm_basis` and `fm_evaluator` objects (version
  `0.1.7.9003`)
- Add
  [`fm_manifold_get()`](https://inlabru-org.github.io/fmesher/reference/fm_manifold.md)
  generic method to extract manifold information from general objects,
  so that external objects can implement their own manifold information
  storage (version `0.1.7.9005`)
- Add
  [`fm_crs_plot()`](https://inlabru-org.github.io/fmesher/reference/fm_crs_plot.md)
  method for plotting `fm_crs` objects with optional graticules and
  Tissot indicatrices (version `0.1.7.9010`)

### New method options

- Add `full` argument to
  [`fm_basis()`](https://inlabru-org.github.io/fmesher/reference/fm_basis.md)
  to toggle between matrix and full `fm_basis` object output (version
  `0.1.7.9002`)
- Add `loc` plotting option to `geom_fm.fm_mesh_2d` and modify the
  `ggplot` mapping interface for interior and boundary segments in the
  same method (version `0.1.7.9009`)
- Add `format="loc"` argument to
  [`fm_as_sfc.fm_mesh_2d()`](https://inlabru-org.github.io/fmesher/reference/fm_as_sfc.md)
  for converting mesh node coordinates to `sfc_POINT` format (version
  `0.1.7.9009`)

### Minor updates

- Update documentation and vector coordinate inputs to
  [`fm_lattice_2d()`](https://inlabru-org.github.io/fmesher/reference/fm_lattice_2d.md)
  to clarify input interpretation and ensure correct boundary
  orientation (version `0.1.7.9001`)
- Add some length unit handling to
  [`fm_crs_bounds()`](https://inlabru-org.github.io/fmesher/reference/fm_crs_wkt.md)
  (version `0.1.7.9008`)
- Add control argument `max_batch_size` to
  [`fm_bary.fm_mesh_2d()`](https://inlabru-org.github.io/fmesher/reference/fm_bary.md),
  that can be supplied via
  [`fm_basis()`](https://inlabru-org.github.io/fmesher/reference/fm_basis.md),
  for optional override of the default maximal batch calculation size,
  see [\#14](https://github.com/inlabru-org/fmesher/issues/14) (version
  `0.1.7.9011`)

### Deprecation updates

- Remove `sp` objects from `fmexample` data. Use
  [`fmexample_sp()`](https://inlabru-org.github.io/fmesher/reference/fmexample_sp.md)
  to access them if needed (version `0.1.7.9004`)
- Move `sp` dependency to Suggests, and remove `inlabru` dependency
  (version `0.1.7.9006`)
- Further `sp` use protection (version `0.1.7.9007`)
- Increased deprecation warning and error messages for old unsupported
  methods

## fmesher 0.1.7

CRAN release: 2024-07-01

- Fix regression bug in
  [`fm_manifold()`](https://inlabru-org.github.io/fmesher/reference/fm_manifold.md)
  that made it ignore all but the first given type options. See
  [\#16](https://github.com/inlabru-org/fmesher/issues/16) (version
  \`0.1.6.9001)
- Fix `plot.fm_mesh_2d` vectorisation bug (version `0.1.6.9002`)
- Add new
  [`fm_subdivide()`](https://inlabru-org.github.io/fmesher/reference/fm_subdivide.md)
  method for `fm_mesh_2d` meshes (version `0.1.6.9003`)

## fmesher 0.1.6

CRAN release: 2024-06-14

- Fix for hiding away-facing triangle edges in
  [`plot.fm_mesh_2d()`](https://inlabru-org.github.io/fmesher/reference/plot.fm_mesh_2d.md)
  and
  [`lines.fm_segm()`](https://inlabru-org.github.io/fmesher/reference/plot.fm_segm.md).
- Fix duplicated
  [`fm_unify_coords.sf()`](https://inlabru-org.github.io/fmesher/reference/fm_unify_coords.md)
  method. Thanks to Pablo Paccioretti for debugging and reporting, see
  [\#13](https://github.com/inlabru-org/fmesher/issues/13) (version
  `0.1.5.9001`)
- Use batched
  [`fm_bary.fm_mesh_2d()`](https://inlabru-org.github.io/fmesher/reference/fm_bary.md)
  computations for data sizes of `2e5` and upwards. This improves
  performance of
  [`fm_basis()`](https://inlabru-org.github.io/fmesher/reference/fm_basis.md)/`inla.spde.make.A()`
  for large data sets, see
  [\#14](https://github.com/inlabru-org/fmesher/issues/14) (version
  `0.1.5.9002`)
- Automatically convert raw `sfc` samplers objects to `sf` objects in
  [`fm_int.list()`](https://inlabru-org.github.io/fmesher/reference/fm_int.md)
  (version `0.1.5.9003`)
- Detect and warn about unsupported use of `fm_segm` objects in
  [`fm_int.list()`](https://inlabru-org.github.io/fmesher/reference/fm_int.md)
  (version `0.1.5.9004`)
- Add
  [`fm_basis()`](https://inlabru-org.github.io/fmesher/reference/fm_basis.md)
  and
  [`fm_fem()`](https://inlabru-org.github.io/fmesher/reference/fm_fem.md)
  support for `fm_tensor` function spaces (version `0.1.5.9005`)
- Add
  [`fm_CRS()`](https://inlabru-org.github.io/fmesher/reference/fm_CRS_sp.md)
  support for `terra` objects (version `0.1.5.9006`)

## fmesher 0.1.5

CRAN release: 2023-12-20

- Fix bug in
  [`fm_matern_sample()`](https://inlabru-org.github.io/fmesher/reference/fm_gmrf.md)
  that prevented the use of a non-NULL `loc` argument. For earlier
  versions, the workaround is to make a separate call to
  [`fm_evaluate()`](https://inlabru-org.github.io/fmesher/reference/fm_evaluate.md)
  after calling
  [`fm_matern_sample()`](https://inlabru-org.github.io/fmesher/reference/fm_gmrf.md).
- Improve speed of
  [`fm_block_log_shift()`](https://inlabru-org.github.io/fmesher/reference/fm_block.md)
  by an order of magnitude for multi-block cases.
- Fix bug in `plot.fm_mesh_2d` for meshes with 2D coordinate storage
  (version `0.1.4.9002`)

## fmesher 0.1.4

CRAN release: 2023-10-28

- Work around `std::get<variant>()` lack of support for MacOS `< 10.14`.

## fmesher 0.1.3

CRAN release: 2023-10-18

- Fix bug in
  [`fm_rcdt_2d_inla()`](https://inlabru-org.github.io/fmesher/reference/fm_rcdt_2d.md)
  that gave different defaults for the `extend` and `refine` arguments
  when equal to `TRUE`, that should be treated the same as
  [`list()`](https://rdrr.io/r/base/list.html) (version `0.1.2.9001`)
- Fix bug in
  [`fm_rcdt_2d_inla()`](https://inlabru-org.github.io/fmesher/reference/fm_rcdt_2d.md)
  to properly update `n` basis counter element when removing unused
  vertices. (version `0.1.2.9002`)
- Fix bug giving a spurious warning about inconsistent `is.bnd` for
  empty `fm_segm` objects, and inconsistent `grp` vector lengths.
  (version `0.1.2.9003`)
- Convert some of the old potentially unsafe C++ pointer methods to type
  safe C++17 features

## fmesher 0.1.2

CRAN release: 2023-08-25

- Various documentation improvements, in particular for INLA
  compatibility
- Modify mesh refinement tests to directly check the refinement criteria
  instead of the specific mesh result, to account for differences in
  floating point behaviour on M1 processors.
- Modify tests of non-mesh-generation features to use precomputed meshes
  or meshes with stable properties
- Protect against invalid `tv` inputs
- Revert from `\text{}` to `\textrm{}`, as AMS extensions are only
  supported from R 4.2.2
  (<https://www.stats.bris.ac.uk/R/doc/manuals/r-devel/R-exts.pdf>
  2023-08-24, page 90), and CRAN oldrel for macOS is 4.2.0, not 4.2.3

## fmesher 0.1.1

CRAN release: 2023-08-18

- Simplify LICENSE information in the R package DESCRIPTION
- Fix example for legacy `plot_PolySet` method
- Fix C++ unused-variable warning detected by clang-tidy version 16.0.6
- Add workaround for
  [`sf::st_buffer`](https://r-spatial.github.io/sf/reference/geos_unary.html)
  not supporting negative distance for `s2` with longlat coordinates.
  Fixes [\#5](https://github.com/inlabru-org/fmesher/issues/5)

## fmesher 0.1.0

- Full conversion of the fmesher interface from the standalone binary in
  R-INLA (<https://www.r-inla.org>) to an interface powered by `Rcpp`
  (<https://cran.r-project.org/package=Rcpp>)
- Tools for 2D and 1D function spaces, see
  [`fm_mesh_2d()`](https://inlabru-org.github.io/fmesher/reference/fm_mesh_2d.md)
  and
  [`fm_mesh_1d()`](https://inlabru-org.github.io/fmesher/reference/fm_mesh_1d.md)
- Tools for CRS handling see
  [`fm_crs()`](https://inlabru-org.github.io/fmesher/reference/fm_crs.md)
  and
  [`fm_transform()`](https://inlabru-org.github.io/fmesher/reference/fm_transform.md)
- Plotting support for base graphics, `rgl`, and `ggplot2`
- Added basic GMRF tools, see
  [`fm_matern_precision()`](https://inlabru-org.github.io/fmesher/reference/fm_gmrf.md)

## fmesher 0.0.9

- Basic fmesher library I/O interface
- Added a `NEWS.md` file to track changes to the package.
