# Converting legacy INLA mesh code to fmesher

## Deprecation of old methods and class names

Great effort has been taken to preserve backwards compatibility as far
as practical, with in particular the old `inla.mesh`, `inla.mesh.1d`,
and `inla.mesh.segment` object classes given fallback methods that carry
out methods for the new `fm_mesh_2d`, `fm_mesh_1d`, and `fm_segm`
classes. Starting in **November 2024** however, some of this direct
fallback support is being phased out, so that old stored objects may
need to be explicitly converted to `fmesher` objects using for example
[`fm_as_fm()`](https://inlabru-org.github.io/fmesher/reference/fm_as_fm.md).
New code, in particular in packages that use `fmesher` objects, should
use the new interface methods, and replace references to `inla.mesh`,
`inla.mesh.1d`, and `inla.mesh.segment` with `fm_mesh_2d`, `fm_mesh_1d`,
and `fm_segm`, respectively, as the old class names will eventually be
dropped from the mesh classes. This in particular applies to S3 method
naming and class checking, where `inherits(mesh, "inla.mesh")` must be
replaced with `inherits(mesh, "fm_mesh_2d")` in order to work in the
future.

## Deprecation warnings

Three `INLA::inla.getOption()`/`INLA::inla.setOption()` options control
the transition behaviour of the `INLA` package use of `fmesher`.

- `fmesher.evolution`, integer:

  - `1L` uses the intermediate `fm_*` methods in `fmesher` that were
    already available via inlabru from 2.8.0, but calls the `INLA`
    built-in `fmesher` standalone programme for mesh construction and
    related operations. (From `INLA` version 23.06.29)
  - `2L` uses the full range of `fmesher` package methods, and does not
    use the standalone `fmesher` programme. (From `INLA` around version
    23.08.20, and is now the only supported option.)

- `fmesher.evolution.warn`, logical: When `TRUE`, `INLA` will show
  deprecation methods for all the methods in `INLA` that are being
  replaced by `fmesher` package methods. When `FALSE`, no warnings will
  be shown. Set this option to `TRUE` if you want to update your own
  code, but keep it at `FALSE` when you need to run existing code
  without changing it.

- `fmesher.evolution.verbosity`, character: Either “soft”, “warn”, or
  “stop”, indicating the minimum warning level when
  `fmesher.evolution.warn` is `TRUE`. Set this to “warn” or “stop” to
  get more immediate feedback when testing conversion of old code,
  e.g. in package testing.

Packages using mesh methods should use these options in package testing,
e.g. in `tests/testthat/setup.R`:

``` r
if (requireNamespace("INLA", quietly = TRUE)) {
  INLA::inla.setOption(fmesher.evolution = 2L)
  INLA::inla.setOption(fmesher.evolution.warn = TRUE)
  INLA::inla.setOption(fmesher.evolution.verbosity = "warn")
}
```

### Compatibility

An important change is that the handling of mesh `crs` information is
now more flexible, but stricter in the sense that user code should make
no assumption about how the information is stored in the mesh, and
should therefore avoid the direct `mesh$crs` access, and instead use the
[`fm_crs()`](https://inlabru-org.github.io/fmesher/reference/fm_crs.md)
and
[`fm_CRS()`](https://inlabru-org.github.io/fmesher/reference/fm_CRS_sp.md)
access methods, depending on what type of CRS object is needed. The
ideal way to specify crs information is in the initial mesh creation
call. If the crs needs to be explicitly assigned a new value, use the
`fm_crs(mesh) <- crs` assignment method.

For mesh and curve creation, the
[`fm_rcdt_2d_inla()`](https://inlabru-org.github.io/fmesher/reference/fm_rcdt_2d.md),
[`fm_mesh_2d_inla()`](https://inlabru-org.github.io/fmesher/reference/fm_mesh_2d.md),
and
[`fm_nonconvex_hull_inla()`](https://inlabru-org.github.io/fmesher/reference/fm_nonconvex_hull_inla.md)
methods will keep the interface syntax used by `inla.mesh.create()`,
`inla.mesh.2d()`, and `inla.nonconvex.hull()` functions, respectively,
whereas the
[`fm_rcdt_2d()`](https://inlabru-org.github.io/fmesher/reference/fm_rcdt_2d.md),
[`fm_mesh_2d()`](https://inlabru-org.github.io/fmesher/reference/fm_mesh_2d.md),
and
[`fm_nonconvex_hull()`](https://inlabru-org.github.io/fmesher/reference/fm_nonconvex_hull.md)
interfaces may change in the future.

## Mesh construction

| INLA                                              | fmesher                                                                                                                                                                                                                                                                            | Comments                                                                                                                                        |
|:--------------------------------------------------|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:------------------------------------------------------------------------------------------------------------------------------------------------|
| `inla.mesh.create()`                              | [`fm_rcdt_2d()`](https://inlabru-org.github.io/fmesher/reference/fm_rcdt_2d.html) , [`fm_rcdt_2d_inla()`](https://inlabru-org.github.io/fmesher/reference/fm_rcdt_2d.html)                                                                                                         |                                                                                                                                                 |
| `inla.mesh.2d()`                                  | [`fm_mesh_2d()`](https://inlabru-org.github.io/fmesher/reference/fm_mesh_2d.html) , [`fm_mesh_2d_inla()`](https://inlabru-org.github.io/fmesher/reference/fm_mesh_2d.html)                                                                                                         |                                                                                                                                                 |
| `inla.delaunay()`                                 | [`fm_delaunay_2d()`](https://inlabru-org.github.io/fmesher/reference/fm_rcdt_2d.html)                                                                                                                                                                                              |                                                                                                                                                 |
| `inla.mesh.1d()`                                  | [`fm_mesh_1d()`](https://inlabru-org.github.io/fmesher/reference/fm_mesh_1d.html)                                                                                                                                                                                                  |                                                                                                                                                 |
| `inla.mesh.lattice()`                             | [`fm_lattice_2d()`](https://inlabru-org.github.io/fmesher/reference/fm_lattice_2d.html)                                                                                                                                                                                            |                                                                                                                                                 |
| `inla.mesh.segment()`                             | [`fm_segm()`](https://inlabru-org.github.io/fmesher/reference/fm_segm.html)                                                                                                                                                                                                        |                                                                                                                                                 |
| `inla.nonconvex.hull()`                           | [`fm_nonconvex_hull()`](https://inlabru-org.github.io/fmesher/reference/fm_nonconvex_hull.html), [`fm_extensions()`](https://inlabru-org.github.io/fmesher/reference/fm_nonconvex_hull.html) , [`fm_simplify()`](https://inlabru-org.github.io/fmesher/reference/fm_simplify.html) | Use `format = "fm"` to get `fm_segm` output from [`fm_nonconvex_hull()`](https://inlabru-org.github.io/fmesher/reference/fm_nonconvex_hull.md). |
| `inla.contour.segment()`, `inla.simplify.curve()` | [`fm_simplify_helper()`](https://inlabru-org.github.io/fmesher/reference/fm_simplify_helper.html) , [`fm_segm_contour_helper()`](https://inlabru-org.github.io/fmesher/reference/fm_segm_contour_helper.html)                                                                      |                                                                                                                                                 |
| `inla.mesh.components()`                          | [`fm_mesh_components()`](https://inlabru-org.github.io/fmesher/reference/fmesher-deprecated.html)                                                                                                                                                                                  |                                                                                                                                                 |
| `NA`                                              | [`fm_subdivide()`](https://inlabru-org.github.io/fmesher/reference/fm_subdivide.html)                                                                                                                                                                                              |                                                                                                                                                 |
| `NA`                                              | [`fm_hexagon_lattice()`](https://inlabru-org.github.io/fmesher/reference/fm_hexagon_lattice.html)                                                                                                                                                                                  | Creates points for equilateral triangles inside a polygon domain.                                                                               |

## Location, basis and function evaluation

| INLA                    | fmesher                                                                                                                                                                                                                                                                                                                                                                                                               | inlabru                                                                                                                                                                                                                                                                                 |
|:------------------------|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `inla.mesh.query()`     | N/A                                                                                                                                                                                                                                                                                                                                                                                                                   | [`No fmesher equivalent in version <= 0.5.0::No fmesher equivalent in version <= 0.5.0No fmesher equivalent in version <= 0.5.0`](NA)                                                                                                                                                   |
| `inla.mesh.projector()` | [`fm_evaluator()`](https://inlabru-org.github.io/fmesher/reference/fm_evaluate.html)                                                                                                                                                                                                                                                                                                                                  |                                                                                                                                                                                                                                                                                         |
| `inla.mesh.project()`   | [`fm_evaluate()`](https://inlabru-org.github.io/fmesher/reference/fm_evaluate.html)                                                                                                                                                                                                                                                                                                                                   |                                                                                                                                                                                                                                                                                         |
| `inla.spde.make.A()`    | [`fm_bary()`](https://inlabru-org.github.io/fmesher/reference/fm_bary.html) , [`fm_basis()`](https://inlabru-org.github.io/fmesher/reference/fm_basis.html) , [`fm_row_kron()`](https://inlabru-org.github.io/fmesher/reference/fm_row_kron.html), [`fm_block()`](https://inlabru-org.github.io/fmesher/reference/fm_block.html) , [`fm_block_eval()`](https://inlabru-org.github.io/fmesher/reference/fm_block.html) | [`inlabru::bm_multi()`](https://inlabru-org.github.io/inlabru/reference/bm_multi.html) , [`inlabru::ibm_jacobian()`](https://inlabru-org.github.io/inlabru/reference/ibm_jacobian.html), [`inlabru::bm_aggregate()`](https://inlabru-org.github.io/inlabru/reference/bm_aggregate.html) |
| `inla.mesh.deriv()`     | [`fm_basis()`](https://inlabru-org.github.io/fmesher/reference/fm_basis.html)                                                                                                                                                                                                                                                                                                                                         |                                                                                                                                                                                                                                                                                         |

## Finite element methods

| INLA                                     | fmesher                                                                                 | Comments                                                                                                                                                                                                  |
|:-----------------------------------------|:----------------------------------------------------------------------------------------|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `inla.mesh.fem()` , `inla.mesh.1d.fem()` | [`fm_fem()`](https://inlabru-org.github.io/fmesher/reference/fm_fem.html)               |                                                                                                                                                                                                           |
| `NA`                                     | [`fm_matern_precision()`](https://inlabru-org.github.io/fmesher/reference/fm_gmrf.html) |                                                                                                                                                                                                           |
| `NA`                                     | [`fm_matern_sample()`](https://inlabru-org.github.io/fmesher/reference/fm_gmrf.html)    | Convenience function that combines [`fm_matern_precision()`](https://inlabru-org.github.io/fmesher/reference/fm_gmrf.md) and [`fm_sample()`](https://inlabru-org.github.io/fmesher/reference/fm_gmrf.md). |
| `NA`                                     | [`fm_covariance()`](https://inlabru-org.github.io/fmesher/reference/fm_gmrf.html)       | Basic helper function for computing covariances between different locations. Can produce sparse inverses like `inla.qinv()`, but currently (version 0.1.1) only by a ‘brute force’ method.                |
|                                          | [`fm_qinv()`](https://inlabru-org.github.io/fmesher/reference/fm_qinv.html)             | Produce sparse inverses like `inla.qinv()`, but currently (version 0.2.0.9010) by an R implementation.                                                                                                    |
|                                          | [`fm_sample()`](https://inlabru-org.github.io/fmesher/reference/fm_gmrf.html)           | Basic sampling method, like `inla.qsample()`                                                                                                                                                              |

## Printing

| INLA                  | fmesher                                                                                                                                                                                                                                                                          | Comments               |
|:----------------------|:---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:-----------------------|
| `summary.inla.mesh()` | [`print.fm_mesh_2d()`](https://inlabru-org.github.io/fmesher/reference/fmesher-print.html), [`print.fm_segm()`](https://inlabru-org.github.io/fmesher/reference/fmesher-print.html) , [`print.fm_mesh_1d()`](https://inlabru-org.github.io/fmesher/reference/fmesher-print.html) | Use `print(mesh)` etc. |

## CRS information and coordinate transformations

| INLA                 | fmesher                                                                                                                                                          | Comments                                                                                                                                                                                                                                                                                    |
|:---------------------|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------|:--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `inla.spTransform()` | [`fm_transform()`](https://inlabru-org.github.io/fmesher/reference/fm_transform.html)                                                                            |                                                                                                                                                                                                                                                                                             |
| `mesh$crs`           | [`fm_crs(mesh)`](https://inlabru-org.github.io/fmesher/reference/fm_crs.html) , [`fm_CRS(mesh)`](https://inlabru-org.github.io/fmesher/reference/fm_CRS_sp.html) | The crs may now be stored in different formats; use [`fm_crs()`](https://inlabru-org.github.io/fmesher/reference/fm_crs.md) for `sf` format, and [`fm_CRS()`](https://inlabru-org.github.io/fmesher/reference/fm_CRS_sp.md) for `sp` format. `fmesher` will attempt to convert when needed. |
| `mesh$crs<-`         | [`fm_crs(mesh)<-`](https://inlabru-org.github.io/fmesher/reference/fm_crs.html)                                                                                  | Direct assignment of crs information should be avoided, but is allowed as long as its compatible with the actual mesh coordinates.                                                                                                                                                          |

## Plotting

| INLA                                   | fmesher                                                                                                                                                                                    | inlabru                                                                        | Comments                                                                                                                                   |
|:---------------------------------------|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:-------------------------------------------------------------------------------|:-------------------------------------------------------------------------------------------------------------------------------------------|
| `No ggplot support`                    | [`geom_fm(data = mesh)`](https://inlabru-org.github.io/fmesher/reference/geom_fm.html), [`geom_fm(data = segm)`](https://inlabru-org.github.io/fmesher/reference/geom_fm.html)             | [`inlabru::gg(mesh)`](https://inlabru-org.github.io/inlabru/reference/gg.html) | Use `ggplot() + geom_fm(data = mesh)` and `inlabru::gg()` methods                                                                          |
| `plot.inla.mesh(rgl = FALSE)`          | [`plot.fm_mesh_2d()`](https://inlabru-org.github.io/fmesher/reference/plot.fm_mesh_2d.html) , [`lines.fm_mesh_2d()`](https://inlabru-org.github.io/fmesher/reference/plot.fm_mesh_2d.html) |                                                                                | Use [`plot()`](https://rspatial.github.io/terra/reference/plot.html) or [`lines()`](https://rspatial.github.io/terra/reference/lines.html) |
| `lines.inla.mesh.segment(rgl = FALSE)` | [`plot.fm_segm()`](https://inlabru-org.github.io/fmesher/reference/plot.fm_segm.html) , [`lines.fm_segm()`](https://inlabru-org.github.io/fmesher/reference/plot.fm_segm.html)             |                                                                                | Use [`plot()`](https://rspatial.github.io/terra/reference/plot.html) or [`lines()`](https://rspatial.github.io/terra/reference/lines.html) |
| `plot.inla.mesh(rgl = TRUE)`           | [`plot_rgl()`](https://inlabru-org.github.io/fmesher/reference/plot_rgl.html) , [`lines_rgl()`](https://inlabru-org.github.io/fmesher/reference/plot_rgl.html)                             |                                                                                |                                                                                                                                            |
| `lines.inla.mesh.segment(rgl = TRUE)`  | [`plot_rgl()`](https://inlabru-org.github.io/fmesher/reference/plot_rgl.html) , [`lines_rgl()`](https://inlabru-org.github.io/fmesher/reference/plot_rgl.html)                             |                                                                                |                                                                                                                                            |
