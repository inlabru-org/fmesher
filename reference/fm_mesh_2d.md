# Make a 2D mesh object

Make a 2D mesh object

## Usage

``` r
fm_mesh_2d(...)

fm_mesh_2d_inla(
  loc = NULL,
  loc.domain = NULL,
  offset = NULL,
  n = NULL,
  boundary = NULL,
  interior = NULL,
  max.edge = NULL,
  min.angle = NULL,
  cutoff = 1e-12,
  max.n.strict = NULL,
  max.n = NULL,
  plot.delay = NULL,
  crs = NULL,
  ...
)
```

## Arguments

- ...:

  Currently passed on to `fm_mesh_2d_inla`

- loc:

  Matrix of point locations to be used as initial triangulation nodes.
  Can alternatively be a `sf`, `sfc`, `SpatialPoints` or
  `SpatialPointsDataFrame` object.

- loc.domain:

  Matrix of point locations used to determine the domain extent. Can
  alternatively be a `SpatialPoints` or `SpatialPointsDataFrame` object.

- offset:

  The automatic extension distance. One or two values, for an inner and
  an optional outer extension. If negative, interpreted as a factor
  relative to the approximate data diameter (default=-0.10???)

- n:

  The number of initial nodes in the automatic extensions (default=16)

- boundary:

  one or more (as list) of
  [`fm_segm()`](https://inlabru-org.github.io/fmesher/reference/fm_segm.md)
  objects, or objects supported by
  [`fm_as_segm()`](https://inlabru-org.github.io/fmesher/reference/fm_as_segm.md)

- interior:

  one object supported by
  [`fm_as_segm()`](https://inlabru-org.github.io/fmesher/reference/fm_as_segm.md),
  or (from version `0.2.0.9016`) a list of such objects. If a list, the
  objects are joined into a single object.

- max.edge:

  The largest allowed triangle edge length. One or two values.

- min.angle:

  The smallest allowed triangle angle. One or two values. (Default=21)

- cutoff:

  The minimum allowed distance between points. Point at most as far
  apart as this are replaced by a single vertex prior to the mesh
  refinement step.

- max.n.strict:

  The maximum number of vertices allowed, overriding `min.angle` and
  `max.edge` (default=-1, meaning no limit). One or two values, where
  the second value gives the number of additional vertices allowed for
  the extension.

- max.n:

  The maximum number of vertices allowed, overriding `max.edge` only
  (default=-1, meaning no limit). One or two values, where the second
  value gives the number of additional vertices allowed for the
  extension.

- plot.delay:

  If logical `TRUE` or a negative numeric value, activates displaying
  the result after each step of the multi-step domain extension
  algorithm.

- crs:

  An optional
  [`fm_crs()`](https://inlabru-org.github.io/fmesher/reference/fm_crs.md),
  [`sf::crs`](https://r-spatial.github.io/sf/reference/coerce-methods.html)
  or [`sp::CRS`](https://edzer.github.io/sp/reference/CRS-class.html)
  object

## Value

An `fm_mesh_2d` object.

## Functions

- `fm_mesh_2d_inla()`: Legacy method for `INLA::inla.mesh.2d()` Create a
  triangle mesh based on initial point locations, specified or automatic
  boundaries, and mesh quality parameters.

## INLA compatibility

For mesh and curve creation, the
[`fm_rcdt_2d_inla()`](https://inlabru-org.github.io/fmesher/reference/fm_rcdt_2d.md),
`fm_mesh_2d_inla()`, and
[`fm_nonconvex_hull_inla()`](https://inlabru-org.github.io/fmesher/reference/fm_nonconvex_hull_inla.md)
methods will keep the interface syntax used by
`INLA::inla.mesh.create()`, `INLA::inla.mesh.2d()`, and
`INLA::inla.nonconvex.hull()` functions, respectively, whereas the
[`fm_rcdt_2d()`](https://inlabru-org.github.io/fmesher/reference/fm_rcdt_2d.md),
`fm_mesh_2d()`, and
[`fm_nonconvex_hull()`](https://inlabru-org.github.io/fmesher/reference/fm_nonconvex_hull.md)
interfaces may be different, and potentially change in the future.

## See also

[`fm_rcdt_2d()`](https://inlabru-org.github.io/fmesher/reference/fm_rcdt_2d.md),
`fm_mesh_2d()`,
[`fm_delaunay_2d()`](https://inlabru-org.github.io/fmesher/reference/fm_rcdt_2d.md),
[`fm_nonconvex_hull()`](https://inlabru-org.github.io/fmesher/reference/fm_nonconvex_hull.md),
[`fm_extensions()`](https://inlabru-org.github.io/fmesher/reference/fm_nonconvex_hull.md),
[`fm_refine()`](https://inlabru-org.github.io/fmesher/reference/fm_refine.md)

Other object creation and conversion:
[`fm_as_collect()`](https://inlabru-org.github.io/fmesher/reference/fm_as_collect.md),
[`fm_as_fm()`](https://inlabru-org.github.io/fmesher/reference/fm_as_fm.md),
[`fm_as_lattice_2d()`](https://inlabru-org.github.io/fmesher/reference/fm_as_lattice_2d.md),
[`fm_as_lattice_Nd()`](https://inlabru-org.github.io/fmesher/reference/fm_as_lattice_Nd.md),
[`fm_as_mesh_1d()`](https://inlabru-org.github.io/fmesher/reference/fm_as_mesh_1d.md),
[`fm_as_mesh_2d()`](https://inlabru-org.github.io/fmesher/reference/fm_as_mesh_2d.md),
[`fm_as_mesh_3d()`](https://inlabru-org.github.io/fmesher/reference/fm_as_mesh_3d.md),
[`fm_as_segm()`](https://inlabru-org.github.io/fmesher/reference/fm_as_segm.md),
[`fm_as_sfc()`](https://inlabru-org.github.io/fmesher/reference/fm_as_sfc.md),
[`fm_as_tensor()`](https://inlabru-org.github.io/fmesher/reference/fm_as_tensor.md),
[`fm_collect()`](https://inlabru-org.github.io/fmesher/reference/fm_collect.md),
[`fm_lattice_2d()`](https://inlabru-org.github.io/fmesher/reference/fm_lattice_2d.md),
[`fm_lattice_Nd()`](https://inlabru-org.github.io/fmesher/reference/fm_lattice_Nd.md),
[`fm_mesh_1d()`](https://inlabru-org.github.io/fmesher/reference/fm_mesh_1d.md),
[`fm_segm()`](https://inlabru-org.github.io/fmesher/reference/fm_segm.md),
[`fm_simplify()`](https://inlabru-org.github.io/fmesher/reference/fm_simplify.md),
[`fm_tensor()`](https://inlabru-org.github.io/fmesher/reference/fm_tensor.md)

## Author

Finn Lindgren <Finn.Lindgren@gmail.com>

## Examples

``` r
fm_mesh_2d_inla(boundary = fm_extensions(cbind(2, 1), convex = 1, 2))
#> fm_mesh_2d object:
#>   Manifold:  R2
#>   V / E / T: 32 / 69 / 38
#>   Euler char.:   1
#>   Constraints:   Boundary: 24 boundary edges (1 group: 1), Interior: 0 edges
#>   Bounding box: (1.002964,2.997036) x (0.002963687,1.997036313)
#>   Basis d.o.f.:  32
```
