# Non-convex hull computation

**\[deprecated\]** Legacy method for `INLA::inla.nonconvex.hull()`. Use
[`fm_nonconvex_hull()`](https://inlabru-org.github.io/fmesher/reference/fm_nonconvex_hull.md)
with `method = "fm"` instead, with either `format = "fm"` (for
compatibility with code expecting `fm_segm` output) or `format = "sf"`.

## Usage

``` r
fm_nonconvex_hull_inla(
  x,
  convex = -0.15,
  concave = convex,
  resolution = 40,
  eps = NULL,
  eps_rel = NULL,
  crs = NULL,
  ...
)

fm_nonconvex_hull_inla_basic(
  x,
  convex = -0.15,
  resolution = 40,
  eps = NULL,
  crs = fm_crs(x)
)
```

## Arguments

- x:

  A spatial object

- convex:

  numeric vector; How much to extend

- concave:

  numeric vector; The minimum allowed reentrant curvature. Default equal
  to `convex`

- resolution:

  integer; The internal computation resolution. A warning will be issued
  when this needs to be increased for higher accuracy, with the required
  resolution stated. For `method="fm"` only.

- eps, eps_rel:

  The polygonal curve simplification tolerances used for simplifying the
  resulting boundary curve. See
  [`fm_simplify_helper()`](https://inlabru-org.github.io/fmesher/reference/fm_simplify_helper.md)
  for details. For `method="fm"` only.

- crs:

  Optional crs object for the resulting polygon. Default is `fm_crs(x)`

- ...:

  Unused.

## Value

`fm_nonconvex_hull_inla()` returns an
[fm_segm](https://inlabru-org.github.io/fmesher/reference/fm_segm.md)
object, for compatibility with `inla.nonconvex.hull()`.

## Functions

- `fm_nonconvex_hull_inla_basic()`: Special method
  [`fm_nonconvex_hull_fm()`](https://inlabru-org.github.io/fmesher/reference/fm_nonconvex_hull.md)
  method for `concave = 0`. Requires
  [`splancs::nndistF()`](https://rsbivand.github.io/splancs/reference/nndistF.html).

## INLA compatibility

For mesh and curve creation, the
[`fm_rcdt_2d_inla()`](https://inlabru-org.github.io/fmesher/reference/fm_rcdt_2d.md),
[`fm_mesh_2d_inla()`](https://inlabru-org.github.io/fmesher/reference/fm_mesh_2d.md),
and `fm_nonconvex_hull_inla()` methods will keep the interface syntax
used by `INLA::inla.mesh.create()`, `INLA::inla.mesh.2d()`, and
`INLA::inla.nonconvex.hull()` functions, respectively, whereas the
[`fm_rcdt_2d()`](https://inlabru-org.github.io/fmesher/reference/fm_rcdt_2d.md),
[`fm_mesh_2d()`](https://inlabru-org.github.io/fmesher/reference/fm_mesh_2d.md),
and
[`fm_nonconvex_hull()`](https://inlabru-org.github.io/fmesher/reference/fm_nonconvex_hull.md)
interfaces may be different, and potentially change in the future.

## See also

[`fm_nonconvex_hull()`](https://inlabru-org.github.io/fmesher/reference/fm_nonconvex_hull.md)

Other nonconvex inla legacy support:
[`fm_segm_contour_helper()`](https://inlabru-org.github.io/fmesher/reference/fm_segm_contour_helper.md),
[`fm_simplify_helper()`](https://inlabru-org.github.io/fmesher/reference/fm_simplify_helper.md)

## Examples

``` r
# New preferred method for "fm_segm" output:
fm_nonconvex_hull(cbind(0, 0), convex = 1, format = "fm")
#> fm_segm object:
#>   32 boundary edges (1 group: 1)
#>   Bounding box = (-0.9992483, 0.9992483) x (-0.9992483, 0.9992483) x (0,0)

# Deprecated:
suppressWarnings(
  fm_nonconvex_hull_inla(cbind(0, 0), convex = 1)
)
#> fm_segm object:
#>   32 boundary edges (1 group: 1)
#>   Bounding box = (-0.9992483, 0.9992483) x (-0.9992483, 0.9992483) x (0,0)
```
