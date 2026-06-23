# Contour segment

Helper from legacy `INLA::inla.contour.segment()`

## Usage

``` r
fm_segm_contour_helper(
  x = seq(0, 1, length.out = nrow(z)),
  y = seq(0, 1, length.out = ncol(z)),
  z,
  nlevels = 10,
  levels = NULL,
  groups = NULL,
  positive = TRUE,
  eps = NULL,
  eps_rel = NULL,
  crs = NULL
)
```

## Arguments

- x, y, z:

  The `x` and `y` coordinates of the grid, and the `z` values

- nlevels:

  Number of contour levels

- levels:

  The contour levels. If `NULL`,
  `pretty(range(z, na.rm = TRUE), nlevels)` is used.

- groups:

  The group values for each contour level. If `NULL`,
  `seq_len(length(levels))` is used.

- positive:

  Logical; if `TRUE`, the contour lines are made to be CCW around
  positive excursions

- eps, eps_rel:

  Polygonal curve simplification tolerances

- crs:

  A coordinate reference system

## Value

An `fm_segm` object

## See also

Other nonconvex inla legacy support:
[`fm_simplify_helper()`](https://inlabru-org.github.io/fmesher/reference/fm_simplify_helper.md),
[`fmesher-deprecated`](https://inlabru-org.github.io/fmesher/reference/fmesher-deprecated.md)

## Examples

``` r
fm_segm_contour_helper(z = matrix(1:16, 4, 4))
#> fm_segm object:
#>   8 interior edges (8 groups: 2 3 4 5 6 7 8 9)
#>   Bounding box = (0,1) x (0,1) x (0,0)
```
