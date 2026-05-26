# Interactive mesh building and diagnostics

Assess the finite element approximation errors in a mesh for interactive
R sessions.

## Usage

``` r
fm_assess(mesh, spatial.range, alpha = 2, dims = NULL)
```

## Arguments

- mesh:

  An
  [fm_mesh_2d](https://inlabru-org.github.io/fmesher/reference/fm_mesh_2d.md)
  object

- spatial.range:

  numeric; the spatial range parameter to use for the assessment

- alpha:

  numeric; A valid
  [`fm_matern_precision()`](https://inlabru-org.github.io/fmesher/reference/fm_gmrf.md)
  `alpha` parameter

- dims:

  2-numeric; the grid size

## Value

An `sf` object with gridded mesh assessment information

## See also

[`fm_mesh_2d()`](https://inlabru-org.github.io/fmesher/reference/fm_mesh_2d.md),
[fm_rcdt_2d](https://inlabru-org.github.io/fmesher/reference/fm_rcdt_2d.md)

## Author

Finn Lindgren <Finn.Lindgren@gmail.com>

## Examples

``` r

bnd <- fm_segm(cbind(
  c(0, 10, 10, 0, 0),
  c(0, 0, 10, 10, 0)
), is.bnd = TRUE)
mesh <- fm_rcdt_2d_inla(boundary = bnd, max.edge = 1)
out <- fm_assess(mesh, spatial.range = 3, alpha = 2)
```
