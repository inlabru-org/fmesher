# Make a lattice object

Construct a lattice grid for
[`fm_mesh_2d()`](https://inlabru-org.github.io/fmesher/reference/fm_mesh_2d.md)

## Usage

``` r
fm_lattice_2d(...)

# Default S3 method
fm_lattice_2d(
  x = seq(0, 1, length.out = 2),
  y = seq(0, 1, length.out = 2),
  z = NULL,
  dims = if (is.matrix(x)) {
     dim(x)
 } else {
     c(length(x), length(y))
 },
  units = NULL,
  crs = NULL,
  ...
)
```

## Arguments

- ...:

  Passed on to submethods

- x:

  vector or grid matrix of x-values. Vector values are sorted before
  use. Matrix input is assumed to be a grid of x-values with the same
  ordering convention of `as.vector(x)` as `rep(x, times = dims[2])` for
  vector input.

- y:

  vector of grid matrix of y-values. Vector values are sorted before
  use. Matrix input is assumed to be a grid of y-values with the same
  ordering convention of `as.vector(y)` as `rep(y, each = dims[1])` for
  vector input.

- z:

  if x is a matrix, a grid matrix of z-values, with the same ordering as
  `x` and `y`. If `x` is a vector, `z` is ignored.

- dims:

  the size of the grid, length 2 vector

- units:

  One of `c("default", "longlat", "longsinlat", "mollweide")` or NULL
  (equivalent to `"default"`).

- crs:

  An optional `fm_crs`,
  [`sf::st_crs`](https://r-spatial.github.io/sf/reference/st_crs.html),
  or [`sp::CRS`](https://edzer.github.io/sp/reference/CRS-class.html)
  object, denoting the CRS info for the x-y grid.

## Value

An `fm_lattice_2d` object with elements

- dims:

  integer vector

- x:

  x-values for original vector input

- y:

  y-values for original vector input

- loc:

  matrix of `(x, y)` values or `(x, y, z)` values. May be altered by
  [`fm_transform()`](https://inlabru-org.github.io/fmesher/reference/fm_transform.md)

- segm:

  `fm_segm` object

- crs:

  `fm_crs` object for `loc`, or `NULL`

- crs0:

  `fm_crs` object for `(x,y)`, or `NULL`

## See also

[`fm_mesh_2d()`](https://inlabru-org.github.io/fmesher/reference/fm_mesh_2d.md)

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
[`fm_lattice_Nd()`](https://inlabru-org.github.io/fmesher/reference/fm_lattice_Nd.md),
[`fm_mesh_1d()`](https://inlabru-org.github.io/fmesher/reference/fm_mesh_1d.md),
[`fm_mesh_2d()`](https://inlabru-org.github.io/fmesher/reference/fm_mesh_2d.md),
[`fm_segm()`](https://inlabru-org.github.io/fmesher/reference/fm_segm.md),
[`fm_simplify()`](https://inlabru-org.github.io/fmesher/reference/fm_simplify.md),
[`fm_tensor()`](https://inlabru-org.github.io/fmesher/reference/fm_tensor.md)

## Author

Finn Lindgren <Finn.Lindgren@gmail.com>

## Examples

``` r
lattice <- fm_lattice_2d(
  seq(0, 1, length.out = 17),
  seq(0, 1, length.out = 10)
)

## Use the lattice "as-is", without refinement:
mesh <- fm_rcdt_2d_inla(lattice = lattice, boundary = lattice$segm)
mesh <- fm_rcdt_2d_inla(lattice = lattice, extend = FALSE)

## Refine the triangulation, with limits on triangle angles and edges:
mesh <- fm_rcdt_2d(
  lattice = lattice,
  refine = list(max.edge = 0.08),
  extend = FALSE
)

## Add an extension around the lattice, but maintain the lattice edges:
mesh <- fm_rcdt_2d(
  lattice = lattice,
  refine = list(max.edge = 0.08),
  interior = lattice$segm
)

## Only add extension:
mesh <- fm_rcdt_2d(lattice = lattice, refine = list(max.edge = 0.08))
```
