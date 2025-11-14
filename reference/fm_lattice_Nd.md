# Lattice grids for N dimensions

Construct an N-dimensional lattice grid

## Usage

``` r
fm_lattice_Nd(x = NULL, ...)

# S3 method for class 'matrix'
fm_lattice_Nd(x = NULL, dims = NULL, values = NULL, ...)

# S3 method for class 'data.frame'
fm_lattice_Nd(x = NULL, ...)

# S3 method for class 'list'
fm_lattice_Nd(x = NULL, dims = NULL, ...)

# S3 method for class 'fm_bbox'
fm_lattice_Nd(x = NULL, dims = NULL, ...)

# S3 method for class '`NULL`'
fm_lattice_Nd(x = NULL, ..., dims = NULL)
```

## Arguments

- x:

  `list`, `data.frame`, `matrix`, `fm_bbox` or `NULL`. If a list of
  vectors, `as.matrix(expand.grid(x))` is used to create a full grid
  coordinates. `data.frame` and `matrix` input is assumed to follow the
  same ordering convention as the output of
  [`expand.grid()`](https://rdrr.io/r/base/expand.grid.html). of length
  N of vectors or grid matrices of coordinate values. List vector values
  are sorted before use.

- ...:

  Passed on to submethods

- dims:

  numeric; the size of the grid of dimension `length(dims)`

- values:

  list of grid axis values

## Value

An `fm_lattice_Nd` object with elements

- dims:

  integer vector

- values:

  the grid coordinate axis values

- loc:

  matrix of constructed grid coordinates

## Methods (by class)

- `` fm_lattice_Nd(`NULL`) ``: Ignores the `NULL` `x` and creates a
  lattice based on `values` (if non-NULL) and `dims` unit hypercube
  lattice grid with `dims` dimensions.

## See also

[`fm_mesh_3d()`](https://inlabru-org.github.io/fmesher/reference/fm_mesh_3d.md)

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
[`fm_mesh_1d()`](https://inlabru-org.github.io/fmesher/reference/fm_mesh_1d.md),
[`fm_mesh_2d()`](https://inlabru-org.github.io/fmesher/reference/fm_mesh_2d.md),
[`fm_segm()`](https://inlabru-org.github.io/fmesher/reference/fm_segm.md),
[`fm_simplify()`](https://inlabru-org.github.io/fmesher/reference/fm_simplify.md),
[`fm_tensor()`](https://inlabru-org.github.io/fmesher/reference/fm_tensor.md)

## Author

Finn Lindgren <Finn.Lindgren@gmail.com>

## Examples

``` r
(lattice <- fm_lattice_Nd(
  list(
    seq(0, 1, length.out = 3),
    seq(0, 1, length.out = 4),
    seq(0, 1, length.out = 2)
  )
))
#> fm_lattice_Nd object:
#>   Manifold:  R3
#>   Dimensions:    3 x 4 x 2
#>   Bounding box: (0,1) x (0,1) x (0,1)
#>   Basis d.o.f.:  24

if (requireNamespace("geometry", quietly = TRUE)) {
  (mesh <- fm_delaunay_3d(lattice$loc))
}
#> fm_mesh_3d object:
#>   Manifold:  R3
#>   V / E / T / Tet:   24 / 84 / 100 / 36
#>   Euler char.:   4
#>   Bounding box: (0,1) x (0,1) x (0,1)
#>   Basis d.o.f.:  24
```
