# Convert objects to `fm_mesh_3d`

Convert objects to `fm_mesh_3d`

## Usage

``` r
fm_as_mesh_3d(x, ...)

fm_as_mesh_3d_list(x, ...)

# S3 method for class 'fm_mesh_3d'
fm_as_mesh_3d(x, ...)
```

## Arguments

- x:

  Object to be converted

- ...:

  Arguments passed on to submethods

## Value

An `fm_mesh_3d` or `fm_mesh_3d_list` object

## Functions

- `fm_as_mesh_3d()`: Convert an object to `fm_mesh_3d`.

- `fm_as_mesh_3d_list()`: Convert each element of a list

## See also

Other object creation and conversion:
[`fm_as_collect()`](https://inlabru-org.github.io/fmesher/reference/fm_as_collect.md),
[`fm_as_fm()`](https://inlabru-org.github.io/fmesher/reference/fm_as_fm.md),
[`fm_as_lattice_2d()`](https://inlabru-org.github.io/fmesher/reference/fm_as_lattice_2d.md),
[`fm_as_lattice_Nd()`](https://inlabru-org.github.io/fmesher/reference/fm_as_lattice_Nd.md),
[`fm_as_mesh_1d()`](https://inlabru-org.github.io/fmesher/reference/fm_as_mesh_1d.md),
[`fm_as_mesh_2d()`](https://inlabru-org.github.io/fmesher/reference/fm_as_mesh_2d.md),
[`fm_as_segm()`](https://inlabru-org.github.io/fmesher/reference/fm_as_segm.md),
[`fm_as_sfc()`](https://inlabru-org.github.io/fmesher/reference/fm_as_sfc.md),
[`fm_as_tensor()`](https://inlabru-org.github.io/fmesher/reference/fm_as_tensor.md),
[`fm_collect()`](https://inlabru-org.github.io/fmesher/reference/fm_collect.md),
[`fm_lattice_2d()`](https://inlabru-org.github.io/fmesher/reference/fm_lattice_2d.md),
[`fm_lattice_Nd()`](https://inlabru-org.github.io/fmesher/reference/fm_lattice_Nd.md),
[`fm_mesh_1d()`](https://inlabru-org.github.io/fmesher/reference/fm_mesh_1d.md),
[`fm_mesh_2d()`](https://inlabru-org.github.io/fmesher/reference/fm_mesh_2d.md),
[`fm_segm()`](https://inlabru-org.github.io/fmesher/reference/fm_segm.md),
[`fm_simplify()`](https://inlabru-org.github.io/fmesher/reference/fm_simplify.md),
[`fm_tensor()`](https://inlabru-org.github.io/fmesher/reference/fm_tensor.md)

## Examples

``` r
(m <- fm_mesh_3d(
  matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0), 4, 3, byrow = TRUE),
  matrix(c(1, 2, 3, 4), 1, 4, byrow = TRUE)
))
#> fm_mesh_3d object:
#>   Manifold:  R3
#>   V / E / T / Tet:   4 / 6 / 4 / 1
#>   Euler char.:   1
#>   Bounding box: (0,1) x (0,1) x (0,1)
#>   Basis d.o.f.:  4
fm_as_mesh_3d_list(list(m))
#> fm_mesh_3d object:
#>   Manifold:  R3
#>   V / E / T / Tet:   4 / 6 / 4 / 1
#>   Euler char.:   1
#>   Bounding box: (0,1) x (0,1) x (0,1)
#>   Basis d.o.f.:  4
```
