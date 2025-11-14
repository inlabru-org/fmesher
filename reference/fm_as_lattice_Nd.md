# Convert objects to `fm_lattice_Nd`

Convert objects to `fm_lattice_Nd`

## Usage

``` r
fm_as_lattice_Nd(...)

fm_as_lattice_Nd_list(x, ...)

# S3 method for class 'fm_lattice_Nd'
fm_as_lattice_Nd(x, ...)
```

## Arguments

- ...:

  Arguments passed on to submethods

- x:

  Object to be converted

## Value

An `fm_lattice_Md` or `fm_lattice_Nd_list` object

## Functions

- `fm_as_lattice_Nd()`: Convert an object to `fm_lattice_Nd`.

- `fm_as_lattice_Nd_list()`: Convert each element of a list

## See also

Other object creation and conversion:
[`fm_as_collect()`](https://inlabru-org.github.io/fmesher/reference/fm_as_collect.md),
[`fm_as_fm()`](https://inlabru-org.github.io/fmesher/reference/fm_as_fm.md),
[`fm_as_lattice_2d()`](https://inlabru-org.github.io/fmesher/reference/fm_as_lattice_2d.md),
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
[`fm_mesh_2d()`](https://inlabru-org.github.io/fmesher/reference/fm_mesh_2d.md),
[`fm_segm()`](https://inlabru-org.github.io/fmesher/reference/fm_segm.md),
[`fm_simplify()`](https://inlabru-org.github.io/fmesher/reference/fm_simplify.md),
[`fm_tensor()`](https://inlabru-org.github.io/fmesher/reference/fm_tensor.md)

## Examples

``` r
(fm_as_lattice_Nd_list(list(
  fm_lattice_Nd(list(1:3, 1:2)),
  fm_lattice_Nd(list(1:4))
)))
#> fm_lattice_Nd object:
#>   Manifold:  R2
#>   Dimensions:    3 x 2
#>   Bounding box: (1,3) x (1,2)
#>   Basis d.o.f.:  6
#> fm_lattice_Nd object:
#>   Manifold:  R1
#>   Dimensions:    4
#>   Bounding box: (1,4)
#>   Basis d.o.f.:  4
```
