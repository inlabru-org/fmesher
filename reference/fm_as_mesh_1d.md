# Convert objects to `fm_segm`

Convert objects to `fm_segm`

## Usage

``` r
fm_as_mesh_1d(x, ...)

fm_as_mesh_1d_list(x, ...)

# S3 method for class 'fm_mesh_1d'
fm_as_mesh_1d(x, ...)

# S3 method for class 'inla.mesh.1d'
fm_as_mesh_1d(x, ...)
```

## Arguments

- x:

  Object to be converted

- ...:

  Arguments passed on to submethods

## Value

An `fm_mesh_1d` or `fm_mesh_1d_list` object

## Functions

- `fm_as_mesh_1d()`: Convert an object to `fm_mesh_1d`.

- `fm_as_mesh_1d_list()`: Convert each element of a list

## See also

Other object creation and conversion:
[`fm_as_collect()`](https://inlabru-org.github.io/fmesher/reference/fm_as_collect.md),
[`fm_as_fm()`](https://inlabru-org.github.io/fmesher/reference/fm_as_fm.md),
[`fm_as_lattice_2d()`](https://inlabru-org.github.io/fmesher/reference/fm_as_lattice_2d.md),
[`fm_as_lattice_Nd()`](https://inlabru-org.github.io/fmesher/reference/fm_as_lattice_Nd.md),
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
fm_as_mesh_1d_list(list(fm_mesh_1d(1:4)))
#> fm_mesh_1d object:
#>   Manifold:  R1
#>   #{knots}:  4
#>   Interval:  (1, 4)
#>   Boundary:  (neumann, neumann)
#>   B-spline degree:   1
#>   Basis d.o.f.:  4
```
