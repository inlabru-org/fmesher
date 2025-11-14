# Convert objects to `fm_mesh_2d`

Convert objects to `fm_mesh_2d`

## Usage

``` r
fm_as_mesh_2d(x, ...)

fm_as_mesh_2d_list(x, ...)

# S3 method for class 'fm_mesh_2d'
fm_as_mesh_2d(x, ...)

# S3 method for class 'inla.mesh'
fm_as_mesh_2d(x, ...)

# S3 method for class 'fm_mesh_3d'
fm_as_mesh_2d(x, ...)

# S3 method for class 'sfg'
fm_as_mesh_2d(x, ...)

# S3 method for class 'sfc_MULTIPOLYGON'
fm_as_mesh_2d(x, ...)

# S3 method for class 'sfc_POLYGON'
fm_as_mesh_2d(x, ...)

# S3 method for class 'sf'
fm_as_mesh_2d(x, ...)
```

## Arguments

- x:

  Object to be converted

- ...:

  Arguments passed on to submethods

## Value

An `fm_mesh_2d` or `fm_mesh_2d_list` object

## Methods (by class)

- `fm_as_mesh_2d(fm_mesh_3d)`: Construct a 2D mesh of the boundary of a
  3D mesh

## Functions

- `fm_as_mesh_2d()`: Convert an object to `fm_mesh_2d`.

- `fm_as_mesh_2d_list()`: Convert each element of a list

## See also

Other object creation and conversion:
[`fm_as_collect()`](https://inlabru-org.github.io/fmesher/reference/fm_as_collect.md),
[`fm_as_fm()`](https://inlabru-org.github.io/fmesher/reference/fm_as_fm.md),
[`fm_as_lattice_2d()`](https://inlabru-org.github.io/fmesher/reference/fm_as_lattice_2d.md),
[`fm_as_lattice_Nd()`](https://inlabru-org.github.io/fmesher/reference/fm_as_lattice_Nd.md),
[`fm_as_mesh_1d()`](https://inlabru-org.github.io/fmesher/reference/fm_as_mesh_1d.md),
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
fm_as_mesh_2d_list(list(fm_mesh_2d(cbind(2, 1))))
#> fm_mesh_2d object:
#>   Manifold:  R2
#>   V / E / T: 9 / 16 / 8
#>   Euler char.:   1
#>   Constraints:   Boundary: 8 boundary edges (1 group: 0), Interior: 0 edges
#>   Bounding box: (1,3) x (-4.440892e-16, 2.000000e+00)
#>   Basis d.o.f.:  9
```
