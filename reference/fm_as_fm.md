# Convert objects to fmesher objects

Used for conversion from general objects (usually `inla.mesh` and other
legacy INLA specific classes) to `fmesher` classes.

## Usage

``` r
fm_as_fm(x, ...)

## S3 method for class 'NULL'
fm_as_fm(x, ...)

# S3 method for class 'fm_mesh_1d'
fm_as_fm(x, ...)

# S3 method for class 'fm_mesh_2d'
fm_as_fm(x, ...)

# S3 method for class 'fm_mesh_3d'
fm_as_fm(x, ...)

# S3 method for class 'fm_tensor'
fm_as_fm(x, ...)

# S3 method for class 'fm_collect'
fm_as_fm(x, ...)

# S3 method for class 'fm_segm'
fm_as_fm(x, ...)

# S3 method for class 'fm_lattice_Nd'
fm_as_fm(x, ...)

# S3 method for class 'fm_lattice_2d'
fm_as_fm(x, ...)

# S3 method for class 'fm_bbox'
fm_as_fm(x, ...)

# S3 method for class 'crs'
fm_as_fm(x, ...)

# S3 method for class 'CRS'
fm_as_fm(x, ...)

# S3 method for class 'fm_crs'
fm_as_fm(x, ...)

# S3 method for class 'inla.CRS'
fm_as_fm(x, ...)

# S3 method for class 'inla.mesh.1d'
fm_as_fm(x, ...)

# S3 method for class 'inla.mesh'
fm_as_fm(x, ...)

# S3 method for class 'inla.mesh.segment'
fm_as_fm(x, ...)

# S3 method for class 'inla.mesh.lattice'
fm_as_fm(x, ...)
```

## Arguments

- x:

  Object to be converted

- ...:

  Arguments forwarded to submethods

## Value

An object of some `fm_*` class

## See also

Other object creation and conversion:
[`fm_as_collect()`](https://inlabru-org.github.io/fmesher/reference/fm_as_collect.md),
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
[`fm_mesh_2d()`](https://inlabru-org.github.io/fmesher/reference/fm_mesh_2d.md),
[`fm_segm()`](https://inlabru-org.github.io/fmesher/reference/fm_segm.md),
[`fm_simplify()`](https://inlabru-org.github.io/fmesher/reference/fm_simplify.md),
[`fm_tensor()`](https://inlabru-org.github.io/fmesher/reference/fm_tensor.md)

## Examples

``` r
fm_as_fm(NULL)
#> NULL
```
