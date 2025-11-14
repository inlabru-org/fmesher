# Convert objects to `fm_lattice_2d`

Convert objects to `fm_lattice_2d`

## Usage

``` r
fm_as_lattice_2d(...)

fm_as_lattice_2d_list(x, ...)

# S3 method for class 'fm_lattice_2d'
fm_as_lattice_2d(x, ...)

# S3 method for class 'inla.mesh.lattice'
fm_as_lattice_2d(x, ...)
```

## Arguments

- ...:

  Arguments passed on to submethods

- x:

  Object to be converted

## Value

An `fm_lattice_2d` or `fm_lattice_2d_list` object

## Functions

- `fm_as_lattice_2d()`: Convert an object to `fm_lattice_2d`.

- `fm_as_lattice_2d_list()`: Convert each element of a list

## See also

Other object creation and conversion:
[`fm_as_collect()`](https://inlabru-org.github.io/fmesher/reference/fm_as_collect.md),
[`fm_as_fm()`](https://inlabru-org.github.io/fmesher/reference/fm_as_fm.md),
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
str(fm_as_lattice_2d_list(list(fm_lattice_2d(), fm_lattice_2d())))
#> List of 2
#>  $ :List of 7
#>   ..$ dims: int [1:2] 2 2
#>   ..$ x   : num [1:2] 0 1
#>   ..$ y   : num [1:2] 0 1
#>   ..$ loc : num [1:4, 1:2] 0 1 0 1 0 0 1 1
#>   ..$ segm:List of 5
#>   .. ..$ loc   : num [1:4, 1:3] 0 1 1 0 0 0 1 1 0 0 ...
#>   .. ..$ idx   : int [1:4, 1:2] 1 2 3 4 2 3 4 1
#>   .. ..$ grp   : int [1:4] 1 2 3 4
#>   .. ..$ is.bnd: logi [1:4] TRUE TRUE TRUE TRUE
#>   .. ..$ crs   : NULL
#>   .. ..- attr(*, "class")= chr [1:2] "fm_segm" "inla.mesh.segment"
#>   ..$ crs : NULL
#>   ..$ crs0: NULL
#>   ..- attr(*, "class")= chr [1:2] "fm_lattice_2d" "inla.mesh.lattice"
#>  $ :List of 7
#>   ..$ dims: int [1:2] 2 2
#>   ..$ x   : num [1:2] 0 1
#>   ..$ y   : num [1:2] 0 1
#>   ..$ loc : num [1:4, 1:2] 0 1 0 1 0 0 1 1
#>   ..$ segm:List of 5
#>   .. ..$ loc   : num [1:4, 1:3] 0 1 1 0 0 0 1 1 0 0 ...
#>   .. ..$ idx   : int [1:4, 1:2] 1 2 3 4 2 3 4 1
#>   .. ..$ grp   : int [1:4] 1 2 3 4
#>   .. ..$ is.bnd: logi [1:4] TRUE TRUE TRUE TRUE
#>   .. ..$ crs   : NULL
#>   .. ..- attr(*, "class")= chr [1:2] "fm_segm" "inla.mesh.segment"
#>   ..$ crs : NULL
#>   ..$ crs0: NULL
#>   ..- attr(*, "class")= chr [1:2] "fm_lattice_2d" "inla.mesh.lattice"
#>  - attr(*, "class")= chr [1:3] "fm_lattice_2d_list" "fm_list" "list"
```
