# Compute finite element matrices

Compute finite element mass and structure matrices

## Usage

``` r
fm_fem(mesh, order = 2, ...)

# S3 method for class 'fm_mesh_1d'
fm_fem(mesh, order = 2, ...)

# S3 method for class 'fm_mesh_2d'
fm_fem(mesh, order = 2, aniso = NULL, ...)

# S3 method for class 'fm_tensor'
fm_fem(mesh, order = 2, ...)

# S3 method for class 'fm_collect'
fm_fem(mesh, order = 2, ...)

# S3 method for class 'fm_mesh_3d'
fm_fem(mesh, order = 2, ...)
```

## Arguments

- mesh:

  [`fm_mesh_1d()`](https://inlabru-org.github.io/fmesher/reference/fm_mesh_1d.md),
  [`fm_mesh_2d()`](https://inlabru-org.github.io/fmesher/reference/fm_mesh_2d.md),
  or other supported mesh class object

- order:

  integer; the maximum operator order

- ...:

  Currently unused

- aniso:

  If non-NULL, a `list(gamma, v)`. Calculates anisotropic structure
  matrices (in addition to the regular) for \\\gamma\\ and \\v\\ for an
  anisotropic operator \\\nabla\cdot H \nabla\\, where \\H=\gamma I + v
  v^\top\\. Currently (2023-08-05) the fields need to be given per
  vertex.

## Value

`fm_fem.fm_mesh_1d`: A list with elements `c0`, `c1`, `g1`, `g2`, etc.
When `mesh$degree == 2`, also `g01`, `g02`, and `g12`.

`fm_fem.fm_mesh_2d`: A list with elements `c0`, `c1`, `g1`, `va`, `ta`,
and more if `order > 1`. When `aniso` is non-NULL, also `g1aniso`
matrices, etc.

`fm_fem.fm_tensor`: A list with elements `cc`, `g1`, `g2`.

`fm_fem.fm_collect`: A list with elements `c0`, `c1`, `g1`, `g2`, etc,
and `cc` (`c0` for every model except `fm_mesh_1d` with `degree=2`, for
which it is `c1`). If the base type for the collection provides `va` and
`ta` values, those are also returned.

`fm_fem.fm_mesh_3d`: A list with elements `c0`, `c1`, `g1`, `g2`, `va`,
`ta`, and more if `order > 2`.

## Examples

``` r
names(fm_fem(fm_mesh_1d(1:4), order = 3))
#> [1] "c0" "c1" "g1" "g2" "g3"
names(fm_fem(fmexample$mesh, order = 3))
#>  [1] "b1" "c0" "c1" "g1" "g2" "g3" "k1" "k2" "k3" "ta" "va"
```
