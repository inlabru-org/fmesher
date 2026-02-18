# Internal helper functions for mesh field evaluation

Methods called internally by
[`fm_basis()`](https://inlabru-org.github.io/fmesher/reference/fm_basis.md)
methods.

## Usage

``` r
fm_basis_mesh_2d(
  mesh,
  loc = NULL,
  weights = NULL,
  derivatives = NULL,
  crs = NULL,
  ...
)

fm_basis_mesh_1d(mesh, loc, weights = NULL, derivatives = NULL, ...)
```

## Arguments

- loc:

  A location/value information object (`numeric`, `matrix`, `sf`,
  `fm_bary`, etc, depending on the class of `x`)

- weights:

  Optional weight vector, one weight for each location

- derivatives:

  logical; If true, also return matrices `dA` and `d2A` for `fm_mesh_1d`
  objects, and `dx`, `dy`, `dz` for `fm_mesh_2d`.

- ...:

  Passed on to submethods

## Value

A `fm_basis` object; a list of evaluator information objects, at least a
matrix `A` and logical vector `ok`.

## Examples

``` r
str(fm_basis_mesh_2d(fmexample$mesh, loc = fmexample$loc))
#> List of 3
#>  $ bary: fm_bary [10 × 2] (S3: fm_bary/tbl_df/tbl/data.frame)
#>   ..$ index: int [1:10] 70 47 393 83 363 326 231 395 395 13
#>   ..$ where: num [1:10, 1:3] 0.33809 0.08736 0.79683 0.03091 0.00355 ...
#>  $ A   :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
#>   .. ..@ i       : int [1:30] 3 9 3 6 0 5 4 2 3 0 ...
#>   .. ..@ p       : int [1:293] 0 0 0 0 0 0 0 0 0 0 ...
#>   .. ..@ Dim     : int [1:2] 10 292
#>   .. ..@ Dimnames:List of 2
#>   .. .. ..$ : NULL
#>   .. .. ..$ : NULL
#>   .. ..@ x       : num [1:30] 0.0309 0.1235 0.4649 0.5614 0.3381 ...
#>   .. ..@ factors : list()
#>  $ ok  : logi [1:10] TRUE TRUE TRUE TRUE TRUE TRUE ...
#>  - attr(*, "class")= chr "fm_basis"
```
