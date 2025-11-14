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
#>   ..$ index: int [1:10] 358 301 413 337 369 221 363 329 329 142
#>   ..$ where: num [1:10, 1:3] 0.0699 0.1367 0.4852 0.0831 0.6245 ...
#>  $ A   :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
#>   .. ..@ i       : int [1:30] 1 9 3 7 8 5 0 6 4 2 ...
#>   .. ..@ p       : int [1:280] 0 0 0 0 0 0 0 0 0 0 ...
#>   .. ..@ Dim     : int [1:2] 10 279
#>   .. ..@ Dimnames:List of 2
#>   .. .. ..$ : NULL
#>   .. .. ..$ : NULL
#>   .. ..@ x       : num [1:30] 0.3065 0.0477 0.0831 0.6947 0.4188 ...
#>   .. ..@ factors : list()
#>  $ ok  : logi [1:10] TRUE TRUE TRUE TRUE TRUE TRUE ...
#>  - attr(*, "class")= chr "fm_basis"
```
