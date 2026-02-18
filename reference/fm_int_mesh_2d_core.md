# Integration scheme for mesh triangle interiors

Integration scheme for mesh triangle interiors

## Usage

``` r
fm_int_mesh_2d_core(mesh, tri_subset = NULL, nsub = NULL)
```

## Arguments

- mesh:

  Mesh on which to integrate

- tri_subset:

  Optional triangle index vector for integration on a subset of the mesh
  triangles (Default `NULL`)

- nsub:

  number of subdivision points along each triangle edge, giving
  `(nsub + 1)^2` proto-integration points used to compute the vertex
  weights (default `nsub=9`, giving 100 integration points for each
  triangle)

## Value

`tibble` with columns `loc`, `weight`, and `bary` with integration
points for the mesh

## Author

Finn Lindgren <Finn.Lindgren@gmail.com>

## Examples

``` r
str(fm_int_mesh_2d_core(fmexample$mesh))
#> tibble [54,700 × 3] (S3: tbl_df/tbl/data.frame)
#>  $ loc   : num [1:54700, 1:3] 0.0482 0.0562 0.0642 0.0722 0.0802 ...
#>  $ weight: num [1:54700] 0.000544 0.000544 0.000544 0.000544 0.000544 ...
#>  $ bary  : fm_bary [54,700 × 2] (S3: fm_bary/tbl_df/tbl/data.frame)
#>   ..$ index: int [1:54700] 1 1 1 1 1 1 1 1 1 1 ...
#>   ..$ where: num [1:54700, 1:3] 0.933 0.833 0.733 0.633 0.533 ...
```
