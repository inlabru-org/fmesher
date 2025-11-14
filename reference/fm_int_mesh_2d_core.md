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

`tibble` with columns `loc` and `weight` with integration points for the
mesh

## Author

Finn Lindgren <Finn.Lindgren@gmail.com>

## Examples

``` r
str(fm_int_mesh_2d_core(fmexample$mesh))
#> tibble [52,700 × 2] (S3: tbl_df/tbl/data.frame)
#>  $ loc   : num [1:52700, 1:3] 1.58 1.57 1.56 1.55 1.53 ...
#>  $ weight: num [1:52700] 0.000641 0.000641 0.000641 0.000641 0.000641 ...
```
