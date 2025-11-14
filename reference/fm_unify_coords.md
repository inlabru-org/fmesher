# Unify coordinates to 3-column matrix

Convert coordinate information to a 3-column matrix. This is mainly an
internal function, and the interface may change.

## Usage

``` r
fm_unify_coords(x, crs = NULL)

## S3 method for class 'NULL'
fm_unify_coords(x, crs = NULL)

# Default S3 method
fm_unify_coords(x, crs = NULL)

# S3 method for class 'Spatial'
fm_unify_coords(x, crs = NULL)

# S3 method for class 'sf'
fm_unify_coords(x, crs = NULL)

# S3 method for class 'sfc'
fm_unify_coords(x, crs = NULL)
```

## Arguments

- x:

  A object with coordinate information

- crs:

  A optional crs object to convert the coordinates to

## Value

A coordinate matrix

## Examples

``` r
fm_unify_coords(fmexample$loc_sf)
#>             [,1]        [,2] [,3]
#>  [1,] -1.2070657 -0.47719270    0
#>  [2,]  0.2774292 -0.99838644    0
#>  [3,]  1.0844412 -0.77625389    0
#>  [4,] -2.3456977  0.06445882    0
#>  [5,]  0.4291247  0.95949406    0
#>  [6,]  0.5060559 -0.11028549    0
#>  [7,] -0.5747400 -0.51100951    0
#>  [8,] -0.5466319 -0.91119542    0
#>  [9,] -0.5644520 -0.83717168    0
#> [10,] -0.8900378  2.41583518    0
```
