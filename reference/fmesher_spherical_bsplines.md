# Rotationally invariant spherical B-splines

Internal C++ method.

Compute rotationally invariant spherical B-splines on the unit sphere

## Usage

``` r
fmesher_spherical_bsplines1(loc, n, degree, uniform)

fmesher_spherical_bsplines(loc, n, degree, uniform)
```

## Arguments

- loc:

  numeric vector/matrix; coordinates of points to locate in the mesh,
  only the z-coordinates are used (`sin(latitude)`)

- n:

  The number of basis functions

- degree:

  The polynomial basis degree

- uniform:

  logical; If `TRUE`, the knots are spaced uniformly by latitude, if
  `FALSE`, the knots are spaced uniformly by `sin(latitude)`

## Value

A matrix of evaluated b-spline basis functions

## Examples

``` r
m <- fm_rcdt_2d(globe = 1)
fmesher_spherical_bsplines(m$loc, n = 3, degree = 2, uniform = FALSE)
#>         [,1]  [,2]   [,3]
#>  [1,] 0.0000 0.000 1.0000
#>  [2,] 0.0625 0.375 0.5625
#>  [3,] 0.0625 0.375 0.5625
#>  [4,] 0.0625 0.375 0.5625
#>  [5,] 0.0625 0.375 0.5625
#>  [6,] 0.0625 0.375 0.5625
#>  [7,] 0.5625 0.375 0.0625
#>  [8,] 0.5625 0.375 0.0625
#>  [9,] 0.5625 0.375 0.0625
#> [10,] 0.5625 0.375 0.0625
#> [11,] 0.5625 0.375 0.0625
#> [12,] 1.0000 0.000 0.0000
fmesher_spherical_bsplines1(m$loc[, 3], n = 3, degree = 2, uniform = FALSE)
#>         [,1]  [,2]   [,3]
#>  [1,] 0.0000 0.000 1.0000
#>  [2,] 0.0625 0.375 0.5625
#>  [3,] 0.0625 0.375 0.5625
#>  [4,] 0.0625 0.375 0.5625
#>  [5,] 0.0625 0.375 0.5625
#>  [6,] 0.0625 0.375 0.5625
#>  [7,] 0.5625 0.375 0.0625
#>  [8,] 0.5625 0.375 0.0625
#>  [9,] 0.5625 0.375 0.0625
#> [10,] 0.5625 0.375 0.0625
#> [11,] 0.5625 0.375 0.0625
#> [12,] 1.0000 0.000 0.0000
```
