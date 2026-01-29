# Globe points

C++ method, may get a stable R interface as `fm_globe_points()` in the
future.

Create points on a globe

## Usage

``` r
fmesher_globe_points(globe)
```

## Arguments

- globe:

  integer; the number of edge subdivision segments, 1 or higher.

## Value

A matrix of points on a unit radius globe

## Examples

``` r
fmesher_globe_points(1)
#>             [,1]          [,2] [,3]
#>  [1,]  0.0000000  0.000000e+00  1.0
#>  [2,]  0.8660254  0.000000e+00  0.5
#>  [3,]  0.2676166  8.236391e-01  0.5
#>  [4,] -0.7006293  5.090370e-01  0.5
#>  [5,] -0.7006293 -5.090370e-01  0.5
#>  [6,]  0.2676166 -8.236391e-01  0.5
#>  [7,]  0.7006293  5.090370e-01 -0.5
#>  [8,] -0.2676166  8.236391e-01 -0.5
#>  [9,] -0.8660254  1.060575e-16 -0.5
#> [10,] -0.2676166 -8.236391e-01 -0.5
#> [11,]  0.7006293 -5.090370e-01 -0.5
#> [12,]  0.0000000  0.000000e+00 -1.0
```
