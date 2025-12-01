# Spherical harmonics

Compute spherical harmonics on the unit sphere

## Usage

``` r
fmesher_spherical_harmonics(loc, max_order, rotationally_symmetric)
```

## Arguments

- loc:

  numeric matrix; coordinates of points to locate in the mesh

- rotationally_symmetric:

  logical; If `TRUE`, only evaluate rotationally invariant basis
  functions

- n:

  integer; the maximum basis order

## Value

A matrix of evaluated spherical harmonic basis functions

## Examples

``` r
m <- fm_rcdt_2d(globe = 1)
fmesher_spherical_bsplines(m$loc, max_order = 2, TRUE)
#> Error in fmesher_spherical_bsplines(m$loc, max_order = 2, TRUE): unused argument (max_order = 2)
fmesher_spherical_bsplines(m$loc, max_order = 2, FALSE)
#> Error in fmesher_spherical_bsplines(m$loc, max_order = 2, FALSE): unused argument (max_order = 2)
```
