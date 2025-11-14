# Detect manifold type

Detect if a 2d object is on "R2", "S2", or "M2"

## Usage

``` r
fm_detect_manifold(x)

fm_crs_detect_manifold(x)

# S3 method for class 'crs'
fm_detect_manifold(x)

# S3 method for class 'CRS'
fm_detect_manifold(x)

# S3 method for class 'numeric'
fm_detect_manifold(x)

# S3 method for class 'matrix'
fm_detect_manifold(x)

# S3 method for class 'fm_mesh_2d'
fm_detect_manifold(x)
```

## Arguments

- x:

  Object to investigate

## Value

A string containing the detected manifold classification

## Functions

- `fm_crs_detect_manifold()`: Detect if a crs is on "R2" or "S2" (if
  `fm_crs_is_geocent(crs)` is `TRUE`). Returns `NA_character_` if the
  crs is NULL or NA.

## Examples

``` r
fm_detect_manifold(1:4)
#> [1] "R1"
fm_detect_manifold(rbind(c(1, 0, 0), c(0, 1, 0), c(1, 1, 0)))
#> [1] "R2"
fm_detect_manifold(rbind(c(1, 0, 0), c(0, 1, 0), c(0, 0, 1)))
#> [1] "S2"
```
