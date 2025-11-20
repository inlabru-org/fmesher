# Add or remove Z/M information

**\[experimental\]** Add and/or remove Z and/or M information from
simple feature geometries.

## Usage

``` r
fm_zm(x, ...)

# S3 method for class 'sf'
fm_zm(x, ...)

# S3 method for class 'sfc'
fm_zm(x, ..., add = NULL, remove = NULL, target = NULL)

# S3 method for class 'list'
fm_zm(x, ..., add = NULL, remove = NULL, target = NULL)

# S3 method for class 'sfg'
fm_zm(x, ..., add = NULL, remove = NULL, target = NULL)

# S3 method for class 'numeric'
fm_zm(x, ..., add = NULL, remove = NULL, target = NULL, input = NULL)

# S3 method for class 'matrix'
fm_zm(x, ..., add = NULL, remove = NULL, target = NULL, input = NULL)

fm_zm_input(x, ...)

# S3 method for class 'sf'
fm_zm_input(x, ...)

# S3 method for class 'sfc'
fm_zm_input(x, ...)

# S3 method for class 'list'
fm_zm_input(x, ...)

# S3 method for class 'sfg'
fm_zm_input(x, ...)

# S3 method for class 'numeric'
fm_zm_input(x, ..., input = NULL)

# S3 method for class 'matrix'
fm_zm_input(x, ..., input = NULL)

fm_zm_target(input, add = NULL, remove = NULL, target = NULL)
```

## Arguments

- x:

  An object to modify

- ...:

  Further arguments passed to methods

- add:

  character; one of `NULL`, `"Z"`, `"M"`, or `"ZM"`. Specifies which
  dimensions to add.

- remove:

  character; one of `NULL`, `"Z"`, `"M"`, or `"ZM"`. Specifies which
  dimensions to remove.

- target:

  character; one of `"XY"`, `"XYZ"`, `"XYM"`, or `"XYZM"`. Specifies the
  target dimension format. If provided, overrides `add` and `remove`.
  When both `add` and `remove` are `NULL`, the default target is the
  smallest format that can hold all the inputs without loss of
  information.

- input:

  character or character vector; one of `NULL`, `"XY"`, `"XYZ"`,
  `"XYM"`, or `"XYZM"`. Specifies the input dimension format. If `NULL`
  (default), the input format is inferred from the number of columns in
  `x` (for matrices/numerics) or from the geometry type (for `sfc`
  objects).

## Value

An object of the same class as `x`, with modified Z/M dimensions.

## Functions

- `fm_zm_input()`: Find the set of distinct XY/XYZ/XYM/XYZM types

- `fm_zm_target()`: Determines the target XY/XYZ/XYM/XYZM format

## See also

[`sf::st_zm()`](https://r-spatial.github.io/sf/reference/st_zm.html)
that supports a subset of these operations.

## Author

Finn Lindgren <Finn.Lindgren@gmail.com>

## Examples

``` r
fm_zm(fmexample$loc_sf, add = "Z")
#> Geometry set for 10 features 
#> Geometry type: POINT
#> Dimension:     XYZ
#> Bounding box:  xmin: -2.345698 ymin: -0.9983864 xmax: 1.084441 ymax: 2.415835
#> z_range:       zmin: 0 zmax: 0
#> CRS:           NA
#> First 5 geometries:
#> POINT Z (-1.207066 -0.4771927 0)
#> POINT Z (0.2774292 -0.9983864 0)
#> POINT Z (1.084441 -0.7762539 0)
#> POINT Z (-2.345698 0.06445882 0)
#> POINT Z (0.4291247 0.9594941 0)

fm_zm_input(fmexample$loc_sf)
#> [1] "XY"

fm_zm_target(c("XY", "XYZ"))
#> [1] "XYZ"
fm_zm_target("XY", add = "Z")
#> [1] "XYZ"
fm_zm_target(c("XY", "XYZM"), remove = "M")
#> [1] "XYZ"
```
