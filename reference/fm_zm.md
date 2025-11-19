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

# S3 method for class 'sfg'
fm_zm(x, ..., add = NULL, remove = NULL, target = NULL)

# S3 method for class 'numeric'
fm_zm(x, ..., add = NULL, remove = NULL, target = NULL, input = NULL)

# S3 method for class 'matrix'
fm_zm(x, ..., add = NULL, remove = NULL, target = NULL, input = NULL)

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

- `fm_zm_target()`: Determines the target target Z/M format

## Author

Finn Lindgren <Finn.Lindgren@gmail.com>
