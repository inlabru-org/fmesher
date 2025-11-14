# Deprecated functions in fmesher

These functions still attempt to do their job, but will be removed in a
future version.

## Usage

``` r
fm_spTransform(x, ...)

# Default S3 method
fm_spTransform(x, crs0 = NULL, crs1 = NULL, passthrough = FALSE, ...)

# S3 method for class 'SpatialPoints'
fm_spTransform(x, CRSobj, passthrough = FALSE, ...)

# S3 method for class 'SpatialPointsDataFrame'
fm_spTransform(x, CRSobj, passthrough = FALSE, ...)

fm_sp2segment(...)
```

## Arguments

- x:

  The object that should be transformed from it's current CRS to a new
  CRS

- ...:

  Potential additional arguments

- crs0:

  The source sp::CRS or inla.CRS object

- crs1:

  The target sp::CRS or inla.CRS object

- passthrough:

  Default is FALSE. Setting to TRUE allows objects with no CRS
  information to be passed through without transformation.

- CRSobj:

  The target sp::CRS or inla.CRS object

## Functions

- `fm_spTransform()`: **\[deprecated\]** (See
  [`fm_transform()`](https://inlabru-org.github.io/fmesher/reference/fm_transform.md)
  instead) Handle transformation of various inla objects according to
  coordinate reference systems of
  [`sp::CRS`](https://edzer.github.io/sp/reference/CRS-class.html) or
  `INLA::inla.CRS` class.

- `fm_spTransform(default)`: The default method handles low level
  transformation of raw coordinates.

- `fm_sp2segment()`: **\[deprecated\]** in favour of
  [`fm_as_segm()`](https://inlabru-org.github.io/fmesher/reference/fm_as_segm.md)

## See also

[`fm_transform()`](https://inlabru-org.github.io/fmesher/reference/fm_transform.md)

## Author

Finn Lindgren <Finn.Lindgren@gmail.com>
