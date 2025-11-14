# Assignment operators for crs information objects

Assigns new crs information.

## Usage

``` r
fm_crs(x) <- value

fm_crs_oblique(x) <- value

## S3 replacement method for class 'NULL'
fm_crs(x) <- value

## S3 replacement method for class 'NULL'
fm_crs_oblique(x) <- value

# S3 method for class 'fm_segm'
fm_crs(x) <- value

# S3 method for class 'fm_list'
fm_crs(x) <- value

# S3 method for class 'fm_mesh_2d'
fm_crs(x) <- value

# S3 method for class 'fm_collect'
fm_crs(x) <- value

# S3 method for class 'fm_lattice_2d'
fm_crs(x) <- value

# S3 method for class 'sf'
fm_crs(x) <- value

# S3 method for class 'sfg'
fm_crs(x) <- value

# S3 method for class 'sfc'
fm_crs(x) <- value

# S3 method for class 'Spatial'
fm_crs(x) <- value

# S3 method for class 'crs'
fm_crs_oblique(x) <- value

# S3 method for class 'CRS'
fm_crs_oblique(x) <- value

# S3 method for class 'fm_CRS'
fm_crs_oblique(x) <- value

# S3 method for class 'fm_crs'
fm_crs_oblique(x) <- value

# S3 method for class 'fm_segm'
fm_crs_oblique(x) <- value

# S3 method for class 'fm_mesh_2d'
fm_crs_oblique(x) <- value

# S3 method for class 'fm_collect'
fm_crs_oblique(x) <- value

# S3 method for class 'fm_lattice_2d'
fm_crs_oblique(x) <- value

# S3 method for class 'inla.CRS'
fm_crs_oblique(x) <- value
```

## Arguments

- x:

  Object to assign crs information to

- value:

  For `fm_crs<-()`, object supported by `fm_crs(value)`.

  For `fm_crs_oblique<-()`, `NA` or a numeric vector, see the `oblique`
  argument for
  [`fm_crs()`](https://inlabru-org.github.io/fmesher/reference/fm_crs.md).
  For assignment, `NULL` is treated as `NA`.

## Value

The modified object

## Functions

- `fm_crs(x) <- value`: Automatically converts the input value with
  `fm_crs(value)`, `fm_crs(value, oblique = NA)`, `fm_CRS(value)`, or
  `fm_CRS(value, oblique = NA)`, depending on the type of `x`.

- `fm_crs_oblique(x) <- value`: Assigns new `oblique` information.

## See also

[`fm_crs()`](https://inlabru-org.github.io/fmesher/reference/fm_crs.md)

## Examples

``` r
x <- fm_segm()
fm_crs(x) <- fm_crs("+proj=longlat")
fm_crs(x)$proj4string
#> [1] "+proj=longlat +datum=WGS84 +no_defs"
```
