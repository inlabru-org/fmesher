# Check if a crs is NULL or NA

Methods of checking whether various kinds of CRS objects are `NULL` or
`NA`. Logically equivalent to either `is.na(fm_crs(x))` or
`is.na(fm_crs(x, oblique = NA))`, but with a short-cut pre-check for
`is.null(x)`.

## Usage

``` r
fm_crs_is_null(x, crsonly = FALSE)

# S3 method for class 'fm_crs'
is.na(x)
```

## Arguments

- x:

  An object supported by `fm_crs(x)`

- crsonly:

  For crs objects with extended functionality, such as
  [`fm_crs()`](https://inlabru-org.github.io/fmesher/reference/fm_crs.md)
  objects with `oblique` information, `crsonly = TRUE` only checks the
  plain CRS part.

## Value

logical

## Functions

- `fm_crs_is_null()`: Check if an object is or has `NULL` or `NA` CRS
  information. If not `NULL`, `is.na(fm_crs(x))` is returned. This
  allows the input to be e.g. a proj4string or epsg number, since the
  default
  [`fm_crs()`](https://inlabru-org.github.io/fmesher/reference/fm_crs.md)
  method passes its argument on to
  [`sf::st_crs()`](https://r-spatial.github.io/sf/reference/st_crs.html).

- `is.na(fm_crs)`: Check if a `fm_crs` has `NA` crs information and `NA`
  obliqueness

## See also

[`fm_crs()`](https://inlabru-org.github.io/fmesher/reference/fm_crs.md),
[fm_CRS()](https://inlabru-org.github.io/fmesher/reference/fm_crs.md),
[`fm_crs_is_identical()`](https://inlabru-org.github.io/fmesher/reference/fm_crs_is_identical.md)

## Examples

``` r
fm_crs_is_null(NULL)
#> [1] TRUE
fm_crs_is_null(27700)
#> [1] FALSE
fm_crs_is_null(fm_crs())
#> [1] TRUE
fm_crs_is_null(fm_crs(27700))
#> [1] FALSE
fm_crs_is_null(fm_crs(oblique = c(1, 2, 3, 4)))
#> [1] FALSE
fm_crs_is_null(fm_crs(oblique = c(1, 2, 3, 4)), crsonly = TRUE)
#> [1] TRUE
fm_crs_is_null(fm_crs(27700, oblique = c(1, 2, 3, 4)))
#> [1] FALSE
fm_crs_is_null(fm_crs(27700, oblique = c(1, 2, 3, 4)), crsonly = TRUE)
#> [1] FALSE
```
