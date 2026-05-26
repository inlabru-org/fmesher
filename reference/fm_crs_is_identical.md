# Check if two CRS objects are identical

Check if two CRS objects are identical

## Usage

``` r
fm_crs_is_identical(crs0, crs1, crsonly = FALSE)
```

## Arguments

- crs0, crs1:

  Two
  [`sf::crs`](https://r-spatial.github.io/sf/reference/coerce-methods.html),
  [`sp::CRS`](https://edzer.github.io/sp/reference/CRS-class.html),
  `fm_crs` or `inla.CRS` objects to be compared.

- crsonly:

  logical. If `TRUE` and any of `crs0` and `crs1` are `fm_crs` or
  `inla.CRS` objects, extract and compare only the
  [`sf::crs`](https://r-spatial.github.io/sf/reference/coerce-methods.html)
  or [`sp::CRS`](https://edzer.github.io/sp/reference/CRS-class.html)
  aspects. Default: `FALSE`

## Value

logical, indicating if the two crs objects are identical in the
specified sense (see the `crsonly` argument)

## See also

[`fm_crs()`](https://inlabru-org.github.io/fmesher/reference/fm_crs.md),
[fm_CRS()](https://inlabru-org.github.io/fmesher/reference/fm_crs.md),
[`fm_crs_is_null()`](https://inlabru-org.github.io/fmesher/reference/fm_crs_is_null.md)

## Examples

``` r

crs0 <- crs1 <- fm_crs("longlat_globe")
fm_crs_oblique(crs1) <- c(0, 90)
print(c(
  fm_crs_is_identical(crs0, crs0),
  fm_crs_is_identical(crs0, crs1),
  fm_crs_is_identical(crs0, crs1, crsonly = TRUE)
))
#> [1]  TRUE FALSE  TRUE
```
