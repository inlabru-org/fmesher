# Show expanded CRS arguments

**\[deprecated\]** Wrappers for
[`sp::CRS`](https://edzer.github.io/sp/reference/CRS-class.html) and
`inla.CRS` objects to handle the coordinate reference system argument
string. These methods should no longer be used with PROJ6/rgdal3; see
[`fm_wkt()`](https://inlabru-org.github.io/fmesher/reference/fm_crs_wkt.md)
and
[`fm_proj4string()`](https://inlabru-org.github.io/fmesher/reference/fm_crs_wkt.md)
for a new approach.

## Usage

``` r
fm_CRS_as_list(x, ...)

fm_list_as_CRS(x, ...)

fm_CRSargs(x, ...)

fm_list_as_CRSargs(x, ...)

fm_CRSargs_as_list(x, ...)
```

## Arguments

- x:

  An [`sp::CRS`](https://edzer.github.io/sp/reference/CRS-class.html) or
  `inla.CRS` object (for `fm_CRSargs` and `fm_CRS_as_list`), a character
  string (for `fm_CRSargs_as_list`), or a list (for `fm_list_as_CRS` and
  `fm_list_as_CRSargs`).

- ...:

  Additional arguments passed on to other methods.

## Value

For `fm_CRSargs` and `fm_list_as_CRSargs`, a character string with
PROJ.4 arguments.

For `fm_CRS_as_list` and `fm_CRSargs_as_list`, a list of name/value
pairs.

For `fm_list_as_CRS`, a `CRS` or `inla.CRS` object.

For `fm_list_as_CRSargs()`, a CRS proj4 string for name=value pair list

For `fm_CRSargs_as_list()`, a list of name=value pairs from CRS
proj4string

## See also

[`fm_CRS()`](https://inlabru-org.github.io/fmesher/reference/fm_CRS_sp.md)

## Author

Finn Lindgren <Finn.Lindgren@gmail.com>

## Examples

``` r
if (fm_safe_sp()) {
  crs0 <- fm_CRS("longlat_norm")
  p4s <- fm_proj4string(crs0)
  lst <- fm_CRSargs_as_list(p4s)
  crs1 <- fm_list_as_CRS(lst)
  lst$a <- 2
  crs2 <- fm_CRS(p4s, args = lst)
  print(fm_proj4string(crs0))
  print(fm_proj4string(crs1))
  print(fm_proj4string(crs2))
}
#> [1] "+proj=longlat +R=1 +no_defs"
#> [1] "+proj=longlat +R=1 +no_defs"
#> [1] "+proj=longlat +R=1 +no_defs"
```
