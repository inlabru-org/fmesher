# Construct integration scheme objects

Constructor method for integration scheme objects, allowing default
construction of `.block` information. Primarily meant for internal use,
but can be used to manually create data of the same structure as
[`fm_int()`](https://inlabru-org.github.io/fmesher/reference/fm_int.md)
output.

## Usage

``` r
fm_int_object(
  object,
  blocks = FALSE,
  weight = NULL,
  name = NULL,
  override = FALSE
)
```

## Arguments

- object:

  An object representing integration points; either a data.frame-like
  object, or a vector/list of coordinates or other location reference
  objects.

- blocks:

  logical; if `TRUE`, set per-element `.block` indices. If `FALSE`
  (default), set a common block, `1L`.

- weight:

  Optional weight variable; if `NULL`, all weights are set to 1.

- name:

  character; name of the integration domain.

- override:

  logical; If `name` is non-NULL and `override=TRUE` for sf object, the
  current `sf_column` is renamed to `name`.

## Value

A tibble or sf/tibble object

## See also

[`fm_int()`](https://inlabru-org.github.io/fmesher/reference/fm_int.md)

## Examples

``` r
fm_int_object(1:4, blocks = TRUE, weight = c(1, 2, 1, 3), name = "z")
#> # A tibble: 4 × 4
#>       z weight .block .block_origin[,"z"]
#>   <int>  <dbl>  <int>               <int>
#> 1     1      1      1                   1
#> 2     2      2      2                   2
#> 3     3      1      3                   3
#> 4     4      3      4                   4
```
