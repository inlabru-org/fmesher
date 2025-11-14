# Multi-domain sampler integration

Combine integration over different domains

## Usage

``` r
fm_int_multi_sampler(domain, samplers, ...)
```

## Arguments

- domain:

  A list of named domains

- samplers:

  A named list of samplers

- ...:

  Passed on to each
  [`fm_int()`](https://inlabru-org.github.io/fmesher/reference/fm_int.md)
  call.

## Value

An object with integration points and weights

## Examples

``` r
fm_int_multi_sampler(
  domain = list(x = fm_mesh_1d(1:4), y = 11:12),
  samplers = tibble::tibble(
    x = rbind(c(1, 3), c(2, 4)),
    y = c(12, 11)
  )
)
#> # A tibble: 10 × 5
#>        x     y .block weight .block_origin[,"x"] [,"y"]
#>    <dbl> <int>  <int>  <dbl>               <int>  <int>
#>  1   1      12      1  0.167                   1      1
#>  2   2      12      1  0.333                   1      1
#>  3   3      12      1  0.167                   1      1
#>  4   1.5    12      1  0.667                   1      1
#>  5   2.5    12      1  0.667                   1      1
#>  6   2      11      2  0.167                   2      2
#>  7   3      11      2  0.333                   2      2
#>  8   4      11      2  0.167                   2      2
#>  9   2.5    11      2  0.667                   2      2
#> 10   3.5    11      2  0.667                   2      2
```
