# Print method for `fm_basis`

Prints information for an
[fm_basis](https://inlabru-org.github.io/fmesher/reference/fm_basis.md)
object.

## Usage

``` r
# S3 method for class 'fm_basis'
print(x, ..., prefix = "")
```

## Arguments

- x:

  [`fm_basis()`](https://inlabru-org.github.io/fmesher/reference/fm_basis.md)
  object

- ...:

  Unused

- prefix:

  a prefix to be used for each line. Default is an empty string.

## Value

`invisible(x)`

## See also

[`fm_basis()`](https://inlabru-org.github.io/fmesher/reference/fm_basis.md)

## Examples

``` r
print(fm_basis(fmexample$mesh, fmexample$loc, full = TRUE))
#> fm_basis object
#>   Projection matrix (A): 10-by-292
#>   Valid evaluations (ok): 10 out of 10
#>   Additional information: bary
```
