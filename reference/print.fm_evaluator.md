# Print method for `fm_evaluator()`

Prints information for an
[fm_evaluator](https://inlabru-org.github.io/fmesher/reference/fm_evaluate.md)
object.

## Usage

``` r
# S3 method for class 'fm_evaluator'
print(x, ...)
```

## Arguments

- x:

  [`fm_evaluator()`](https://inlabru-org.github.io/fmesher/reference/fm_evaluate.md)
  object

- ...:

  Unused

## Value

`invisible(x)`

## See also

[`fm_evaluator()`](https://inlabru-org.github.io/fmesher/reference/fm_evaluate.md)

## Examples

``` r
print(fm_evaluator(fmexample$mesh, fmexample$loc))
#> fm_evaluator object
#>   proj:
#>     fm_basis object
#>       Projection matrix (A): 10-by-279
#>       Valid evaluations (ok): 10 out of 10
#>       Additional information: bary
#>   Additional evaluator information: x, y, lattice, loc, crs
```
