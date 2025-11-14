# Query if points are inside a mesh

Queries whether each input point is within a mesh or not.

## Usage

``` r
fm_is_within(x, y, ...)
```

## Arguments

- x:

  A set of points/locations of a class supported by
  `fm_basis(y, loc = x, ..., full = TRUE)`

- y:

  An
  [fm_mesh_2d](https://inlabru-org.github.io/fmesher/reference/fm_mesh_2d.md)
  or other class supported by `fm_basis(y, loc = x, ..., full = TRUE)`

- ...:

  Passed on to
  [`fm_basis()`](https://inlabru-org.github.io/fmesher/reference/fm_basis.md)

## Value

A logical vector

## Examples

``` r
all(fm_is_within(fmexample$loc, fmexample$mesh))
#> [1] TRUE
```
