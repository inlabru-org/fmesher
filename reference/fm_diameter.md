# Diameter bound for a geometric object

Find an upper bound to the convex hull of a point set or function space

## Usage

``` r
fm_diameter(x, ...)

# S3 method for class 'matrix'
fm_diameter(x, manifold = NULL, ...)

# S3 method for class 'sf'
fm_diameter(x, ...)

# S3 method for class 'sfg'
fm_diameter(x, ...)

# S3 method for class 'sfc'
fm_diameter(x, ...)

# S3 method for class 'fm_lattice_2d'
fm_diameter(x, ...)

# S3 method for class 'fm_mesh_1d'
fm_diameter(x, ...)

# S3 method for class 'fm_mesh_2d'
fm_diameter(x, ...)

# S3 method for class 'fm_segm'
fm_diameter(x, ...)

# S3 method for class 'fm_mesh_3d'
fm_diameter(x, ...)

# S3 method for class 'fm_tensor'
fm_diameter(x, ..., multi = FALSE)

# S3 method for class 'fm_collect'
fm_diameter(x, ..., multi = FALSE)

# S3 method for class 'fm_list'
fm_diameter(x, ...)
```

## Arguments

- x:

  A point set as an \\n\times d\\ matrix, or an `fm_mesh_2d`/`1d`/`sf`
  related object.

- ...:

  Additional parameters passed on to the submethods.

- manifold:

  Character string specifying the manifold type. Default for `matrix`
  input is to treat the point set with Euclidean \\\mathbb{R}^d\\
  metrics. Use `manifold="S2"` for great circle distances on a sphere
  centred at the origin.

- multi:

  logical; For multi-domain spaces (e.g.
  [fm_tensor](https://inlabru-org.github.io/fmesher/reference/fm_tensor.md)
  and
  [fm_collect](https://inlabru-org.github.io/fmesher/reference/fm_collect.md)),
  if `TRUE`, return a vector of diameter bounds for each domain. If
  `FALSE` (the default), return a single diameter bound, by taking the
  maximum of the individual bounds.

## Value

A scalar, upper bound for the diameter of the convex hull of the point
set. For multi-domain spaces (e.g.
[`fm_tensor()`](https://inlabru-org.github.io/fmesher/reference/fm_tensor.md)
and
[`fm_collect()`](https://inlabru-org.github.io/fmesher/reference/fm_collect.md)),
a vector of upper bounds for each domain is returned.

## Methods (by class)

- `fm_diameter(fm_tensor)`: Returns either a single diameter bound
  (default), or a vector of sub-domain bounds; see the `multi` argument.

- `fm_diameter(fm_collect)`: Returns either a single diameter bound
  (default), or a vector of sub-domain bounds; see the `multi` argument.

## Author

Finn Lindgren <Finn.Lindgren@gmail.com>

## Examples

``` r

fm_diameter(matrix(c(0, 1, 1, 0, 0, 0, 1, 1), 4, 2))
#> [1] 1.414214
```
