# Finite element matrix computation

Construct finite element structure matrices

## Usage

``` r
fmesher_fem(mesh_loc, mesh_tv, fem_order_max, aniso, options)
```

## Arguments

- mesh_loc:

  numeric matrix; mesh vertex coordinates

- mesh_tv:

  3-column integer matrix with 0-based vertex indices for each triangle

- fem_order_max:

  integer; the highest operator order to compute

- aniso:

  If non-NULL, a `list(gamma, v)`. Calculates anisotropic structure
  matrices (in addition to the regular) for \\\gamma\\ and \\v\\ for an
  anisotropic operator \\\nabla\cdot H \nabla\\, where \\H=\gamma I + v
  v^\top\\. Currently (2023-08-05) the fields need to be given per
  vertex.

- options:

  list of triangulation options (`sphere_tolerance`)

## Value

A list of matrices

## Examples

``` r
m <- fmesher_rcdt(list(cet_margin = 1), matrix(0, 1, 2))
b <- fmesher_fem(m$s, m$tv, fem_order_max = 2, aniso = NULL, options = list())
```
