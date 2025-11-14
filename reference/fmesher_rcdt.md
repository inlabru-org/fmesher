# Refined Constrained Delaunay Triangulation

(...)

## Usage

``` r
fmesher_rcdt(
  options,
  loc,
  tv = NULL,
  boundary = NULL,
  interior = NULL,
  boundary_grp = NULL,
  interior_grp = NULL
)
```

## Arguments

- options:

  list of triangulation options

- loc:

  numeric matrix; initial points to include

- tv:

  3-column integer matrix with 0-based vertex indices for each triangle

- boundary:

  2-column integer matrix with 0-based vertex indices for each boundary
  edge constraint

- interior:

  2-column integer matrix with 0-based vertex indices for each interior
  edge constraint

- boundary_grp:

  integer vector with group labels

- interior_grp:

  integer vector with group labels

## Value

A list of information objects for a generated triangulation

## Examples

``` r
m <- fmesher_rcdt(list(cet_margin = 1), matrix(0, 1, 2))
```
