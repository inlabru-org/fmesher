# 3D tetrahedralisation storage

Internal C++ method.

(...)

## Usage

``` r
fmesher_mesh3d(options, loc, tv)
```

## Arguments

- options:

  list of triangulation options

- loc:

  numeric matrix; initial points to include

- tv:

  4-column integer matrix with 0-based vertex indices for each triangle

## Value

A list of information objects for a generated tetrahedralisation

## Examples

``` r
m <- fmesher_mesh3d(list(),
                    matrix(c(1,0,0,0,1,0,0,0,1,0,0,0), 4, 3, byrow=TRUE),
                    matrix(c(0,1,2,3), 1, 4, byrow=TRUE))
```
