# Barycentric coordinate computation

Internal C++ method.

Locate points and compute triangular barycentric coordinates

## Usage

``` r
fmesher_bary3d(mesh_loc, mesh_tv, loc, options)
```

## Arguments

- mesh_loc:

  numeric matrix; mesh vertex coordinates

- mesh_tv:

  3-column integer matrix with 0-based vertex indices for each triangle

- loc:

  numeric matrix; coordinates of points to locate in the mesh

- options:

  list of triangulation options

## Value

A list with vector `index` (tetra index) and matrix `where` (4-column
barycentric matrix)

## Examples

``` r
m <- fmesher_mesh3d(list(cet_margin = 1),
                    matrix(rnorm(15), 5, 3),
                    matrix(c(0,1,2,3), 1, 4))
b <- fmesher_bary3d(m$loc,
                    m$tv,
                    matrix(c(0.5, 0.5, 0.5), 1, 3),
                    list())
```
