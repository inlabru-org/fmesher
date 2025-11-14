# Check which mesh triangles are inside a polygon

Wrapper for the
[`sf::st_contains()`](https://r-spatial.github.io/sf/reference/geos_binary_pred.html)
(previously
[`sp::over()`](https://edzer.github.io/sp/reference/over.html)) method
to find triangle centroids or vertices inside `sf` or `sp` polygon
objects

## Usage

``` r
fm_contains(x, y, ...)

# S3 method for class 'Spatial'
fm_contains(x, y, ...)

# S3 method for class 'sf'
fm_contains(x, y, ...)

# S3 method for class 'sfc'
fm_contains(x, y, ..., type = c("centroid", "vertex"))
```

## Arguments

- x:

  geometry (typically an `sf` or
  [`sp::SpatialPolygons`](https://edzer.github.io/sp/reference/SpatialPolygons.html)
  object) for the queries

- y:

  an
  [`fm_mesh_2d()`](https://inlabru-org.github.io/fmesher/reference/fm_mesh_2d.md)
  object

- ...:

  Passed on to other methods

- type:

  the query type; either `'centroid'` (default, for triangle centroids),
  or `'vertex'` (for mesh vertices)

## Value

List of vectors of triangle indices (when `type` is `'centroid'`) or
vertex indices (when `type` is `'vertex'`). The list has one entry per
row of the `sf` object. Use `unlist(fm_contains(...))` if the combined
union is needed.

## Author

Haakon Bakka, <bakka@r-inla.org>, and Finn Lindgren
<Finn.Lindgren@gmail.com>

## Examples

``` r
# Create a polygon and a mesh
obj <- sf::st_sfc(
  sf::st_polygon(
    list(rbind(
      c(0, 0),
      c(50, 0),
      c(50, 50),
      c(0, 50),
      c(0, 0)
    ))
  ),
  crs = fm_crs("longlat_globe")
)
mesh <- fm_rcdt_2d_inla(globe = 2, crs = fm_crs("sphere"))

## 2 vertices found in the polygon
fm_contains(obj, mesh, type = "vertex")
#> Sparse geometry binary predicate list of length 1, where the predicate
#> was `contains'
#>  1: 8, 17

## 3 triangles found in the polygon
fm_contains(obj, mesh)
#> Sparse geometry binary predicate list of length 1, where the predicate
#> was `contains'
#>  1: 15, 18, 35

## Multiple transformations can lead to slightly different results
## due to edge cases:
## 4 triangles found in the polygon
fm_contains(
  obj,
  fm_transform(mesh, crs = fm_crs("mollweide_norm"))
)
#> Sparse geometry binary predicate list of length 1, where the predicate
#> was `contains'
#>  1: 7, 15, 18, 35
```
