# Convert a 3D mesh to a 3D rgl triangulation

Extracts a matrix of coordinates of triangles, suitable for passing to
[`rgl::triangles3d()`](https://dmurdoch.github.io/rgl/dev/reference/primitives.html).

## Usage

``` r
as.triangles3d.fm_mesh_3d(obj, subset = NULL, ...)
```

## Arguments

- obj:

  An `fm_mesh_3d` object

- subset:

  Character string specifying which triangles to extract. Either "all"
  (default) or "boundary".

- ...:

  Currently unused

## Value

A 3-column matrix of coordinates of triangles, suitable for passing to
[`rgl::triangles3d()`](https://dmurdoch.github.io/rgl/dev/reference/primitives.html).

## Examples

``` r
# Protect against unavailable rgl device by only running interactively
if (interactive() &&
  requireNamespace("geometry", quietly = TRUE) &&
  requireNamespace("rgl", quietly = TRUE)) {
  (m <- fm_delaunay_3d(matrix(rnorm(30), 10, 3)))
  rgl::open3d()
  rgl::triangles3d(rgl::as.triangles3d(m, "boundary"), col = "blue")
  rgl::axes3d()
}
```
