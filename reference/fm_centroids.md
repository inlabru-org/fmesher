# Extract triangle centroids from an `fm_mesh_2d`

Computes the centroids of the triangles of an
[`fm_mesh_2d()`](https://inlabru-org.github.io/fmesher/reference/fm_mesh_2d.md)
object.

## Usage

``` r
fm_centroids(x, format = NULL)
```

## Arguments

- x:

  An `fm_mesh_2d` object.

- format:

  character; `"sf"`, `"df"`, `"sp"`

## Value

An `sf`, `data.frame`, or `SpatialPointsDataFrame` object, with the
vertex coordinates, and a `.triangle` column with the triangle indices.

## See also

[`fm_vertices()`](https://inlabru-org.github.io/fmesher/reference/fm_vertices.md)

## Author

Finn Lindgren <Finn.Lindgren@gmail.com>

## Examples

``` r
if (require("ggplot2", quietly = TRUE)) {
  vrt <- fm_centroids(fmexample$mesh, format = "sf")
  ggplot() +
    geom_sf(data = fm_as_sfc(fmexample$mesh)) +
    geom_sf(data = vrt, color = "red")
}

```
