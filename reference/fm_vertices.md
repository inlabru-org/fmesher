# Extract vertex locations from an `fm_mesh_2d`

Extracts the vertices of an `fm_mesh_2d` object.

## Usage

``` r
fm_vertices(x, format = NULL)
```

## Arguments

- x:

  An `fm_mesh_2d` object.

- format:

  character; `"sf"`, `"df"`, `"sp"`

## Value

An `sf`, `data.frame`, or `SpatialPointsDataFrame` object, with the
vertex coordinates, and a `.vertex` column with the vertex indices.

## See also

[`fm_centroids()`](https://inlabru-org.github.io/fmesher/reference/fm_centroids.md)

## Author

Finn Lindgren <Finn.Lindgren@gmail.com>

## Examples

``` r
if (require("ggplot2", quietly = TRUE)) {
  vrt <- fm_vertices(fmexample$mesh, format = "sf")
  ggplot() +
    geom_sf(data = fm_as_sfc(fmexample$mesh)) +
    geom_sf(data = vrt, color = "red")
}

```
