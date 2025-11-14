# Low level triangulation mesh plotting

Plots a triangulation mesh using `rgl`.

## Usage

``` r
plot_rgl(x, ...)

lines_rgl(x, ..., add = TRUE)

# S3 method for class 'fm_segm'
lines_rgl(
  x,
  loc = NULL,
  col = NULL,
  colors = c("black", "blue", "red", "green"),
  ...,
  add = TRUE
)

# S3 method for class 'fm_mesh_2d'
plot_rgl(
  x,
  col = "white",
  color.axis = NULL,
  color.n = 512,
  color.palette = cm.colors,
  color.truncate = FALSE,
  alpha = NULL,
  lwd = 1,
  specular = "black",
  draw.vertices = TRUE,
  draw.edges = TRUE,
  draw.faces = TRUE,
  draw.segments = draw.edges,
  size = 2,
  edge.color = rgb(0.3, 0.3, 0.3),
  t.sub = seq_len(nrow(x$graph$tv)),
  visibility = "",
  S = deprecated(),
  add = FALSE,
  ...
)

# S3 method for class 'fm_segm'
plot_rgl(x, ..., add = FALSE)

# S3 method for class 'fm_segm_list'
plot_rgl(x, ...)

# S3 method for class 'fm_segm_list'
lines_rgl(x, ...)
```

## Arguments

- x:

  A
  [`fm_mesh_2d()`](https://inlabru-org.github.io/fmesher/reference/fm_mesh_2d.md)
  object

- ...:

  Additional parameters passed to and from other methods.

- add:

  If `TRUE`, adds to the current plot instead of starting a new one.

- loc:

  Point locations to be used if `x$loc` is `NULL`.

- col:

  Segment color specification.

- colors:

  Colors to cycle through if `col` is `NULL`.

- color.axis:

  The min/max limit values for the color mapping.

- color.n:

  The number of colors to use in the color palette.

- color.palette:

  A color palette function.

- color.truncate:

  If `TRUE`, truncate the colors at the color axis limits.

- alpha:

  Transparency/opaqueness values. See `rgl.material`.

- lwd:

  Line width for edges. See `rgl.material`.

- specular:

  Specular color. See `rgl.material`.

- draw.vertices:

  If `TRUE`, draw triangle vertices.

- draw.edges:

  If `TRUE`, draw triangle edges.

- draw.faces:

  If `TRUE`, draw triangles.

- draw.segments:

  If `TRUE`, draw boundary and interior constraint edges more
  prominently.

- size:

  Size for vertex points.

- edge.color:

  Edge color specification.

- t.sub:

  Optional triangle index subset to be drawn.

- visibility:

  If "front" only display mesh faces with normal pointing towards the
  camera.

- S:

  Deprecated.

## Value

An rgl device identifier, invisibly.

## See also

[`plot.fm_mesh_2d()`](https://inlabru-org.github.io/fmesher/reference/plot.fm_mesh_2d.md)

## Author

Finn Lindgren <Finn.Lindgren@gmail.com>

## Examples

``` r
# \donttest{
if (interactive() && require("rgl")) {
  mesh <- fm_rcdt_2d(globe = 10)
  plot_rgl(mesh, col = mesh$loc[, 1])
}
# }
```
