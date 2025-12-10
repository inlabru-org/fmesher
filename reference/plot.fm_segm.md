# Draw `fm_segm` objects.

Draws a
[`fm_segm()`](https://inlabru-org.github.io/fmesher/reference/fm_segm.md)
object with generic or `rgl` graphics.

## Usage

``` r
# S3 method for class 'fm_segm'
plot(x, ..., add = FALSE)

# S3 method for class 'fm_segm'
lines(
  x,
  loc = NULL,
  col = NULL,
  colors = c("black", "blue", "red", "green"),
  add = TRUE,
  xlim = NULL,
  ylim = NULL,
  asp = 1,
  axes = FALSE,
  xlab = "",
  ylab = "",
  visibility = "front",
  rgl = deprecated(),
  ...
)

# S3 method for class 'fm_segm_list'
plot(x, ...)

# S3 method for class 'fm_segm_list'
lines(x, ...)
```

## Arguments

- x:

  An
  [`fm_segm()`](https://inlabru-org.github.io/fmesher/reference/fm_segm.md)
  object.

- ...:

  Additional parameters, passed on to graphics methods.

- add:

  If `TRUE`, add to the current plot, otherwise start a new plot.

- loc:

  Point locations to be used if `x$loc` is `NULL`.

- col:

  Segment color specification.

- colors:

  Colors to cycle through if `col` is `NULL`.

- xlim, ylim:

  X and Y axis limits for a new plot.

- asp:

  Aspect ratio for new plots. Default 1.

- axes:

  logical; whether axes should be drawn on the plot. Default FALSE.

- xlab, ylab:

  character; labels for the axes.

- visibility:

  If "front" only display mesh faces with normal pointing towards the
  camera.

- rgl:

  **\[deprecated\]** since `0.5.0.9000` in favour of the
  [`plot_rgl()`](https://inlabru-org.github.io/fmesher/reference/plot_rgl.md)
  and
  [`lines_rgl()`](https://inlabru-org.github.io/fmesher/reference/plot_rgl.md)
  methods. If `TRUE`, use `rgl` for plotting.

## Value

None

## See also

[`fm_segm()`](https://inlabru-org.github.io/fmesher/reference/fm_segm.md),
[plot.fm_mesh_2d](https://inlabru-org.github.io/fmesher/reference/plot.fm_mesh_2d.md)

## Author

Finn Lindgren <Finn.Lindgren@gmail.com>

## Examples

``` r
plot(fm_segm(fmexample$mesh, boundary = TRUE))
lines(fm_segm(fmexample$mesh, boundary = FALSE), col = 2)

```
