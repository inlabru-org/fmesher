# Split lines at triangle edges

Compute intersections between line segments and triangle edges, and
filter out segment of length zero.

## Usage

``` r
fm_split_lines(mesh, ...)

# S3 method for class 'fm_mesh_2d'
fm_split_lines(mesh, segm, ...)
```

## Arguments

- mesh:

  An
  [fm_mesh_2d](https://inlabru-org.github.io/fmesher/reference/fm_mesh_2d.md)
  object

- ...:

  Unused.

- segm:

  An
  [`fm_segm()`](https://inlabru-org.github.io/fmesher/reference/fm_segm.md)
  object with segments to be split

## Value

An
[`fm_segm()`](https://inlabru-org.github.io/fmesher/reference/fm_segm.md)
object with the same crs as the mesh, with an added field `origin`, that
for each new segment gives the originator index into to original `segm`
object for each new line segment.

## Author

Finn Lindgren <Finn.Lindgren@gmail.com>

## Examples

``` r
mesh <- fm_mesh_2d(
  boundary = fm_segm(
    rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1)),
    is.bnd = TRUE
  )
)
splitter <- fm_segm(rbind(c(0.8, 0.2), c(0.2, 0.8)))
segm_split <- fm_split_lines(mesh, splitter)

plot(mesh)
lines(splitter)
points(segm_split$loc)
```
