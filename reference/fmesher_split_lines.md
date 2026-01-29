# Split lines at triangle edges

Internal C++ method.

Split a sequence of line segments at triangle edges

## Usage

``` r
fmesher_split_lines(mesh_loc, mesh_tv, loc, idx, options)
```

## Arguments

- mesh_loc:

  numeric matrix; mesh vertex coordinates

- mesh_tv:

  3-column integer matrix with 0-based vertex indices for each triangle

- loc:

  numeric coordinate matrix

- idx:

  2-column integer matrix

- options:

  list of triangulation options (`sphere_tolerance`)

## Value

A list of line splitting information objects

## See also

[`fm_split_lines()`](https://inlabru-org.github.io/fmesher/reference/fm_split_lines.md)

## Examples

``` r
mesh <- fm_mesh_2d(
  boundary = fm_segm(rbind(c(0,0), c(1,0), c(1,1), c(0, 1)), is.bnd = TRUE)
)
splitter <- fm_segm(rbind(c(0.8, 0.2), c(0.2, 0.8)))
segm_split <- fm_split_lines(mesh, splitter)
```
