# Compute areas and edge lengths

Internal C++ method.

Compute triangle areas, edge lengths, and vertex-associated areas

## Usage

``` r
fmesher_sizes_mesh2d(mesh_loc, mesh_tv, options)
```

## Arguments

- mesh_loc:

  numeric matrix; mesh vertex coordinates

- mesh_tv:

  3-column integer matrix with 0-based vertex indices for each triangle

- options:

  list of triangulation options (`sphere_tolerance`)

## Value

A list of `face`, `face_edge`, and `vertex`

## See also

[`fm_sizes()`](https://inlabru-org.github.io/fmesher/reference/fm_sizes.md)

## Examples

``` r
mesh <- fm_mesh_2d(
  boundary = fm_segm(rbind(c(0,0), c(1,0), c(1,1), c(0, 1)), is.bnd = TRUE)
)
sz <- fm_sizes(mesh)
summary(sz$face)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>     0.5     0.5     0.5     0.5     0.5     0.5 
```
