# fm_sizes

**\[experimental\]** Compute effective sizes of faces/cells and vertices
in a mesh

## Usage

``` r
fm_sizes(...)

# S3 method for class 'fm_mesh_2d'
fm_sizes(mesh, ...)

# S3 method for class 'fm_mesh_3d'
fm_sizes(mesh, ...)
```

## Arguments

- ...:

  Passed on to submethods

- mesh:

  object of a supported mesh class

## Value

A `list` with elements `face` and `vertex` for 2D meshes, or `cell` and
`vertex` for 3D meshes. The elements are vectors of effective sizes of
the faces/cells and vertices, respectively. For 2D meshes, also
`face_edge`, a matrix with one row per triangle and 3 columns, with edge
lengths for the edge opposing each triangle vertex.

## Examples

``` r
str(fm_sizes(fmexample$mesh))
#> List of 3
#>  $ face     : num [1:527] 0.0641 0.0799 0.0794 0.0734 0.0383 ...
#>  $ face_edge: num [1:527, 1:3] 0.495 0.48 0.408 0.45 0.417 ...
#>  $ vertex   : num [1:279] 0.343 0.542 0.427 0.564 0.439 ...
```
