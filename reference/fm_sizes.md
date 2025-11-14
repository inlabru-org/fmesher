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
the faces/cells and vertices, respectively.

## Examples

``` r
str(fm_sizes(fmexample$mesh))
#> List of 2
#>  $ face  : num [1:527] 0.0214 0.0266 0.0265 0.0245 0.0128 ...
#>  $ vertex: num [1:279] 0.114 0.181 0.142 0.188 0.146 ...
```
