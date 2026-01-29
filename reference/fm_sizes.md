# fm_sizes

**\[experimental\]** Compute effective sizes of faces/cells and vertices
in a mesh

## Usage

``` r
fm_sizes(...)

# S3 method for class 'fm_mesh_2d'
fm_sizes(mesh, ..., method = "R")

# S3 method for class 'fm_mesh_3d'
fm_sizes(mesh, ...)
```

## Arguments

- ...:

  Passed on to submethods

- mesh:

  object of a supported mesh class

- method:

  character; "R" or "Rcpp". For "S2" manifolds, the "Rcpp" method is
  always used. The "R" method is currently faster, due to the cost of
  building internal data structures in the C++ code.

## Value

A `list` with elements of simplex size information. For 2D meshes:

- `face`:

  Vector with the area of each triangle

- `vertex`:

  Vector with the triangle area apportioned to each vertex

- `face_edge`:

  A matrix with one row per triangle and 3 columns, with edge lengths
  for the edge opposing each triangle vertex.

For 3D meshes:

- `cell`:

  Vector with the volume of each tetrahedron

- `vertex`:

  Vector with the tetrahedron volume apportioned to each vertex

- `cell_face`:

  A matrix with one row per cell and 4 columns, with triangle areas for
  the triangle opposing each tetrahedron vertex.

- `cell_edge`:

  A matrix with one row per cell and 4 columns, with edge lengths for
  the edge anchored at each vertex, pointing to the next vertex in the
  internal ordering.

## Examples

``` r
str(fm_sizes(fmexample$mesh))
#> List of 3
#>  $ face     : num [1:527] 0.0641 0.0799 0.0794 0.0734 0.0383 ...
#>  $ face_edge: num [1:527, 1:3] 0.495 0.48 0.408 0.45 0.417 ...
#>  $ vertex   : num [1:279] 0.343 0.542 0.427 0.564 0.439 ...
```
