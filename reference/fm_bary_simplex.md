# Extract Simplex information for Barycentric coordinates

Extract the simplex vertex information for a combination of a mesh and
[fm_bary](https://inlabru-org.github.io/fmesher/reference/fm_bary.md)
coordinates.

## Usage

``` r
fm_bary_simplex(mesh, bary = NULL, ...)

# S3 method for class 'fm_mesh_2d'
fm_bary_simplex(mesh, bary = NULL, ...)

# S3 method for class 'fm_mesh_3d'
fm_bary_simplex(mesh, bary = NULL, ...)

# S3 method for class 'fm_mesh_1d'
fm_bary_simplex(mesh, bary = NULL, ...)

# S3 method for class 'fm_lattice_2d'
fm_bary_simplex(mesh, bary = NULL, ...)

# S3 method for class 'fm_lattice_Nd'
fm_bary_simplex(mesh, bary = NULL, ...)
```

## Arguments

- mesh:

  A mesh object, e.g.
  [fm_mesh_2d](https://inlabru-org.github.io/fmesher/reference/fm_mesh_2d.md)
  or
  [fm_mesh_1d](https://inlabru-org.github.io/fmesher/reference/fm_mesh_1d.md).

- bary:

  An
  [fm_bary](https://inlabru-org.github.io/fmesher/reference/fm_bary.md)
  object. If NULL, return the full simplex information for the mesh.

- ...:

  Further arguments potentially used by sub-methods.

## Value

A matrix of vertex indices, one row per point in `bary`.

## Methods (by class)

- `fm_bary_simplex(fm_mesh_2d)`: Extract the triangle vertex indices for
  a 2D mesh

- `fm_bary_simplex(fm_mesh_3d)`: Extract the tetrahedron vertex indices
  for a 3D mesh

- `fm_bary_simplex(fm_mesh_1d)`: Extract the edge vertex indices for a
  1D mesh

- `fm_bary_simplex(fm_lattice_2d)`: Extract the cell vertex indices for
  a 2D lattice

- `fm_bary_simplex(fm_lattice_Nd)`: Extract the cell vertex indices for
  a ND lattice

## See also

[`fm_bary()`](https://inlabru-org.github.io/fmesher/reference/fm_bary.md),
[`fm_bary_loc()`](https://inlabru-org.github.io/fmesher/reference/fm_bary_loc.md)

## Examples

``` r
bary <- fm_bary(fmexample$mesh, fmexample$loc_sf)
fm_bary_simplex(fmexample$mesh, bary)
#>       [,1] [,2] [,3]
#>  [1,]  172  149  103
#>  [2,]  130   84  144
#>  [3,]  197  114  200
#>  [4,]   88  125  162
#>  [5,]  199  213  112
#>  [6,]  170   97  132
#>  [7,]  174  111  175
#>  [8,]  175  215   89
#>  [9,]  175  215   89
#> [10,]  214   85  168
(m <- fm_mesh_3d(
  matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0), 4, 3, byrow = TRUE),
  matrix(c(1, 2, 3, 4), 1, 4, byrow = TRUE)
))
#> fm_mesh_3d object:
#>   Manifold:  R3
#>   V / E / T / Tet:   4 / 6 / 4 / 1
#>   Euler char.:   1
#>   Bounding box: (0,1) x (0,1) x (0,1)
#>   Basis d.o.f.:  4
(bary <- fm_bary(m, rbind(
  cbind(0.1, 0.2, 0.3),
  cbind(-0.1, 0.2, 0.3)
)))
#> # A tibble: 2 × 2
#>   index where[,1]  [,2]  [,3]  [,4]
#> * <int>     <dbl> <dbl> <dbl> <dbl>
#> 1     1       0.1   0.2   0.3   0.4
#> 2    NA      NA    NA    NA    NA  
fm_bary_simplex(m, bary)
#>      [,1] [,2] [,3] [,4]
#> [1,]    1    2    3    4
#> [2,]   NA   NA   NA   NA
mesh1 <- fm_mesh_1d(1:4)
(bary1 <- fm_bary(mesh1, seq(0, 5, by = 0.5)))
#> # A tibble: 11 × 2
#>    index where[,1]  [,2]
#>  * <int>     <dbl> <dbl>
#>  1     1       2    -1  
#>  2     1       1.5  -0.5
#>  3     1       1     0  
#>  4     1       0.5   0.5
#>  5     2       1     0  
#>  6     2       0.5   0.5
#>  7     3       1     0  
#>  8     3       0.5   0.5
#>  9     3       0     1  
#> 10     3      -0.5   1.5
#> 11     3      -1     2  
(bary1 <- fm_bary(mesh1, seq(0, 5, by = 0.5), restricted = TRUE))
#> # A tibble: 11 × 2
#>    index where[,1]  [,2]
#>  * <int>     <dbl> <dbl>
#>  1     1       1     0  
#>  2     1       1     0  
#>  3     1       1     0  
#>  4     1       0.5   0.5
#>  5     2       1     0  
#>  6     2       0.5   0.5
#>  7     3       1     0  
#>  8     3       0.5   0.5
#>  9     3       0     1  
#> 10     3       0     1  
#> 11     3       0     1  
fm_bary_simplex(mesh1, bary1)
#>       [,1] [,2]
#>  [1,]    1    2
#>  [2,]    1    2
#>  [3,]    1    2
#>  [4,]    1    2
#>  [5,]    2    3
#>  [6,]    2    3
#>  [7,]    3    4
#>  [8,]    3    4
#>  [9,]    3    4
#> [10,]    3    4
#> [11,]    3    4
m <- fm_lattice_2d(x = 1:3, y = 1:4)
bary <- fm_bary(m, cbind(1.5, 3.2))
fm_bary_simplex(m, bary)
#>      [,1] [,2] [,3] [,4]
#> [1,]    7    8   11   10
m <- fm_lattice_Nd(list(x = 1:3, y = 1:4, z = 1:2))
(bary <- fm_bary(m, cbind(1.5, 3.2, 1.5)))
#> # A tibble: 1 × 2
#>   index where[,1]  [,2]   [,3]   [,4]  [,5]  [,6]   [,7]   [,8]
#> * <int>     <dbl> <dbl>  <dbl>  <dbl> <dbl> <dbl>  <dbl>  <dbl>
#> 1     5       0.2   0.2 0.0500 0.0500   0.2   0.2 0.0500 0.0500
(fm_bary_simplex(m, bary))
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8]
#> [1,]    7    8   10   11   19   20   22   23
fm_bary_loc(m, bary)
#>      [,1] [,2] [,3]
#> [1,]  1.5  3.2  1.5
```
