# Extract Euclidean Sgeometry from Barycentric coordinates

Extract the Euclidean coordinates for location identified by an
[fm_bary](https://inlabru-org.github.io/fmesher/reference/fm_bary.md)
object. This acts as the inverse of
[`fm_bary()`](https://inlabru-org.github.io/fmesher/reference/fm_bary.md).

## Usage

``` r
fm_bary_loc(mesh, bary = NULL, ..., format = NULL)

# S3 method for class 'fm_mesh_2d'
fm_bary_loc(mesh, bary = NULL, ..., format = NULL)

# S3 method for class 'fm_mesh_3d'
fm_bary_loc(mesh, bary = NULL, ..., format = NULL)

# S3 method for class 'fm_mesh_1d'
fm_bary_loc(mesh, bary = NULL, ..., format = NULL)

# S3 method for class 'fm_lattice_2d'
fm_bary_loc(mesh, bary = NULL, ..., format = NULL)

# S3 method for class 'fm_lattice_Nd'
fm_bary_loc(mesh, bary = NULL, ..., format = NULL)
```

## Arguments

- mesh:

  A mesh object, e.g.
  [fm_mesh_2d](https://inlabru-org.github.io/fmesher/reference/fm_mesh_2d.md)
  or
  [fm_mesh_1d](https://inlabru-org.github.io/fmesher/reference/fm_mesh_1d.md).

- bary:

  An `fm_bary` object. If `NULL`, return the mesh nodes is the mesh
  class supports it, otherwise gives an error.

- ...:

  Further arguments potentially used by sub-methods.

- format:

  Optional format for the output. If `NULL`, the output format is
  determined by the default for the mesh object.

## Value

Output format depends on the mesh `class`.

## Methods (by class)

- `fm_bary_loc(fm_mesh_2d)`: Extract points on a triangle mesh.
  Implemented formats are `"matrix"` (default) and `"sf"`.

- `fm_bary_loc(fm_mesh_3d)`: Extract points on a tetrahedron mesh.
  Implemented format is `"matrix"` (default).

- `fm_bary_loc(fm_mesh_1d)`: Extract points on a 1D mesh. Implemented
  formats are `"numeric"` (default).

- `fm_bary_loc(fm_lattice_2d)`: Extract points on a 2D lattice.
  Implemented formats are `"matrix"` (default) and `"sf"`.

- `fm_bary_loc(fm_lattice_Nd)`: Extract points on a ND lattice.

## See also

[`fm_bary()`](https://inlabru-org.github.io/fmesher/reference/fm_bary.md),
[`fm_bary_simplex()`](https://inlabru-org.github.io/fmesher/reference/fm_bary_simplex.md)

## Examples

``` r
head(fm_bary_loc(fmexample$mesh))
#>           [,1]       [,2] [,3]
#> [1,] -4.991262  1.4785425    0
#> [2,] -5.287832  0.6496857    0
#> [3,] -5.331027 -0.2295705    0
#> [4,] -5.088120 -1.1511736    0
#> [5,] -4.710559 -1.7810970    0
#> [6,] -4.191253 -2.3004024    0
bary <- fm_bary(fmexample$mesh, fmexample$loc_sf)
fm_bary_loc(fmexample$mesh, bary, format = "matrix")
#>             [,1]        [,2] [,3]
#>  [1,] -1.2070657 -0.47719270    0
#>  [2,]  0.2774292 -0.99838644    0
#>  [3,]  1.0844412 -0.77625389    0
#>  [4,] -2.3456977  0.06445882    0
#>  [5,]  0.4291247  0.95949406    0
#>  [6,]  0.5060559 -0.11028549    0
#>  [7,] -0.5747400 -0.51100951    0
#>  [8,] -0.5466319 -0.91119542    0
#>  [9,] -0.5644520 -0.83717168    0
#> [10,] -0.8900378  2.41583518    0
fm_bary_loc(fmexample$mesh, bary, format = "sf")
#> Simple feature collection with 10 features and 0 fields
#> Geometry type: POINT
#> Dimension:     XYZ
#> Bounding box:  xmin: -2.345698 ymin: -0.9983864 xmax: 1.084441 ymax: 2.415835
#> CRS:           NA
#>                          geometry
#> 1  POINT Z (-1.207066 -0.47719...
#> 2  POINT Z (0.2774292 -0.99838...
#> 3  POINT Z (1.084441 -0.776253...
#> 4  POINT Z (-2.345698 0.064458...
#> 5  POINT Z (0.4291247 0.959494...
#> 6  POINT Z (0.5060559 -0.11028...
#> 7  POINT Z (-0.57474 -0.511009...
#> 8  POINT Z (-0.5466319 -0.9111...
#> 9  POINT Z (-0.564452 -0.83717...
#> 10 POINT Z (-0.8900378 2.41583...
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
fm_bary_loc(m, bary)
#>      [,1] [,2] [,3]
#> [1,]  0.1  0.2  0.3
#> [2,]   NA   NA   NA
mesh1 <- fm_mesh_1d(1:4)
fm_bary_loc(mesh1)
#> [1] 1 2 3 4
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
fm_bary_loc(mesh1, bary1)
#>  [1] 0.0 0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0
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
fm_bary_loc(mesh1, bary1)
#>  [1] 1.0 1.0 1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.0 4.0
fm_basis(mesh1, bary1)
#> 11 x 4 sparse Matrix of class "dgCMatrix"
#>                      
#>  [1,] 1.0 0.0 .   .  
#>  [2,] 1.0 0.0 .   .  
#>  [3,] 1.0 0.0 .   .  
#>  [4,] 0.5 0.5 .   .  
#>  [5,] .   1.0 0.0 .  
#>  [6,] .   0.5 0.5 .  
#>  [7,] .   .   1.0 0.0
#>  [8,] .   .   0.5 0.5
#>  [9,] .   .   0.0 1.0
#> [10,] .   .   0.0 1.0
#> [11,] .   .   0.0 1.0
(bary1 <- fm_bary(mesh1, bary1, method = "nearest"))
#> # A tibble: 11 × 2
#>    index where[,1]
#>  * <int>     <dbl>
#>  1     1         1
#>  2     1         1
#>  3     1         1
#>  4     1         1
#>  5     2         1
#>  6     2         1
#>  7     3         1
#>  8     3         1
#>  9     4         1
#> 10     4         1
#> 11     4         1
fm_bary_loc(mesh1, bary1)
#>  [1] 1 1 1 1 2 2 3 3 4 4 4
fm_basis(mesh1, bary1)
#> 11 x 4 sparse Matrix of class "dgCMatrix"
#>              
#>  [1,] 1 0 . .
#>  [2,] 1 0 . .
#>  [3,] 1 0 . .
#>  [4,] 1 0 . .
#>  [5,] . 1 0 .
#>  [6,] . 1 0 .
#>  [7,] . . 1 0
#>  [8,] . . 1 0
#>  [9,] . . 0 1
#> [10,] . . 0 1
#> [11,] . . 0 1
(bary1 <- fm_bary(mesh1, bary1, method = "linear"))
#> # A tibble: 11 × 2
#>    index where[,1]  [,2]
#>  * <int>     <dbl> <dbl>
#>  1     1         1     0
#>  2     1         1     0
#>  3     1         1     0
#>  4     1         1     0
#>  5     2         1     0
#>  6     2         1     0
#>  7     3         1     0
#>  8     3         1     0
#>  9     3         0     1
#> 10     3         0     1
#> 11     3         0     1
fm_bary_loc(mesh1, bary1)
#>  [1] 1 1 1 1 2 2 3 3 4 4 4
fm_basis(mesh1, bary1)
#> 11 x 4 sparse Matrix of class "dgCMatrix"
#>              
#>  [1,] 1 0 . .
#>  [2,] 1 0 . .
#>  [3,] 1 0 . .
#>  [4,] 1 0 . .
#>  [5,] . 1 0 .
#>  [6,] . 1 0 .
#>  [7,] . . 1 0
#>  [8,] . . 1 0
#>  [9,] . . 0 1
#> [10,] . . 0 1
#> [11,] . . 0 1
m <- fm_lattice_2d(x = 1:3, y = 1:4)
head(fm_bary_loc(m))
#>      [,1] [,2]
#> [1,]    1    1
#> [2,]    2    1
#> [3,]    3    1
#> [4,]    1    2
#> [5,]    2    2
#> [6,]    3    2
(bary <- fm_bary(m, cbind(1.5, 3.2)))
#> # A tibble: 1 × 2
#>   index where[,1]  [,2]  [,3]  [,4]
#> * <int>     <dbl> <dbl> <dbl> <dbl>
#> 1     5       0.4   0.4 0.100 0.100
fm_bary_loc(m, bary, format = "matrix")
#>      [,1] [,2]
#> [1,]  1.5  3.2
fm_bary_loc(m, bary, format = "sf")
#> Simple feature collection with 1 feature and 0 fields
#> Geometry type: POINT
#> Dimension:     XY
#> Bounding box:  xmin: 1.5 ymin: 3.2 xmax: 1.5 ymax: 3.2
#> CRS:           NA
#>          geometry
#> 1 POINT (1.5 3.2)
m <- fm_lattice_Nd(list(x = 1:3, y = 1:4, z = 1:2))
head(fm_bary_loc(m))
#>      x y z
#> [1,] 1 1 1
#> [2,] 2 1 1
#> [3,] 3 1 1
#> [4,] 1 2 1
#> [5,] 2 2 1
#> [6,] 3 2 1
(bary <- fm_bary(m, cbind(1.5, 3.2, 1.5)))
#> # A tibble: 1 × 2
#>   index where[,1]  [,2]   [,3]   [,4]  [,5]  [,6]   [,7]   [,8]
#> * <int>     <dbl> <dbl>  <dbl>  <dbl> <dbl> <dbl>  <dbl>  <dbl>
#> 1     5       0.2   0.2 0.0500 0.0500   0.2   0.2 0.0500 0.0500
fm_bary_loc(m, bary)
#>      [,1] [,2] [,3]
#> [1,]  1.5  3.2  1.5
```
