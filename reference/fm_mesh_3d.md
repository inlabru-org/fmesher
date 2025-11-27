# Construct a 3D tetrahedralisation

Constructs a 3D tetrahedralisation object.

## Usage

``` r
fm_mesh_3d(loc = NULL, tv = NULL, ...)

fm_delaunay_3d(loc, ...)
```

## Arguments

- loc:

  Input coordinates that should be part of the mesh. Can be a matrix,
  `sf`, `sfc`, `SpatialPoints`, or other object supported by
  [`fm_unify_coords()`](https://inlabru-org.github.io/fmesher/reference/fm_unify_coords.md).

- tv:

  Tetrahedron indices, as a N-by-4 index vector into `loc`

- ...:

  Currently unused.

## Value

An `fm_mesh_3d` object

## Functions

- `fm_delaunay_3d()`: Construct a plain Delaunay triangulation in 3D.
  Requires the `geometry` package.

## Examples

``` r
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
(m <- fm_delaunay_3d(matrix(rnorm(30), 10, 3)))
#> fm_mesh_3d object:
#>   Manifold:  R3
#>   V / E / T / Tet:   10 / 34 / 44 / 19
#>   Euler char.:   1
#>   Bounding box: (-1.431271, 1.382911) x (-1.433321, 1.460110) x (-2.177576, 1.893360)
#>   Basis d.o.f.:  10
```
