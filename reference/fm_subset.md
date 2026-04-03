# Extract a subset of a mesh

**\[experimental\]** (from version `0.5.0.9003`) Constructs a new mesh
based on a subset of the triangles of an existing mesh. The current
version drops any edge constraint information from the mesh.

## Usage

``` r
fm_subset(mesh, t_sub)
```

## Arguments

- mesh:

  an mesh to subset

- t_sub:

  triangle or tetrahedron indices.

## Value

A subset mesh.

## Author

Finn Lindgren <Finn.Lindgren@gmail.com>

## Examples

``` r
mesh_sub <- fm_subset(fmexample$mesh, 1:100)
mesh_sub
#> fm_mesh_2d object:
#>   Manifold:  R2
#>   V / E / T: 173 / 267 / 100
#>   Euler char.:   6
#>   Constraints:   Boundary: 234 boundary edges (1 group: 0), Interior: 0 edges
#>   Bounding box: (-4.644074, 4.004812) x (-3.997839, 3.275186)
#>   Basis d.o.f.:  173
plot(mesh_sub)


if (requireNamespace("geometry", quietly = TRUE)) {
  m <- fm_delaunay_3d(matrix(rnorm(30), 10, 3))
  print(m)
  print(fm_subset(m, seq_len(min(5, nrow(m$graph$tv)))))
}
#> fm_mesh_3d object:
#>   Manifold:  R3
#>   V / E / T / Tet:   10 / 36 / 50 / 23
#>   Euler char.:   1
#>   Bounding box: (-1.307042, 2.325888) x (-1.064937, 2.113168) x (-1.2234204, 0.9558412)
#>   Basis d.o.f.:  10
#> fm_mesh_3d object:
#>   Manifold:  R3
#>   V / E / T / Tet:   8 / 19 / 17 / 5
#>   Euler char.:   1
#>   Bounding box: (-1.307042, 2.325888) x (-1.064937, 2.113168) x (-1.2234204, 0.9558412)
#>   Basis d.o.f.:  8
```
