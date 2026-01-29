# Subdivide triangles

Internal C++ method.

Subdivide a mesh with congruent and anti-congruent subtriangles

## Usage

``` r
fmesher_subdivide(
  mesh_loc,
  mesh_tv,
  mesh_boundary,
  mesh_interior,
  subdivisions,
  options
)
```

## Arguments

- mesh_loc:

  numeric matrix; mesh vertex coordinates

- mesh_tv:

  3-column integer matrix with 0-based vertex indices for each triangle

- mesh_boundary:

  2-column integer matrix with 0-based vertex indices for boundary
  constraints, currently ignored

- mesh_interior:

  2-column integer matrix with 0-based vertex indices for interior
  constraints, currently ignored

- subdivisions:

  integer; number of new points along each edge.

- options:

  list of triangulation options (`sphere_tolerance`)

## Value

A list of new `loc` and `tv` information, and `bary_index` and
`bary_where` containing the fm_bary information for the new points. Can
be e.g. used to construct an interpolation mapping matrix from the old
to new mesh.

## See also

[`fm_subdivide()`](https://inlabru-org.github.io/fmesher/reference/fm_subdivide.md)

## Examples

``` r
mesh <- fm_mesh_2d(
  boundary = fm_segm(rbind(c(0,0), c(1,0), c(1,1), c(0, 1)), is.bnd = TRUE)
)
new_mesh <- fm_subdivide(mesh, n = 3)
plot(new_mesh, edge.color = 2)
plot(mesh, add = TRUE, edge.color = 1)
```
