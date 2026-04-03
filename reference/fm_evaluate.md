# Methods for projecting to/from mesh objects

Calculate evaluation information and/or evaluate a function defined on a
mesh or function space.

## Usage

``` r
fm_evaluate(...)

# Default S3 method
fm_evaluate(mesh, field, ...)

# S3 method for class 'fm_evaluator'
fm_evaluate(projector, field, ...)

# S3 method for class 'fm_basis'
fm_evaluate(basis, field, ...)

fm_evaluator(...)

# Default S3 method
fm_evaluator(...)

# S3 method for class 'fm_mesh_3d'
fm_evaluator(mesh, loc = NULL, lattice = NULL, dims = NULL, ...)

# S3 method for class 'fm_mesh_2d'
fm_evaluator(mesh, loc = NULL, lattice = NULL, crs = NULL, ...)

# S3 method for class 'fm_mesh_1d'
fm_evaluator(mesh, loc = NULL, xlim = mesh$interval, dims = 100, ...)

fm_evaluator_lattice(mesh, ...)

# Default S3 method
fm_evaluator_lattice(mesh, dims = 100, ...)

# S3 method for class 'fm_bbox'
fm_evaluator_lattice(mesh, dims = 100, ...)

# S3 method for class 'fm_mesh_2d'
fm_evaluator_lattice(
  mesh,
  xlim = NULL,
  ylim = NULL,
  dims = c(100, 100),
  projection = NULL,
  crs = NULL,
  ...
)
```

## Arguments

- ...:

  Additional arguments passed on to methods.

- mesh:

  An
  [fm_mesh_1d](https://inlabru-org.github.io/fmesher/reference/fm_mesh_1d.md),
  [fm_mesh_2d](https://inlabru-org.github.io/fmesher/reference/fm_mesh_2d.md),
  or other object supported by a sub-method.

- field:

  Basis function weights, one per mesh basis function, describing the
  function to be evaluated at the projection locations

- projector:

  An `fm_evaluator` object.

- basis:

  An
  [fm_basis](https://inlabru-org.github.io/fmesher/reference/fm_basis.md)
  object.

- loc:

  Projection locations. Can be a matrix, `SpatialPoints`,
  `SpatialPointsDataFrame`, `sf`, `sfc`, or `sfg` object.

- lattice:

  An
  [`fm_lattice_2d()`](https://inlabru-org.github.io/fmesher/reference/fm_lattice_2d.md)
  object.

- dims:

  Lattice dimensions.

- crs:

  An optional CRS or inla.CRS object associated with `loc` and/or
  `lattice`.

- xlim:

  X-axis limits for a lattice. For R2 meshes, defaults to covering the
  domain.

- ylim:

  Y-axis limits for a lattice. For R2 meshes, defaults to covering the
  domain.

- projection:

  One of `c("default", "longlat", "longsinlat", "mollweide")`.

## Value

A vector or matrix of the evaluated function

An `fm_evaluator` object

## Methods (by class)

- `fm_evaluate(default)`: The default method calls
  `proj = fm_evaluator(mesh, ...)`, followed by
  `fm_evaluate(proj, field)`.

## Functions

- `fm_evaluate()`: Returns the field function evaluated at the locations
  determined by an `fm_evaluator` object.
  `fm_evaluate(mesh, field = field, ...)` is a shortcut to
  `fm_evaluate(fm_evaluator(mesh, ...), field = field)`.

- `fm_evaluator()`: Returns an `fm_evaluator` list object with
  evaluation information. The `proj` element is a `fm_basis` object,
  containing (at least) a mapping matrix `A` and a logical vector `ok`,
  that indicates which locations were mappable to the input mesh. For
  `fm_mesh_2d` input, `proj` also contains a `bary`
  [fm_bary](https://inlabru-org.github.io/fmesher/reference/fm_bary.md)
  object, with the barycentric coordinates within the triangle each
  input location falls in.

- `fm_evaluator(default)`: The default method calls `fm_basis` and
  creates a basic `fm_evaluator` object

- `fm_evaluator(fm_mesh_3d)`: The `...` arguments are passed on to
  `fm_evaluator_lattice()` if no `loc` or `lattice` is provided.

- `fm_evaluator(fm_mesh_2d)`: The `...` arguments are passed on to
  `fm_evaluator_lattice()` if no `loc` or `lattice` is provided.

- `fm_evaluator_lattice()`: Create a lattice object by default covering
  the input mesh.

- `fm_evaluator_lattice(default)`: Creates an
  [`fm_lattice_2d()`](https://inlabru-org.github.io/fmesher/reference/fm_lattice_2d.md)
  object, by default covering the input mesh.

- `fm_evaluator_lattice(fm_bbox)`: Creates an
  [`fm_lattice_Nd()`](https://inlabru-org.github.io/fmesher/reference/fm_lattice_Nd.md)
  object, by default covering the input mesh.

- `fm_evaluator_lattice(fm_mesh_2d)`: Creates an
  [`fm_lattice_2d()`](https://inlabru-org.github.io/fmesher/reference/fm_lattice_2d.md)
  object, by default covering the input mesh.

## See also

[`fm_mesh_2d()`](https://inlabru-org.github.io/fmesher/reference/fm_mesh_2d.md),
[`fm_mesh_1d()`](https://inlabru-org.github.io/fmesher/reference/fm_mesh_1d.md),
[`fm_lattice_2d()`](https://inlabru-org.github.io/fmesher/reference/fm_lattice_2d.md)

## Author

Finn Lindgren <Finn.Lindgren@gmail.com>

## Examples

``` r
n <- 20
loc <- matrix(runif(n * 2), n, 2)
mesh <- fm_rcdt_2d_inla(loc, refine = list(max.edge = 0.05))
proj <- fm_evaluator(mesh)
field <- cos(mesh$loc[, 1] * 2 * pi * 3) * sin(mesh$loc[, 2] * 2 * pi * 7)
image(proj$x, proj$y, fm_evaluate(proj, field))


# \donttest{
# ## Plotting with inlabru; this feature is not yet in fmesher::geom_fm()
# if (require("inlabru")) {
#  ggplot() +
#    gg(data = mesh, col = field)
# }
# }
```
