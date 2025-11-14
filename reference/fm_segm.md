# Make a spatial segment object

Make a spatial segment object

## Usage

``` r
fm_segm(...)

# Default S3 method
fm_segm(loc = NULL, idx = NULL, grp = NULL, is.bnd = TRUE, crs = NULL, ...)

# S3 method for class 'fm_segm'
fm_segm(..., grp = NULL, grp.default = 0L, is.bnd = NULL)

# S3 method for class 'fm_segm_list'
fm_segm(x, grp = NULL, grp.default = 0L, ...)

fm_segm_join(x, grp = NULL, grp.default = 0L, is.bnd = NULL)

fm_segm_split(x, grp = NULL, grp.default = 0L)

# S3 method for class 'inla.mesh.segment'
fm_segm(..., grp.default = 0)

# S3 method for class 'fm_mesh_2d'
fm_segm(x, boundary = TRUE, grp = NULL, ...)

fm_is_bnd(x)

fm_is_bnd(x) <- value
```

## Arguments

- ...:

  Passed on to submethods

- loc:

  Matrix of point locations, or `SpatialPoints`, or `sf`/`sfc` point
  object.

- idx:

  Segment index sequence vector or index pair matrix. The indices refer
  to the rows of `loc`. If `loc==NULL`, the indices will be interpreted
  as indices into the point specification supplied to
  [`fm_rcdt_2d()`](https://inlabru-org.github.io/fmesher/reference/fm_rcdt_2d.md).
  If `is.bnd==TRUE`, defaults to linking all the points in `loc`, as
  `c(1:nrow(loc),1L)`, otherwise `1:nrow(loc)`.

- grp:

  When joining segments, use these group labels for segments instead of
  the original group labels.

- is.bnd:

  `TRUE` if the segments are boundary segments, otherwise `FALSE`.

- crs:

  An optional
  [`fm_crs()`](https://inlabru-org.github.io/fmesher/reference/fm_crs.md),
  [`sf::st_crs()`](https://r-spatial.github.io/sf/reference/st_crs.html)
  or [`sp::CRS()`](https://edzer.github.io/sp/reference/CRS-class.html)
  object

- grp.default:

  If `grp.default` is `NULL`, use these group labels for segments with
  NULL group.

- x:

  Mesh to extract segments from

- boundary:

  logical; if `TRUE`, extract the boundary segments, otherwise interior
  constrain segments.

- value:

  logical

## Value

An `fm_segm` or `fm_segm_list` object

## Methods (by class)

- `fm_segm(fm_segm)`: Join multiple `fm_segm` objects into a single
  `fm_segm` object. If `is.bnd` is non-NULL, it overrides the input
  segment information. Otherwise, it checks if the inputs are
  consistent.

- `fm_segm(fm_segm_list)`: Join `fm_segm` objects from a `fm_segm_list`
  into a single `fm_segm` object. Equivalent to `fm_segm_join(x)`

- `fm_segm(fm_mesh_2d)`: Extract the boundary or interior segments of a
  2d mesh. If `grp` is non-NULL, extracts only segments matching the
  matching the set of groups given by `grp`.

## Functions

- `fm_segm()`: Create a new `fm_segm` object.

- `fm_segm_join()`: Join multiple `fm_segm` objects into a single
  `fm_segm` object. If `is.bnd` is non-NULL, it overrides the segment
  information. Otherwise it checks for consistency.

- `fm_segm_split()`: Split an `fm_segm` object by `grp` into an
  `fm_segm_list` object, optionally keeping only some groups.

## See also

Other object creation and conversion:
[`fm_as_collect()`](https://inlabru-org.github.io/fmesher/reference/fm_as_collect.md),
[`fm_as_fm()`](https://inlabru-org.github.io/fmesher/reference/fm_as_fm.md),
[`fm_as_lattice_2d()`](https://inlabru-org.github.io/fmesher/reference/fm_as_lattice_2d.md),
[`fm_as_lattice_Nd()`](https://inlabru-org.github.io/fmesher/reference/fm_as_lattice_Nd.md),
[`fm_as_mesh_1d()`](https://inlabru-org.github.io/fmesher/reference/fm_as_mesh_1d.md),
[`fm_as_mesh_2d()`](https://inlabru-org.github.io/fmesher/reference/fm_as_mesh_2d.md),
[`fm_as_mesh_3d()`](https://inlabru-org.github.io/fmesher/reference/fm_as_mesh_3d.md),
[`fm_as_segm()`](https://inlabru-org.github.io/fmesher/reference/fm_as_segm.md),
[`fm_as_sfc()`](https://inlabru-org.github.io/fmesher/reference/fm_as_sfc.md),
[`fm_as_tensor()`](https://inlabru-org.github.io/fmesher/reference/fm_as_tensor.md),
[`fm_collect()`](https://inlabru-org.github.io/fmesher/reference/fm_collect.md),
[`fm_lattice_2d()`](https://inlabru-org.github.io/fmesher/reference/fm_lattice_2d.md),
[`fm_lattice_Nd()`](https://inlabru-org.github.io/fmesher/reference/fm_lattice_Nd.md),
[`fm_mesh_1d()`](https://inlabru-org.github.io/fmesher/reference/fm_mesh_1d.md),
[`fm_mesh_2d()`](https://inlabru-org.github.io/fmesher/reference/fm_mesh_2d.md),
[`fm_simplify()`](https://inlabru-org.github.io/fmesher/reference/fm_simplify.md),
[`fm_tensor()`](https://inlabru-org.github.io/fmesher/reference/fm_tensor.md)

## Examples

``` r
fm_segm(rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1)), is.bnd = FALSE)
#> fm_segm object:
#>   3 interior edges
#>   Bounding box = (0,1) x (0,1) x (0,0)
fm_segm(rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1)), is.bnd = TRUE)
#> fm_segm object:
#>   4 boundary edges
#>   Bounding box = (0,1) x (0,1) x (0,0)

fm_segm_join(fmexample$boundary_fm)
#> fm_segm object:
#>   72 boundary edges (1 group: 1)
#>   Bounding box = (-5.331027, 4.061656) x (-3.998161, 5.415609) x (0,0)

fm_segm(fmexample$mesh, boundary = TRUE)
#> fm_segm object:
#>   29 boundary edges (1 group: 1)
#>   Bounding box = (-5.331027, 4.061656) x (-3.998161, 5.415609) x (0,0)
fm_segm(fmexample$mesh, boundary = FALSE)
#> fm_segm object:
#>   54 interior edges (1 group: 1)
#>   Bounding box = (-3.344418, 2.076846) x (-1.995602, 3.404937) x (0,0)
```
