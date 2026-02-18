# Convert objects to `fm_segm`

Convert objects to `fm_segm`

## Usage

``` r
fm_as_segm(x, ...)

fm_as_segm_list(x, ...)

# S3 method for class 'fm_segm'
fm_as_segm(x, ...)

# S3 method for class 'inla.mesh.segment'
fm_as_segm(x, ...)

# S3 method for class 'sfg'
fm_as_segm(x, ...)

# S3 method for class 'sfc_POINT'
fm_as_segm(x, reverse = FALSE, grp = NULL, is.bnd = TRUE, ...)

# S3 method for class 'sfc_LINESTRING'
fm_as_segm(x, join = TRUE, grp = NULL, reverse = FALSE, ...)

# S3 method for class 'sfc_MULTILINESTRING'
fm_as_segm(x, join = TRUE, grp = NULL, reverse = FALSE, ...)

# S3 method for class 'sfc_POLYGON'
fm_as_segm(x, join = TRUE, grp = NULL, ...)

# S3 method for class 'sfc_MULTIPOLYGON'
fm_as_segm(x, join = TRUE, grp = NULL, ...)

# S3 method for class 'sfc_GEOMETRY'
fm_as_segm(x, grp = NULL, join = TRUE, ...)

# S3 method for class 'sf'
fm_as_segm(x, ...)

# S3 method for class 'matrix'
fm_as_segm(
  x,
  reverse = FALSE,
  grp = NULL,
  is.bnd = FALSE,
  crs = NULL,
  closed = FALSE,
  ...
)

# S3 method for class 'SpatialPoints'
fm_as_segm(x, reverse = FALSE, grp = NULL, is.bnd = TRUE, closed = FALSE, ...)

# S3 method for class 'SpatialPointsDataFrame'
fm_as_segm(x, ...)

# S3 method for class 'Line'
fm_as_segm(x, reverse = FALSE, grp = NULL, crs = NULL, ...)

# S3 method for class 'Lines'
fm_as_segm(x, join = TRUE, grp = NULL, crs = NULL, ...)

# S3 method for class 'SpatialLines'
fm_as_segm(x, join = TRUE, grp = NULL, ...)

# S3 method for class 'SpatialLinesDataFrame'
fm_as_segm(x, ...)

# S3 method for class 'SpatialPolygons'
fm_as_segm(x, join = TRUE, grp = NULL, ...)

# S3 method for class 'SpatialPolygonsDataFrame'
fm_as_segm(x, ...)

# S3 method for class 'Polygons'
fm_as_segm(x, join = TRUE, crs = NULL, grp = NULL, ...)

# S3 method for class 'Polygon'
fm_as_segm(x, crs = NULL, ...)
```

## Arguments

- x:

  Object to be converted.

- ...:

  Arguments passed on to submethods

- reverse:

  logical; When TRUE, reverse the order of the input points. Default
  `FALSE`

- grp:

  if non-null, should be an integer vector of grouping labels for one
  for each segment. Default `NULL`

- is.bnd:

  logical; if `TRUE`, set the boundary flag for the segments. Default
  `TRUE`

- join:

  logical; if `TRUE`, join input segments with common vertices. Default
  `TRUE`

- crs:

  A crs object

- closed:

  logical; whether to treat a point sequence as a closed polygon.
  Default: `FALSE`

## Value

An `fm_segm` or `fm_segm_list` object

## Functions

- `fm_as_segm()`: Convert an object to `fm_segm`.

- `fm_as_segm_list()`: Convert each element, making a `fm_segm_list`
  object

## See also

[`c.fm_segm()`](https://inlabru-org.github.io/fmesher/reference/fm_segm_list.md),
[`c.fm_segm_list()`](https://inlabru-org.github.io/fmesher/reference/fm_segm_list.md),
[`[.fm_segm_list()`](https://inlabru-org.github.io/fmesher/reference/fm_segm_list.md)

Other object creation and conversion:
[`fm_as_collect()`](https://inlabru-org.github.io/fmesher/reference/fm_as_collect.md),
[`fm_as_fm()`](https://inlabru-org.github.io/fmesher/reference/fm_as_fm.md),
[`fm_as_lattice_2d()`](https://inlabru-org.github.io/fmesher/reference/fm_as_lattice_2d.md),
[`fm_as_lattice_Nd()`](https://inlabru-org.github.io/fmesher/reference/fm_as_lattice_Nd.md),
[`fm_as_mesh_1d()`](https://inlabru-org.github.io/fmesher/reference/fm_as_mesh_1d.md),
[`fm_as_mesh_2d()`](https://inlabru-org.github.io/fmesher/reference/fm_as_mesh_2d.md),
[`fm_as_mesh_3d()`](https://inlabru-org.github.io/fmesher/reference/fm_as_mesh_3d.md),
[`fm_as_sfc()`](https://inlabru-org.github.io/fmesher/reference/fm_as_sfc.md),
[`fm_as_tensor()`](https://inlabru-org.github.io/fmesher/reference/fm_as_tensor.md),
[`fm_collect()`](https://inlabru-org.github.io/fmesher/reference/fm_collect.md),
[`fm_lattice_2d()`](https://inlabru-org.github.io/fmesher/reference/fm_lattice_2d.md),
[`fm_lattice_Nd()`](https://inlabru-org.github.io/fmesher/reference/fm_lattice_Nd.md),
[`fm_mesh_1d()`](https://inlabru-org.github.io/fmesher/reference/fm_mesh_1d.md),
[`fm_mesh_2d()`](https://inlabru-org.github.io/fmesher/reference/fm_mesh_2d.md),
[`fm_segm()`](https://inlabru-org.github.io/fmesher/reference/fm_segm.md),
[`fm_simplify()`](https://inlabru-org.github.io/fmesher/reference/fm_simplify.md),
[`fm_tensor()`](https://inlabru-org.github.io/fmesher/reference/fm_tensor.md)

## Examples

``` r
fm_as_segm_list(list(
  fm_segm(fmexample$mesh),
  fm_segm(fmexample$mesh, boundary = FALSE)
))
#> 35 boundary edges (1 group: 1)
#> 55 interior edges (1 group: 1)

(segm <- fm_segm(fmexample$mesh, boundary = FALSE))
#> fm_segm object:
#>   55 interior edges (1 group: 1)
#>   Bounding box = (-3.345586, 2.076840) x (-1.997743, 3.402550) x (0,0)
(segm_sfc <- fm_as_sfc(segm))
#> Geometry set for 1 feature 
#> Geometry type: LINESTRING
#> Dimension:     XYZ
#> Bounding box:  xmin: -3.345586 ymin: -1.997743 xmax: 2.07684 ymax: 3.40255
#> z_range:       zmin: 0 zmax: 0
#> CRS:           NA
#> LINESTRING Z (-3.202599 -0.4432085 0, -3.012083...
(fm_as_segm(segm_sfc))
#> fm_segm object:
#>   55 interior edges (1 group: 1)
#>   Bounding box = (-3.345586, 2.076840) x (-1.997743, 3.402550) x (0,0)
```
