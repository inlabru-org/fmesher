# Conversion methods from mesh related objects to sfc

Conversion methods from mesh related objects to sfc

## Usage

``` r
fm_as_sfc(x, ...)

# S3 method for class 'fm_mesh_2d'
fm_as_sfc(x, ..., format = NULL, multi = FALSE)

# S3 method for class 'fm_segm'
fm_as_sfc(x, ..., multi = FALSE)

# S3 method for class 'fm_segm_list'
fm_as_sfc(x, ...)

# S3 method for class 'sfc'
fm_as_sfc(x, ...)

# S3 method for class 'sf'
fm_as_sfc(x, ...)
```

## Arguments

- x:

  An object to be coerced/transformed/converted into another class

- ...:

  Arguments passed on to other methods

- format:

  One of "mesh", "int", "bnd", or "loc". Default "mesh".

- multi:

  logical; if `TRUE`, attempt to a
  `sfc_MULTIPOLYGON/LINESTRING/POINT/GEOMETRYCOLLECTION`, otherwise a
  set of `sfc_POLYGON/LINESTRING/POINT`. Default `FALSE`

## Value

- `fm_as_sfc`: An `sfc_MULTIPOLYGON/LINESTRING/POINT/GEOMETRYCOLLECTION`
  or `sfc_POLYGON/LINESTRING/POINT` object

## Methods (by class)

- `fm_as_sfc(fm_mesh_2d)`: **\[experimental\]**

- `fm_as_sfc(fm_segm)`: **\[experimental\]**

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
fm_as_sfc(fmexample$mesh)
#> Geometry set for 547 features 
#> Geometry type: POLYGON
#> Dimension:     XYZ
#> Bounding box:  xmin: -5.345477 ymin: -3.997839 xmax: 4.08358 ymax: 5.415519
#> z_range:       zmin: 0 zmax: 0
#> CRS:           NA
#> First 5 geometries:
#> POLYGON Z ((0.03596474 1.965644 0, 0.1158377 1....
#> POLYGON Z ((1.504866 -0.622266 0, 1.353357 -0.8...
#> POLYGON Z ((0.7029864 -1.896195 0, 0.6649962 -1...
#> POLYGON Z ((-4.242024 -2.255549 0, -3.13738 -2....
#> POLYGON Z ((-4.321078 -0.5669298 0, -3.847328 -...
fm_as_sfc(fmexample$mesh, multi = TRUE)
#> Geometry set for 1 feature 
#> Geometry type: MULTIPOLYGON
#> Dimension:     XYZ
#> Bounding box:  xmin: -5.345477 ymin: -3.997839 xmax: 4.08358 ymax: 5.415519
#> z_range:       zmin: 0 zmax: 0
#> CRS:           NA
#> MULTIPOLYGON Z (((0.03596474 1.965644 0, 0.1158...
fm_as_sfc(fmexample$mesh, format = "loc")
#> Geometry set for 292 features 
#> Geometry type: POINT
#> Dimension:     XYZ
#> Bounding box:  xmin: -5.345477 ymin: -3.997839 xmax: 4.08358 ymax: 5.415519
#> z_range:       zmin: 0 zmax: 0
#> CRS:           NA
#> First 5 geometries:
#> POINT Z (-5.180541 -0.8956538 0)
#> POINT Z (-4.924679 -1.465076 0)
#> POINT Z (-4.644074 -1.860313 0)
#> POINT Z (-4.242024 -2.255549 0)
#> POINT Z (-3.620028 -2.650785 0)

# Boundary edge conversion to polygons is supported from version 0.4.0.9002:
fm_as_sfc(fmexample$mesh, format = "bnd")
#> Geometry set for 1 feature 
#> Geometry type: POLYGON
#> Dimension:     XYZ
#> Bounding box:  xmin: -5.345477 ymin: -3.997839 xmax: 4.08358 ymax: 5.415519
#> z_range:       zmin: 0 zmax: 0
#> CRS:           NA
#> POLYGON Z ((-5.180541 -0.8956538 0, -4.924679 -...
```
