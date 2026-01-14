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
#> Geometry set for 527 features 
#> Geometry type: POLYGON
#> Dimension:     XYZ
#> Bounding box:  xmin: -5.331027 ymin: -3.998161 xmax: 4.061656 ymax: 5.415609
#> z_range:       zmin: 0 zmax: 0
#> CRS:           NA
#> First 5 geometries:
#> POLYGON Z ((1.603919 0.1315077 0, 1.477037 0.36...
#> POLYGON Z ((1.929231 -0.2412966 0, 1.426833 -0....
#> POLYGON Z ((0.7169789 -1.03035 0, 0.3472101 -0....
#> POLYGON Z ((-2.348761 -0.5948482 0, -2.384393 -...
#> POLYGON Z ((-3.269508 0.4471134 0, -3.012816 0....
fm_as_sfc(fmexample$mesh, multi = TRUE)
#> Geometry set for 1 feature 
#> Geometry type: MULTIPOLYGON
#> Dimension:     XYZ
#> Bounding box:  xmin: -5.331027 ymin: -3.998161 xmax: 4.061656 ymax: 5.415609
#> z_range:       zmin: 0 zmax: 0
#> CRS:           NA
#> MULTIPOLYGON Z (((1.603919 0.1315077 0, 1.47703...
fm_as_sfc(fmexample$mesh, format = "loc")
#> Geometry set for 279 features 
#> Geometry type: POINT
#> Dimension:     XYZ
#> Bounding box:  xmin: -5.331027 ymin: -3.998161 xmax: 4.061656 ymax: 5.415609
#> z_range:       zmin: 0 zmax: 0
#> CRS:           NA
#> First 5 geometries:
#> POINT Z (-4.991262 1.478543 0)
#> POINT Z (-5.287832 0.6496857 0)
#> POINT Z (-5.331027 -0.2295705 0)
#> POINT Z (-5.08812 -1.151174 0)
#> POINT Z (-4.710559 -1.781097 0)

# Boundary edge conversion to polygons is supported from version 0.4.0.9002:
fm_as_sfc(fmexample$mesh, format = "bnd")
#> Geometry set for 1 feature 
#> Geometry type: POLYGON
#> Dimension:     XYZ
#> Bounding box:  xmin: -5.331027 ymin: -3.998161 xmax: 4.061656 ymax: 5.415609
#> z_range:       zmin: 0 zmax: 0
#> CRS:           NA
#> POLYGON Z ((-4.991262 1.478543 0, -5.287832 0.6...
```
