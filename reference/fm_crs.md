# Obtain coordinate reference system object

Obtain an
[`sf::crs`](https://r-spatial.github.io/sf/reference/coerce-methods.html)
or `fm_crs` object from a spatial object, or convert crs information to
construct a new
[`sf::crs`](https://r-spatial.github.io/sf/reference/coerce-methods.html)
object.

## Usage

``` r
fm_crs(x, ..., units = NULL, oblique = NULL)

fm_crs_oblique(x)

# S3 method for class 'fm_crs'
st_crs(x, ...)

# S3 method for class 'fm_crs'
x$name

# Default S3 method
fm_crs(x, ..., units = NULL, oblique = NULL)

# S3 method for class 'crs'
fm_crs(x, ..., units = NULL, oblique = NULL)

# S3 method for class 'fm_crs'
fm_crs(x, ..., units = NULL, oblique = NULL)

# S3 method for class 'fm_CRS'
fm_crs(x, ..., units = NULL, oblique = NULL)

# S3 method for class 'character'
fm_crs(x, ..., units = NULL, oblique = NULL)

# S3 method for class 'Spatial'
fm_crs(x, ..., units = NULL, oblique = NULL)

# S3 method for class 'SpatVector'
fm_crs(x, ..., units = NULL, oblique = NULL)

# S3 method for class 'SpatRaster'
fm_crs(x, ..., units = NULL, oblique = NULL)

# S3 method for class 'sf'
fm_crs(x, ..., units = NULL, oblique = NULL)

# S3 method for class 'sfc'
fm_crs(x, ..., units = NULL, oblique = NULL)

# S3 method for class 'sfg'
fm_crs(x, ..., units = NULL, oblique = NULL)

# S3 method for class 'fm_mesh_2d'
fm_crs(x, ..., units = NULL, oblique = NULL)

# S3 method for class 'fm_mesh_1d'
fm_crs(x, ..., units = NULL, oblique = NULL)

# S3 method for class 'fm_mesh_3d'
fm_crs(x, ..., units = NULL, oblique = NULL)

# S3 method for class 'fm_tensor'
fm_crs(x, ..., units = NULL, oblique = NULL, .multi = FALSE)

# S3 method for class 'fm_collect'
fm_crs(x, ..., units = NULL, oblique = NULL, .multi = FALSE)

# S3 method for class 'fm_lattice_2d'
fm_crs(x, ..., units = NULL, oblique = NULL)

# S3 method for class 'fm_segm'
fm_crs(x, ..., units = NULL, oblique = NULL)

# S3 method for class 'fm_list'
fm_crs(x, ..., units = NULL, oblique = NULL)

# S3 method for class 'matrix'
fm_crs(x, ..., units = NULL, oblique = NULL)

# S3 method for class 'fm_list'
fm_CRS(x, ..., units = NULL, oblique = NULL)

fm_wkt_predef()

# S3 method for class 'inla.CRS'
fm_crs(x, ..., units = NULL, oblique = NULL)
```

## Arguments

- x:

  Object to convert to `crs` or to extract `crs` information from. If
  `character`, a string suitable for `sf::st_crs(x)`, or the name of a
  predefined `wkt` string from “names(fm_wkt_predef())\`.

- ...:

  Additional parameters. Not currently in use.

- units:

  character; if non-NULL, `fm_length_unit()<-` is called to change the
  length units of the crs object. If `NULL` (default), the length units
  are not changed. (From version `0.3.0.9013`)

- oblique:

  Numeric vector of length at most 4 of rotation angles (in degrees) for
  an oblique projection, all values defaulting to zero. The values
  indicate (longitude, latitude, orientation, orbit), as explained in
  the Details section below. When `oblique` is non-NULL, used to
  override the obliqueness parameters of a `fm_crs` object. When `NA`,
  remove obliqueness from the object, resulting in a return class of
  [`sf::st_crs()`](https://r-spatial.github.io/sf/reference/st_crs.html).
  When `NULL`, pass though any oblique information in the object,
  returning an `fm_crs()` object if needed.

- name:

  element name

- .multi:

  logical; If `TRUE`, return a list of `fm_crs` objects for classes that
  support multiple spaces. Default `FALSE`

## Value

Either an
[`sf::crs`](https://r-spatial.github.io/sf/reference/coerce-methods.html)
object or an `fm_crs` object, depending on if the coordinate reference
system described by the parameters can be expressed with a pure `crs`
object or not.

A `crs` object
([`sf::st_crs()`](https://r-spatial.github.io/sf/reference/st_crs.html))
or a `fm_crs` object. An S3 `fm_crs` object is a list with elements
`crs` and `oblique`.

`fm_wkt_predef` returns a WKT2 string defining a projection

## Details

The first two elements of the `oblique` vector are the (longitude,
latitude) coordinates for the oblique centre point. The third value
(orientation) is a counter-clockwise rotation angle for an observer
looking at the centre point from outside the sphere. The fourth value is
the quasi-longitude (orbit angle) for a rotation along the oblique
observers equator.

Simple oblique: `oblique=c(0, 45)`

Polar: `oblique=c(0, 90)`

Quasi-transversal: `oblique=c(0, 0, 90)`

Satellite orbit viewpoint:
`oblique=c(lon0-time*v1, 0, orbitangle, orbit0+time*v2)`, where `lon0`
is the longitude at which a satellite orbit crosses the equator at
`time=0`, when the satellite is at an angle `orbit0` further along in
its orbit. The orbital angle relative to the equatorial plane is
`orbitangle`, and `v1` and `v2` are the angular velocities of the planet
and the satellite, respectively. Note that "forward" from the
satellite's point of view is "to the right" in the projection.

When `oblique[2]` or `oblique[3]` are non-zero, the resulting projection
is only correct for perfect spheres.

## Methods (by class)

- `fm_crs(fm_tensor)`: By default returns the crs of the first space in
  the tensor product space.

- `fm_crs(fm_collect)`: By default returns the crs of the first space in
  the collection.

- `fm_crs(fm_list)`: returns a list of 'crs' objects, one for each list
  element

## Methods (by generic)

- `st_crs(fm_crs)`: `st_crs(x, ...)` is equivalent to
  `fm_crs(x, oblique = NA, ...)` when `x` is a `fm_crs` object.

- `$`: For a `fm_crs` object `x`, `x$name` calls the accessor method for
  the `crs` object inside it. If `name` is "crs", the internal crs
  object itself is returned. If `name` is "oblique", the internal
  oblique angle parameter vector is returned.

## Functions

- `fm_crs_oblique()`: Return `NA` for object with no oblique
  information, and otherwise a length 4 numeric vector.

- `fm_CRS(fm_list)`: returns a list of 'CRS' objects, one for each list
  element

## See also

[`sf::st_crs()`](https://r-spatial.github.io/sf/reference/st_crs.html),
[`fm_crs_wkt`](https://inlabru-org.github.io/fmesher/reference/fm_crs_wkt.md)

fm_crs_is_null

[`fm_crs<-()`](https://inlabru-org.github.io/fmesher/reference/fm_crs-set.md),
[`fm_crs_oblique<-()`](https://inlabru-org.github.io/fmesher/reference/fm_crs-set.md)

## Author

Finn Lindgren <Finn.Lindgren@gmail.com>

## Examples

``` r
crs1 <- fm_crs("longlat_globe")
crs2 <- fm_crs("lambert_globe")
crs3 <- fm_crs("mollweide_norm")
crs4 <- fm_crs("hammer_globe")
crs5 <- fm_crs("sphere")
crs6 <- fm_crs("globe")
names(fm_wkt_predef())
#>  [1] "hammer_norm"     "lambert_norm"    "longlat_norm"    "mollweide_norm" 
#>  [5] "hammer_globe"    "lambert_globe"   "longlat_globe"   "mollweide_globe"
#>  [9] "sphere"          "globe"          
```
