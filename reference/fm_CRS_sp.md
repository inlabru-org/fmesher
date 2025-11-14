# Create a coordinate reference system object

Creates either a CRS object or an inla.CRS object, describing a
coordinate reference system

## Usage

``` r
fm_CRS(x, ..., units = NULL, oblique = NULL)

# S3 method for class 'fm_CRS'
is.na(x)

# S3 method for class 'crs'
fm_CRS(x, ..., units = NULL, oblique = NULL)

# S3 method for class 'fm_crs'
fm_CRS(x, ..., units = NULL, oblique = NULL)

# S3 method for class 'Spatial'
fm_CRS(x, ..., units = NULL, oblique = NULL)

# S3 method for class 'fm_CRS'
fm_CRS(x, ..., units = NULL, oblique = NULL)

# S3 method for class 'SpatVector'
fm_CRS(x, ..., units = NULL, oblique = NULL)

# S3 method for class 'SpatRaster'
fm_CRS(x, ..., units = NULL, oblique = NULL)

# S3 method for class 'sf'
fm_CRS(x, ..., units = NULL, oblique = NULL)

# S3 method for class 'sfc'
fm_CRS(x, ..., units = NULL, oblique = NULL)

# S3 method for class 'sfg'
fm_CRS(x, ..., units = NULL, oblique = NULL)

# S3 method for class 'fm_mesh_2d'
fm_CRS(x, ..., units = NULL, oblique = NULL)

# S3 method for class 'fm_lattice'
fm_CRS(x, ..., units = NULL, oblique = NULL)

# S3 method for class 'fm_segm'
fm_CRS(x, ..., units = NULL, oblique = NULL)

# S3 method for class 'fm_collect'
fm_CRS(x, ..., units = NULL, oblique = NULL)

# S3 method for class 'matrix'
fm_CRS(x, ..., units = NULL, oblique = NULL)

# S3 method for class 'CRS'
fm_CRS(x, ..., units = NULL, oblique = NULL)

# Default S3 method
fm_CRS(
  x,
  oblique = NULL,
  projargs = NULL,
  doCheckCRSArgs = NULL,
  args = NULL,
  SRS_string = NULL,
  ...,
  units = NULL
)

# S3 method for class 'inla.CRS'
is.na(x)

# S3 method for class 'inla.CRS'
fm_CRS(x, ..., units = NULL, oblique = NULL)
```

## Arguments

- x:

  Object to convert to CRS or to extract CRS information from.

- ...:

  Additional parameters, passed on to sub-methods.

- units:

  character; if non-NULL, `fm_length_unit()<-` is called to change the
  length units of the crs object. If `NULL` (default), the length units
  are not changed. (From version `0.3.0.9013`)

- oblique:

  Vector of length at most 4 of rotation angles (in degrees) for an
  oblique projection, all values defaulting to zero. The values indicate
  (longitude, latitude, orientation, orbit), as explained in the Details
  section for
  [`fm_crs()`](https://inlabru-org.github.io/fmesher/reference/fm_crs.md).

- projargs:

  Either 1) a projection argument string suitable as input to
  [`sp::CRS`](https://edzer.github.io/sp/reference/CRS-class.html),
  or 2) an existing `CRS` object, or 3) a shortcut reference string to a
  predefined projection; run `names(fm_wkt_predef())` for valid
  predefined projections. (projargs is a compatibility parameter that
  can be used for the default `fm_CRS()` method)

- doCheckCRSArgs:

  ignored.

- args:

  An optional list of name/value pairs to add to and/or override the
  PROJ4 arguments in `projargs`. `name=value` is converted to
  `"+name=value"`, and `name=NA` is converted to `"+name"`.

- SRS_string:

  a WKT2 string defining the coordinate system; see
  [`sp::CRS`](https://edzer.github.io/sp/reference/CRS-class.html). This
  takes precedence over `projargs`.

## Value

Either an
[`sp::CRS`](https://edzer.github.io/sp/reference/CRS-class.html) object
or an `inla.CRS` object, depending on if the coordinate reference system
described by the parameters can be expressed with a pure
[`sp::CRS`](https://edzer.github.io/sp/reference/CRS-class.html) object
or not.

An S3 `inla.CRS` object is a list, usually (but not necessarily)
containing at least one element:

- crs :

  The basic
  [`sp::CRS`](https://edzer.github.io/sp/reference/CRS-class.html)
  object

## Details

The first two elements of the `oblique` vector are the (longitude,
latitude) coordinates for the oblique centre point. The third value
(orientation) is a counterclockwise rotation angle for an observer
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

## Functions

- `is.na(fm_CRS)`: Check if a `fm_CRS` has `NA` crs information and `NA`
  obliqueness

- `is.na(inla.CRS)`: Check if a `inla.CRS` has `NA` crs information and
  `NA` obliqueness

## See also

[`fm_crs()`](https://inlabru-org.github.io/fmesher/reference/fm_crs.md),
[`sp::CRS()`](https://edzer.github.io/sp/reference/CRS-class.html),
[`fm_crs_wkt`](https://inlabru-org.github.io/fmesher/reference/fm_crs_wkt.md),
[`fm_crs_is_identical()`](https://inlabru-org.github.io/fmesher/reference/fm_crs_is_identical.md)

## Author

Finn Lindgren <Finn.Lindgren@gmail.com>

## Examples

``` r
if (fm_safe_sp()) {
  crs1 <- fm_CRS("longlat_globe")
  crs2 <- fm_CRS("lambert_globe")
  crs3 <- fm_CRS("mollweide_norm")
  crs4 <- fm_CRS("hammer_globe")
  crs5 <- fm_CRS("sphere")
  crs6 <- fm_CRS("globe")
}
```
