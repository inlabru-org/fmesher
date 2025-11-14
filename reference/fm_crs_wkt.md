# Handling CRS/WKT

Get and set CRS object or WKT string properties.

## Usage

``` r
fm_wkt_is_geocent(wkt)

fm_crs_is_geocent(crs)

fm_wkt_get_ellipsoid_radius(wkt)

fm_crs_get_ellipsoid_radius(crs)

fm_ellipsoid_radius(x)

# Default S3 method
fm_ellipsoid_radius(x)

# S3 method for class 'character'
fm_ellipsoid_radius(x)

fm_wkt_set_ellipsoid_radius(wkt, radius)

fm_ellipsoid_radius(x) <- value

# S3 method for class 'character'
fm_ellipsoid_radius(x) <- value

# S3 method for class 'CRS'
fm_ellipsoid_radius(x) <- value

# S3 method for class 'fm_CRS'
fm_ellipsoid_radius(x) <- value

# S3 method for class 'crs'
fm_ellipsoid_radius(x) <- value

# S3 method for class 'fm_crs'
fm_ellipsoid_radius(x) <- value

fm_crs_set_ellipsoid_radius(crs, radius)

fm_wkt_unit_params()

fm_wkt_get_lengthunit(wkt)

fm_wkt_set_lengthunit(wkt, unit, params = NULL)

fm_crs_get_lengthunit(crs)

fm_crs_set_lengthunit(crs, unit)

fm_length_unit(x)

# Default S3 method
fm_length_unit(x)

# S3 method for class 'character'
fm_length_unit(x)

fm_length_unit(x) <- value

# S3 method for class 'character'
fm_length_unit(x) <- value

# S3 method for class 'CRS'
fm_length_unit(x) <- value

# S3 method for class 'fm_CRS'
fm_length_unit(x) <- value

# S3 method for class 'crs'
fm_length_unit(x) <- value

# S3 method for class 'fm_crs'
fm_length_unit(x) <- value

fm_wkt(crs)

fm_proj4string(crs)

fm_wkt_tree_projection_type(wt)

fm_wkt_projection_type(wkt)

fm_crs_projection_type(crs)

fm_crs_bounds(crs, warn.unknown = FALSE)

# S3 method for class 'inla.CRS'
fm_ellipsoid_radius(x) <- value

# S3 method for class 'inla.CRS'
fm_length_unit(x) <- value
```

## Arguments

- wkt:

  A WKT2 character string

- crs:

  An
  [`sf::crs`](https://r-spatial.github.io/sf/reference/coerce-methods.html),
  [`sp::CRS`](https://edzer.github.io/sp/reference/CRS-class.html),
  `fm_crs` or `inla.CRS` object

- x:

  crs object to extract value from or assign values in

- radius:

  numeric; The new radius value

- value:

  Value to assign

- unit:

  character, name of a unit. Supported names are "metre", "kilometre",
  and the aliases "meter", "m", International metre", "kilometer", and
  "km", as defined by `fm_wkt_unit_params` or the `params` argument.

- params:

  Length unit definitions, in the list format produced by
  `fm_wkt_unit_params()`, Default: NULL, which invokes
  `fm_wkt_unit_params()`

- wt:

  A parsed wkt tree, see
  [`fm_wkt_as_wkt_tree()`](https://inlabru-org.github.io/fmesher/reference/wkt_tree.md)

- warn.unknown:

  logical, default `FALSE`. Produce warning if the shape of the
  projection bounds is unknown.

## Value

For `fm_wkt_unit_params`, a list of named unit definitions

For `fm_wkt_get_lengthunit`, a list of length units used in the wkt
string, excluding the ellipsoid radius unit.

For `fm_wkt_set_lengthunit`, a WKT2 string with altered length units.
Note that the length unit for the ellipsoid radius is unchanged.

For `fm_crs_get_lengthunit`, a list of length units used in the wkt
string, excluding the ellipsoid radius unit.

For `fm_length_unit<-`, a crs object with altered length units. Note
that the length unit for the ellipsoid radius is unchanged.

## Functions

- `fm_wkt()`: Returns a WKT2 string, for any input supported by
  [`fm_crs()`](https://inlabru-org.github.io/fmesher/reference/fm_crs.md).

- `fm_proj4string()`: Returns a proj4 string, for any input supported by
  [`fm_crs()`](https://inlabru-org.github.io/fmesher/reference/fm_crs.md).

- `fm_wkt_tree_projection_type()`: Returns "longlat", "lambert",
  "mollweide", "hammer", "tmerc", or `NULL`

- `fm_wkt_projection_type()`: See `fm_wkt_tree_projection_type`

- `fm_crs_projection_type()`: See `fm_wkt_tree_projection_type`

- `fm_crs_bounds()`: Returns bounds information for a projection, as a
  list with elements `type` ("rectangle" or "ellipse"), `xlim`, `ylim`,
  and `polygon`.

## See also

[`fm_crs()`](https://inlabru-org.github.io/fmesher/reference/fm_crs.md)

## Author

Finn Lindgren <Finn.Lindgren@gmail.com>

## Examples

``` r
c1 <- fm_crs("globe")
fm_length_unit(c1)
#> $kilometre
#> $kilometre[[1]]
#> [1] "\"kilometre\""
#> 
#> $kilometre[[2]]
#> [1] "1000"
#> 
#> $kilometre[[3]]
#> $kilometre[[3]]$label
#> [1] "ID"
#> 
#> $kilometre[[3]]$params
#> $kilometre[[3]]$params[[1]]
#> [1] "\"EPSG\""
#> 
#> $kilometre[[3]]$params[[2]]
#> [1] "9036"
#> 
#> 
#> 
#> 
fm_length_unit(c1) <- "m"
fm_length_unit(c1)
#> $metre
#> $metre[[1]]
#> [1] "\"metre\""
#> 
#> $metre[[2]]
#> [1] "1"
#> 
#> $metre[[3]]
#> $metre[[3]]$label
#> [1] "ID"
#> 
#> $metre[[3]]$params
#> $metre[[3]]$params[[1]]
#> [1] "\"EPSG\""
#> 
#> $metre[[3]]$params[[2]]
#> [1] "9001"
#> 
#> 
#> 
#> 
```
