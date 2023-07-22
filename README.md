
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fmesher:

<!-- badges: start -->

[![inlabru status
badge](https://inlabru-org.r-universe.dev/badges/fmesher)](https://inlabru-org.r-universe.dev)
[![R build
status](https://github.com/inlabru-org/fmesher/workflows/R-CMD-check/badge.svg)](https://github.com/inlabru-org/fmesher/actions)
[![test-coverage](https://github.com/inlabru-org/fmesher/workflows/test-coverage/badge.svg)](https://github.com/inlabru-org/fmesher/actions)
[![Codecov test
coverage](https://codecov.io/gh/inlabru-org/fmesher/branch/devel/graph/badge.svg)](https://codecov.io/gh/inlabru-org/fmesher?branch=devel)
<!-- badges: end -->

fmesher: Triangle Meshes and Other Geometry Tools

## Installation

You can install the latest stable release of fmesher from
[GitHub](https://github.com/):

``` r
# install.packages("remotes")
remotes::install_github("inlabru-org/fmesher", ref = "stable")
```

And the development version with:

``` r
# install.packages("remotes")
remotes::install_github("inlabru-org/fmesher")
```

## Online documentation

<https://inlabru-org.github.io/fmesher/>

## Example

``` r
library(fmesher)
#> The legacy packages maptools, rgdal, and rgeos, underpinning the sp package,
#> which was just loaded, will retire in October 2023.
#> Please refer to R-spatial evolution reports for details, especially
#> https://r-spatial.org/r/2023/05/15/evolution4.html.
#> It may be desirable to make the sf package available;
#> package maintainers should consider adding sf to Suggests:.
#> The sp package is now running under evolution status 2
#>      (status 2 uses the sf package in place of rgdal)
# longlat for a spherical version of the Earth
print(fm_crs("longlat_globe"))
#> Coordinate Reference System:
#>   User input: +proj=longlat +ellps=sphere +no_defs 
#>   wkt:
#> GEOGCRS["unknown",
#>     DATUM["Unknown based on Normal Sphere (r=6370997) ellipsoid",
#>         ELLIPSOID["Normal Sphere (r=6370997)",6370997,0,
#>             LENGTHUNIT["metre",1,
#>                 ID["EPSG",9001]]]],
#>     PRIMEM["Greenwich",0,
#>         ANGLEUNIT["degree",0.0174532925199433],
#>         ID["EPSG",8901]],
#>     CS[ellipsoidal,2],
#>         AXIS["longitude",east,
#>             ORDER[1],
#>             ANGLEUNIT["degree",0.0174532925199433,
#>                 ID["EPSG",9122]]],
#>         AXIS["latitude",north,
#>             ORDER[2],
#>             ANGLEUNIT["degree",0.0174532925199433,
#>                 ID["EPSG",9122]]]]

# longlat for a sphere of radius 1m
print(fm_crs("longlat_norm"))
#> Coordinate Reference System:
#>   User input: +proj=longlat +R=1 +no_defs 
#>   wkt:
#> GEOGCRS["unknown",
#>     DATUM["unknown",
#>         ELLIPSOID["unknown",1,0,
#>             LENGTHUNIT["metre",1,
#>                 ID["EPSG",9001]]]],
#>     PRIMEM["Reference meridian",0,
#>         ANGLEUNIT["degree",0.0174532925199433,
#>             ID["EPSG",9122]]],
#>     CS[ellipsoidal,2],
#>         AXIS["longitude",east,
#>             ORDER[1],
#>             ANGLEUNIT["degree",0.0174532925199433,
#>                 ID["EPSG",9122]]],
#>         AXIS["latitude",north,
#>             ORDER[2],
#>             ANGLEUNIT["degree",0.0174532925199433,
#>                 ID["EPSG",9122]]]]

# A sphere of radius 1m
print(fm_crs("sphere"))
#> Coordinate Reference System:
#>   User input: +proj=geocent +R=1 +units=m +no_defs 
#>   wkt:
#> GEODCRS["unknown",
#>     DATUM["unknown",
#>         ELLIPSOID["unknown",1,0,
#>             LENGTHUNIT["metre",1,
#>                 ID["EPSG",9001]]]],
#>     PRIMEM["Reference meridian",0,
#>         ANGLEUNIT["degree",0.0174532925199433,
#>             ID["EPSG",9122]]],
#>     CS[Cartesian,3],
#>         AXIS["(X)",geocentricX,
#>             ORDER[1],
#>             LENGTHUNIT["metre",1,
#>                 ID["EPSG",9001]]],
#>         AXIS["(Y)",geocentricY,
#>             ORDER[2],
#>             LENGTHUNIT["metre",1,
#>                 ID["EPSG",9001]]],
#>         AXIS["(Z)",geocentricZ,
#>             ORDER[3],
#>             LENGTHUNIT["metre",1,
#>                 ID["EPSG",9001]]]]
```
