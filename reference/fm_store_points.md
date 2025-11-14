# Store points in different formats

Convert a matrix of points into different formats.

## Usage

``` r
fm_store_points(loc, crs = NULL, info = NULL, format = NULL)
```

## Arguments

- loc:

  a coordinate matrix

- crs:

  CRS information to associate with the coordinates

- info:

  An optional data.frame of additional data

- format:

  character; `"sf"`, `"df"`, `"sp"`

## Value

An `sf`, `data.frame`, or `SpatialPointsDataFrame` object, with optional
added information.

## Examples

``` r
fm_store_points(fmexample$loc, format = "sf")
#> Simple feature collection with 10 features and 0 fields
#> Geometry type: POINT
#> Dimension:     XY
#> Bounding box:  xmin: -2.345698 ymin: -0.9983864 xmax: 1.084441 ymax: 2.415835
#> CRS:           NA
#>                         geometry
#> 1   POINT (-1.207066 -0.4771927)
#> 2   POINT (0.2774292 -0.9983864)
#> 3    POINT (1.084441 -0.7762539)
#> 4   POINT (-2.345698 0.06445882)
#> 5    POINT (0.4291247 0.9594941)
#> 6   POINT (0.5060559 -0.1102855)
#> 7    POINT (-0.57474 -0.5110095)
#> 8  POINT (-0.5466319 -0.9111954)
#> 9   POINT (-0.564452 -0.8371717)
#> 10   POINT (-0.8900378 2.415835)
```
