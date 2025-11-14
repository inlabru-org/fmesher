# Project integration points to mesh vertices

Compute information for assigning points to the vertices of the covering
triangle

## Usage

``` r
fm_vertex_projection(points, mesh)
```

## Arguments

- points:

  A `SpatialPointsDataFrame`, `sf`, `tibble`, or `list` object

- mesh:

  An `fm_mesh_2d` object

## Value

`SpatialPointsDataFrame`, `sf`, `tibble`, or `list` of mesh vertices
with projected data attached

## Examples

``` r
head(fm_vertex_projection(list(loc = fmexample$loc), fmexample$mesh))
#> # A tibble: 6 × 5
#>   weight .block .block_origin[,1] .vertex loc[,"X"] [,"Y"] [,"Z"]
#>    <dbl>  <int>             <int>   <int>     <dbl>  <dbl>  <dbl>
#> 1 0.306       1                 1      84     0.294 -1.26       0
#> 2 0.0477      1                 1      85    -1.22   2.45       0
#> 3 0.0831      1                 1      88    -2.50  -0.254      0
#> 4 1.11        1                 1      89    -0.574 -1.00       0
#> 5 0.0253      1                 1      97     0.419 -0.434      0
#> 6 0.0355      1                 1     103    -1.49  -0.250      0
head(fm_vertex_projection(fmexample$loc_sf, fmexample$mesh))
#> Simple feature collection with 6 features and 4 fields
#> Geometry type: POINT
#> Dimension:     XY
#> Bounding box:  xmin: -2.496261 ymin: -1.262638 xmax: 0.4193104 ymax: 2.448375
#> CRS:           NA
#> # A tibble: 6 × 5
#>   weight .block .block_origin[,1] .vertex               geometry
#>    <dbl>  <int>             <int>   <int>                <POINT>
#> 1 0.306       1                 1      84  (0.2936639 -1.262638)
#> 2 0.0477      1                 1      85   (-1.220425 2.448375)
#> 3 0.0831      1                 1      88 (-2.496261 -0.2538804)
#> 4 1.11        1                 1      89 (-0.5743072 -1.002429)
#> 5 0.0253      1                 1      97  (0.4193104 -0.433815)
#> 6 0.0355      1                 1     103 (-1.486897 -0.2497477)
```
