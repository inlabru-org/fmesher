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
#> # A tibble: 6 × 6
#>   weight .block .block_origin[,1] .vertex bary$index $where[,1] $[,2] loc[,"X"]
#>    <dbl>  <int>             <int>   <int>      <int>      <dbl> <dbl>     <dbl>
#> 1 0.0309      1                 1     103         66          1     0    -1.95 
#> 2 0.124       1                 1     107         13          0     0    -1.09 
#> 3 0.465       1                 1     111         50          0     0    -2.41 
#> 4 0.561       1                 1     118         59          0     0    -0.635
#> 5 0.338       1                 1     119         70          1     0    -1.45 
#> 6 0.255       1                 1     123         43          0     1     0.366
#> # ℹ 2 more variables: bary$where[3] <dbl>, loc[2:3] <dbl>
head(fm_vertex_projection(fmexample$loc_sf, fmexample$mesh))
#> Simple feature collection with 6 features and 5 fields
#> Geometry type: POINT
#> Dimension:     XY
#> Bounding box:  xmin: -2.412672 ymin: -0.5847901 xmax: 0.3656907 ymax: 2.127547
#> CRS:           NA
#> # A tibble: 6 × 6
#>   weight .block .block_origin[,1] .vertex bary$index               geometry
#>    <dbl>  <int>             <int>   <int>      <int>                <POINT>
#> 1 0.0309      1                 1     103         66  (-1.954974 0.0878308)
#> 2 0.124       1                 1     107         13   (-1.094063 2.127547)
#> 3 0.465       1                 1     111         50  (-2.412672 0.3138631)
#> 4 0.561       1                 1     118         59 (-0.635088 -0.3472952)
#> 5 0.338       1                 1     119         70  (-1.44505 -0.5847901)
#> 6 0.255       1                 1     123         43  (0.3656907 0.1023152)
#> # ℹ 1 more variable: bary$where <dbl[,3]>
```
