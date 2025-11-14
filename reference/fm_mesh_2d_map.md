# Special coordinate mappings for `fm_mesh_2d` projections.

Calculates coordinate mappings for spherical `fm_mesh_2d` projections.
This is an internal function not intended for general use.

## Usage

``` r
fm_mesh_2d_map(loc, projection = NULL, inverse = TRUE)

fm_mesh_2d_map_lim(loc = NULL, projection = NULL)
```

## Arguments

- loc:

  Coordinates to be mapped.

- projection:

  The projection type. One of `NULL`, "default", "longlat",
  "longsinlat", or "mollweide".

- inverse:

  If `TRUE`, `loc` are map coordinates and coordinates in the spherical
  domain are calculated. If `FALSE`, `loc` are coordinates in the
  spherical domain and the forward map projection is calculated.
  Default: `TRUE`

## Value

For `fm_mesh_2d_map_lim`, a list:

- xlim :

  X axis limits in the map domain

- ylim :

  Y axis limits in the map domain

No attempt is made to find minimal limits for partial spherical domains.

## Functions

- `fm_mesh_2d_map_lim()`: Projection extent limit calculations

## See also

[`fm_evaluator()`](https://inlabru-org.github.io/fmesher/reference/fm_evaluate.md)

## Author

Finn Lindgren <Finn.Lindgren@gmail.com>

## Examples

``` r
(loc <- fm_mesh_2d_map(cbind(20, 10), "longlat"))
#>           [,1]      [,2]      [,3]
#> [1,] 0.9254166 0.3368241 0.1736482
fm_mesh_2d_map(loc, "longlat", inverse = FALSE)
#>      [,1] [,2]
#> [1,]   20   10
```
