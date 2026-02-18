# Add sp data to fmexample

Adds `loc_sp` and `boundary_sp` to
[fmexample](https://inlabru-org.github.io/fmesher/reference/fmexample.md)
for use in `sp` related code examples and tests.

## Usage

``` r
fmexample_sp()
```

## Value

Returns a copy of
[fmexample](https://inlabru-org.github.io/fmesher/reference/fmexample.md)
with `loc_sp` (`SpatialPoints`) and `boundary_sp` (`SpatialPolygons`)
added.

## Examples

``` r
if (fm_safe_sp()) {
  fmexample_sp()
}
#> Error in StopZ(zm): sp supports Z dimension only for POINT and MULTIPOINT.
#> use `st_zm(...)` to coerce to XY dimensions
```
