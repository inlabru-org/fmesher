# Subset integration on a mesh

Integration methods for spatial samplers on `fm_mesh_2d` meshes.

## Usage

``` r
fm_int_mesh_2d(samplers, domain, name = NULL, int.args = NULL, ...)

fm_int_mesh_2d_NULL(samplers, domain, name = NULL, int.args = NULL, ...)

# S3 method for class 'sf'
fm_int_mesh_2d(samplers, domain, name = NULL, int.args = NULL, ...)

# S3 method for class 'sfc_POINT'
fm_int_mesh_2d(
  samplers,
  domain,
  name = NULL,
  int.args = NULL,
  .weight = rep(1, NROW(samplers)),
  ...
)

# S3 method for class 'sfc_MULTIPOINT'
fm_int_mesh_2d(
  samplers,
  domain,
  name = NULL,
  int.args = NULL,
  .weight = rep(1, NROW(samplers)),
  ...
)

# S3 method for class 'sfc_LINESTRING'
fm_int_mesh_2d(
  samplers,
  domain,
  name = NULL,
  int.args = NULL,
  .weight = rep(1, NROW(samplers)),
  ...
)

# S3 method for class 'sfc_MULTILINESTRING'
fm_int_mesh_2d(
  samplers,
  domain,
  name = NULL,
  int.args = NULL,
  .weight = rep(1, NROW(samplers)),
  ...
)

# S3 method for class 'sfc_POLYGON'
fm_int_mesh_2d(
  samplers,
  domain,
  name = NULL,
  int.args = NULL,
  .weight = rep(1, NROW(samplers)),
  ...
)

# S3 method for class 'sfc_MULTIPOLYGON'
fm_int_mesh_2d(
  samplers,
  domain,
  name = NULL,
  int.args = NULL,
  .weight = rep(1, NROW(samplers)),
  ...
)

# S3 method for class 'sfc_GEOMETRY'
fm_int_mesh_2d(
  samplers,
  domain,
  name = NULL,
  int.args = NULL,
  .weight = rep(1, NROW(samplers)),
  ...
)

# S3 method for class 'Spatial'
fm_int_mesh_2d(
  samplers,
  domain,
  name = NULL,
  int.args = NULL,
  format = NULL,
  ...
)

# S3 method for class 'fm_segm'
fm_int_mesh_2d(
  samplers,
  domain,
  name = NULL,
  int.args = NULL,
  format = NULL,
  ...
)
```

## Arguments

- samplers:

  For single domain `fm_int` methods, an object specifying one or more
  subsets of the domain, and optional weighting in a `weight` variable.
  For `fm_int.list`, a list of sampling definitions, where data frame
  elements may contain information for multiple domains, in which case
  each row represent a separate tensor product integration subspace.

- domain:

  Functional space specification; single domain or a named list of
  domains

- name:

  For single-domain methods, the variable name to use for the
  integration points. Default 'x'

- int.args:

  List of arguments passed to line and integration methods.

  - `method`: "stable" (to aggregate integration weights onto mesh
    nodes) or "direct" (to construct a within triangle/segment
    integration scheme without aggregating onto mesh nodes)

  - `nsub1`, `nsub2`: integers controlling the number of internal
    integration points before aggregation. Points per triangle:
    `(nsub2+1)^2`. Points per knot segment: `nsub1`

- ...:

  Additional arguments passed on to other methods

- format:

  character; determines the output format, as either "sf" (default for
  `fm_mesh_2d` when the sampler is `NULL`), "numeric" (default for
  `fm_mesh_1d`), "bary", or "sp". When `NULL`, determined by the domain
  and sampler types.

## Value

A `list`, `sf`, or `Spatial` object with point coordinate information
and additional columns `weight` and `.block`

## Methods (by class)

- `fm_int_mesh_2d(sf)`: `sf` integration

- `fm_int_mesh_2d(sfc_POINT)`: `sfc_POINT` integration

- `fm_int_mesh_2d(sfc_MULTIPOINT)`: `sfc_MULTIPOINT` integration

- `fm_int_mesh_2d(sfc_LINESTRING)`: `sfc_LINESTRING` integration

- `fm_int_mesh_2d(sfc_MULTILINESTRING)`: `sfc_MULTILINESTRING`
  integration

- `fm_int_mesh_2d(sfc_POLYGON)`: `sfc_POLYGON` integration

- `fm_int_mesh_2d(sfc_MULTIPOLYGON)`: `sfc_MULTIPOLYGON` integration

- `fm_int_mesh_2d(sfc_GEOMETRY)`: `sfc_GEOMERY` integration

- `fm_int_mesh_2d(Spatial)`: `Spatial` integration

- `fm_int_mesh_2d(fm_segm)`: `fm_segm` integration

## Functions

- `fm_int_mesh_2d_NULL()`: Full domain integration

## Examples

``` r
str(fm_int_mesh_2d(samplers = NULL, domain = fmexample$mesh))
#> sf [279 × 4] (S3: sf/tbl_df/tbl/data.frame)
#>  $ weight       : num [1:279] 0.343 0.542 0.427 0.564 0.439 ...
#>  $ .block       : int [1:279] 1 1 1 1 1 1 1 1 1 1 ...
#>  $ geometry     :sfc_POINT of length 279; first list element:  'XYZ' num [1:3] -4.99 1.48 0
#>  $ .block_origin: int [1:279, 1] 1 1 1 1 1 1 1 1 1 1 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : NULL
#>  - attr(*, "sf_column")= chr "geometry"
#>  - attr(*, "agr")= Factor w/ 3 levels "constant","aggregate",..: NA NA NA
#>   ..- attr(*, "names")= chr [1:3] "weight" ".block" ".block_origin"
```
