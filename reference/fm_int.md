# Multi-domain integration

Construct integration points on tensor product spaces

## Usage

``` r
fm_int(domain, samplers = NULL, ...)

# S3 method for class 'list'
fm_int(domain, samplers = NULL, ...)

# S3 method for class 'numeric'
fm_int(domain, samplers = NULL, name = "x", ...)

# S3 method for class 'character'
fm_int(domain, samplers = NULL, name = "x", ...)

# S3 method for class 'factor'
fm_int(domain, samplers = NULL, name = "x", ...)

# S3 method for class 'SpatRaster'
fm_int(domain, samplers = NULL, name = "x", ...)

# S3 method for class 'fm_lattice_2d'
fm_int(domain, samplers = NULL, name = "x", ...)

# S3 method for class 'fm_mesh_1d'
fm_int(
  domain,
  samplers = NULL,
  name = "x",
  int.args = NULL,
  format = NULL,
  ...
)

# S3 method for class 'fm_mesh_2d'
fm_int(
  domain,
  samplers = NULL,
  name = NULL,
  int.args = NULL,
  format = NULL,
  ...
)
```

## Arguments

- domain:

  Functional space specification; single domain or a named list of
  domains

- samplers:

  For single domain `fm_int` methods, an object specifying one or more
  subsets of the domain, and optional weighting in a `weight` variable.
  For `fm_int.list`, a list of sampling definitions, where data frame
  elements may contain information for multiple domains, in which case
  each row represent a separate tensor product integration subspace.

- ...:

  Additional arguments passed on to other methods

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

- format:

  character; determines the output format, as either "sf" (default for
  `fm_mesh_2d` when the sampler is `NULL`), "numeric" (default for
  `fm_mesh_1d`), "bary", or "sp". When `NULL`, determined by the domain
  and sampler types.

## Value

A `tibble`, `sf`, or `SpatialPointsDataFrame` of 1D and 2D integration
points, including a `weight` column, a`.block` column, and a matrix
column `.block_origin`. The `.block` column is used to identify the
integration blocks defined by the samplers. The `.block_origin` collects
the original subdomain block information for tensor product blocks.

## Methods (by class)

- `fm_int(list)`: Multi-domain integration

- `fm_int(numeric)`: Discrete double or integer space integration

- `fm_int(character)`: Discrete character space integration

- `fm_int(factor)`: Discrete factor space integration

- `fm_int(SpatRaster)`: `SpatRaster` integration. Not yet implemented.

- `fm_int(fm_lattice_2d)`: `fm_lattice_2d` integration. Not yet
  implemented.

- `fm_int(fm_mesh_1d)`: `fm_mesh_1d` integration. Supported samplers:

  - `NULL` for integration over the entire domain;

  - A vector defining points for summation (up to `0.5.0`, length 2
    vectors were interpreted as intervals. From 0.6.0 intervals must be
    specified as rows of a 2-column matrix);

  - A 2-column matrix with a single interval in each row;

  - A list of such vectors or matrices

  - A tibble with a named column containing a vector/matrix/list as
    above, and optionally a `weight` column.

- `fm_int(fm_mesh_2d)`: `fm_mesh_2d` integration. Any sampler class with
  an associated
  [`fm_int_mesh_2d()`](https://inlabru-org.github.io/fmesher/reference/fm_int_mesh_2d.md)
  method is supported.

## Examples

``` r
# Integration on the interval (2, 3.5) with Simpson's rule
ips <- fm_int(fm_mesh_1d(0:4), samplers = cbind(2, 3.5))
plot(ips$x, ips$weight)


# Create integration points for the two intervals [0,3] and [5,10]
ips <- fm_int(
  fm_mesh_1d(0:10),
  rbind(c(0, 3), c(5, 10))
)
plot(ips$x, ips$weight)


# Convert a 1D mesh into integration points
mesh <- fm_mesh_1d(seq(0, 10, by = 1))
ips <- fm_int(mesh, name = "time")
plot(ips$time, ips$weight)


if (require("ggplot2", quietly = TRUE)) {
  #' Integrate on a 2D mesh with polygon boundary subset
  ips <- fm_int(fmexample$mesh, fmexample$boundary_sf[[1]])
  ggplot() +
    geom_sf(data = fm_as_sfc(fmexample$mesh, multi = TRUE), alpha = 0.5) +
    geom_sf(data = fmexample$boundary_sf[[1]], fill = "red", alpha = 0.5) +
    geom_sf(data = ips, aes(size = weight)) +
    scale_size_area()
}


# Individual sampling points:
(ips <- fm_int(0:10, c(0, 3, 5, 6, 10)))
#> # A tibble: 5 × 4
#>       x weight .block .block_origin[,"x"]
#>   <int>  <dbl>  <int>               <int>
#> 1     0      1      1                   1
#> 2     3      1      2                   2
#> 3     5      1      3                   3
#> 4     6      1      4                   4
#> 5    10      1      5                   5
# Sampling blocks:
(ips <- fm_int(0:10, list(c(0, 3), c(5, 6, 10))))
#> # A tibble: 5 × 4
#>       x weight .block .block_origin[,"x"]
#>   <int>  <dbl>  <int>               <int>
#> 1     0      1      1                   1
#> 2     3      1      1                   1
#> 3     5      1      2                   2
#> 4     6      1      2                   2
#> 5    10      1      2                   2

# Continuous integration on intervals
ips <- fm_int(
  fm_mesh_1d(0:10, boundary = "cyclic"),
  rbind(c(0, 3), c(5, 10))
)
plot(ips$x, ips$weight)

```
