# (Blockwise) cross product of integration points

Calculates the groupwise cross product of integration points in
different dimensions and multiplies their weights accordingly. If the
object defining points in a particular dimension has no weights attached
to it all weights are assumed to be 1.

## Usage

``` r
fm_cprod(..., na.rm = NULL, .blockwise = FALSE)
```

## Arguments

- ...:

  `tibble`, `data.frame`, `sf`, or `SpatialPointsDataFrame` objects,
  each one usually obtained by a call to an
  [`fm_int()`](https://inlabru-org.github.io/fmesher/reference/fm_int.md)
  method.

- na.rm:

  logical; if `TRUE`, the rows with weight `NA` from the non-overlapping
  full_join will be removed; if `FALSE`, set the undefined weights to
  `NA`. If `NULL` (default), act as `TRUE`, but warn if any elements
  needed removing.

- .blockwise:

  logical; if `FALSE`, computes full tensor product integration. If
  `TRUE`, computes within-block tensor product integration (used
  internally by
  [`fm_int()`](https://inlabru-org.github.io/fmesher/reference/fm_int.md)).
  Default `FALSE`

## Value

A `data.frame`, `sf`, or `SpatialPointsDataFrame` of multidimensional
integration points and their weights

## Examples

``` r
if (require("ggplot2")) {
  # Create integration points in dimension 'myDim' and 'myDiscreteDim'
  ips1 <- fm_int(fm_mesh_1d(1:20),
    rbind(c(0, 3), c(3, 8)),
    name = "myDim"
  )
  ips2 <- fm_int(domain = c(1, 2, 4), name = "myDiscreteDim")

  # Calculate the cross product
  ips <- fm_cprod(ips1, ips2)

  # Plot the integration points
  ggplot(ips) +
    geom_point(aes(myDim, myDiscreteDim, size = weight), stroke = 0) +
    scale_size_area()
}

```
