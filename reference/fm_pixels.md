# Generate lattice points covering a mesh

Generate `terra`, `sf`, or `sp` lattice locations

## Usage

``` r
fm_pixels(
  mesh,
  dims = c(150, 150),
  xlim = NULL,
  ylim = NULL,
  mask = TRUE,
  format = "sf",
  minimal = TRUE
)
```

## Arguments

- mesh:

  An `fm_mesh_2d` object

- dims:

  A length 2 integer vector giving the dimensions of the target lattice.

- xlim, ylim:

  Length 2 numeric vectors of x- and y- axis limits. Defaults taken from
  the range of the mesh or mask; see `minimal`.

- mask:

  If logical and TRUE, remove pixels that are outside the mesh. If
  `mask` is an `sf` or `Spatial` object, only return pixels covered by
  this object.

- format:

  character; "sf", "terra" or "sp"

- minimal:

  logical; if `TRUE` (default), the default range is determined by the
  minimum of the ranges of the mesh and mask, otherwise only the mesh.

## Value

`sf`, `SpatRaster`, or `SpatialPixelsDataFrame` covering the mesh or
mask.

## Author

Finn Lindgren <Finn.Lindgren@gmail.com>

## Examples

``` r
if (require("ggplot2", quietly = TRUE)) {
  dims <- c(50, 50)
  pxl <- fm_pixels(
    fmexample$mesh,
    dims = dims,
    mask = fmexample$boundary_sf[[1]],
    minimal = TRUE
  )
  pxl$val <- rnorm(NROW(pxl)) +
    fm_evaluate(fmexample$mesh, pxl, field = 2 * fmexample$mesh$loc[, 1])
  ggplot() +
    geom_tile(
      data = pxl,
      aes(geometry = geometry, fill = val),
      stat = "sf_coordinates"
    ) +
    geom_sf(data = fm_as_sfc(fmexample$mesh), alpha = 0.2)
}


# \donttest{
if (require("ggplot2", quietly = TRUE) &&
  require("terra", quietly = TRUE) &&
  require("tidyterra", quietly = TRUE)) {
  pxl <- fm_pixels(fmexample$mesh,
    dims = c(50, 50), mask = fmexample$boundary_sf[[1]],
    format = "terra"
  )
  pxl$val <- rnorm(NROW(pxl) * NCOL(pxl))
  pxl <-
    terra::mask(
      pxl,
      mask = pxl$.mask,
      maskvalues = c(FALSE, NA),
      updatevalue = NA
    )
  ggplot() +
    geom_spatraster(data = pxl, aes(fill = val)) +
    geom_sf(data = fm_as_sfc(fmexample$mesh), alpha = 0.2)
}
#> terra 1.9.34
#> 
#> Attaching package: ‘terra’
#> The following object is masked from ‘package:patchwork’:
#> 
#>     area
#> 
#> Attaching package: ‘tidyterra’
#> The following object is masked from ‘package:Matrix’:
#> 
#>     expand
#> The following object is masked from ‘package:stats’:
#> 
#>     filter

# }
```
