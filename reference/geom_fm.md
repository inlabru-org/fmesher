# ggplot2 geomes for fmesher related objects

**\[experimental\]**

`geom_fm` is a generic function for generating geomes from various kinds
of `fmesher` objects, e.g. `fm_segm` and `fm_mesh_2d`. The function
invokes particular methods which depend on the
[class](https://rdrr.io/r/base/class.html) of the `data` argument.
Requires the `ggplot2` package.

Note: `geom_fm` is not yet a "proper" `ggplot2` geom method; the
interface may therefore change in the future.

## Usage

``` r
geom_fm(mapping = NULL, data = NULL, ...)

# S3 method for class 'fm_mesh_2d'
geom_fm(
  mapping = NULL,
  data = NULL,
  ...,
  mappings = NULL,
  defs = NULL,
  crs = NULL
)

# S3 method for class 'fm_segm'
geom_fm(mapping = NULL, data = NULL, ..., crs = NULL)

# S3 method for class 'fm_mesh_1d'
geom_fm(
  mapping = NULL,
  data = NULL,
  ...,
  mappings = NULL,
  defs = NULL,
  xlim = NULL,
  basis = TRUE,
  knots = TRUE,
  derivatives = FALSE,
  weights = NULL
)
```

## Arguments

- mapping:

  [`ggplot2::aes()`](https://ggplot2.tidyverse.org/reference/aes.html)
  mapping information.

- data:

  an object for which to generate a geom.

- ...:

  Arguments passed on to the geom method.

- mappings, defs:

  optional lists of `aes` mappings and non-`aes` settings. For
  `fm_mesh_2d`, the non-triangle parts of the mesh, named "int" for
  interior constraint edges, "bnd" for boundary edges, and "loc" for the
  vertices. For `fm_mesh_1d`, the elements are "knots" and "fun".

- crs:

  Optional crs to transform the object to before plotting.

- xlim:

  numeric 2-vector; specifies the interval for which to compute
  functions. Default is `data$interval`

- basis:

  logical; if `TRUE` (default), show the spline basis functions

- knots:

  logical; if `TRUE` (default), show the spline knot locations

- derivatives:

  logical; if `TRUE` (not default), draw first order derivatives instead
  of function values

- weights:

  numeric vector; if provided, draw weighted basis functions and the
  resulting weighted sum.

## Value

A combination of `ggplot2` geoms.

## Methods (by class)

- `geom_fm(fm_mesh_2d)`: Converts an
  [`fm_mesh_2d()`](https://inlabru-org.github.io/fmesher/reference/fm_mesh_2d.md)
  object to `sf` with
  [`fm_as_sfc()`](https://inlabru-org.github.io/fmesher/reference/fm_as_sfc.md)
  and uses `geom_sf` to visualize the triangles and edges.

  The mesh vertices are only plotted if `mappings$loc` or `defs$loc` is
  non-`NULL`, e.g. `defs = list(loc = list())`. Default argument
  settings:

      ... = linewidth = 0.25, color = "grey" # default for triangle mapping
      defs = list(
        int = list(linewidth = 0.5, color = "blue"),
        bnd = list(linewidth = 1, color = "black", alpha = 0),
        loc = list(size = 1, color = "red")
      )

- `geom_fm(fm_segm)`: Converts an
  [`fm_segm()`](https://inlabru-org.github.io/fmesher/reference/fm_segm.md)
  object to `sf` with
  [`fm_as_sfc()`](https://inlabru-org.github.io/fmesher/reference/fm_as_sfc.md)
  and uses `geom_sf` to visualize it.

- `geom_fm(fm_mesh_1d)`: Evaluates and plots the basis functions defined
  by an
  [`fm_mesh_1d()`](https://inlabru-org.github.io/fmesher/reference/fm_mesh_1d.md)
  object.

## Examples

``` r
ggplot() +
  geom_fm(data = fmexample$mesh)

m <- fm_mesh_2d(
  cbind(10, 20),
  boundary = fm_extensions(cbind(10, 20), c(25, 65)),
  max.edge = c(4, 10),
  crs = fm_crs("+proj=longlat")
)
ggplot() +
  geom_fm(data = m)
#> st_as_s2(): dropping Z and/or M coordinate

ggplot() +
  geom_fm(data = m, defs = list(loc = list()))
#> st_as_s2(): dropping Z and/or M coordinate

ggplot() +
  geom_fm(data = m, crs = fm_crs("epsg:27700"))

# \donttest{
# Compute a mesh vertex based function on a different grid
px <- fm_pixels(
  fm_transform(m, fm_crs("mollweide_globe")),
  dims = c(50, 50) # Speed up the example by lowering the resolution
)
px$fun <- fm_evaluate(m,
  loc = px,
  field = sin(m$loc[, 1] / 5) * sin(m$loc[, 2] / 5)
)
ggplot() +
  geom_tile(aes(geometry = geometry, fill = fun),
    data = px,
    stat = "sf_coordinates"
  ) +
  geom_fm(
    data = m, alpha = 0.2, linewidth = 0.05,
    crs = fm_crs("mollweide_globe")
  )

# }
m1 <- fm_segm(rbind(c(1, 2), c(4, 3), c(2, 4)), is.bnd = TRUE)
m2 <- fm_segm(rbind(c(2, 2), c(3, 4), c(2, 3)), is.bnd = FALSE)
ggplot() +
  geom_fm(data = m1) +
  geom_fm(data = m2)

m <- fm_mesh_1d(
  c(1, 2, 3, 5, 7),
  boundary = c("dirichlet", "neumann"),
  degree = 2
)
ggplot() +
  geom_fm(data = m)
```
