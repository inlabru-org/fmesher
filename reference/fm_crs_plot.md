# Plot CRS and fm_crs objects

**\[experimental\]** Plot the outline of a `crs` or
[`fm_crs()`](https://inlabru-org.github.io/fmesher/reference/fm_crs.md)
projection, with optional graticules (transformed parallels and
meridians) and Tissot indicatrices.

## Usage

``` r
fm_crs_plot(
  x,
  xlim = NULL,
  ylim = NULL,
  outline = TRUE,
  graticule = c(15, 15, 45),
  tissot = c(30, 30, 30),
  asp = 1,
  add = FALSE,
  eps = 0.05,
  ...
)

fm_crs_graticule(
  x,
  by = c(15, 15, 45),
  add = FALSE,
  do.plot = TRUE,
  eps = 0.05,
  ...
)

fm_crs_tissot(
  x,
  by = c(30, 30, 30),
  add = FALSE,
  do.plot = TRUE,
  eps = 0.05,
  diff.eps = 0.01,
  ...
)
```

## Arguments

- x:

  A `crs` or
  [`fm_crs()`](https://inlabru-org.github.io/fmesher/reference/fm_crs.md)
  object.

- xlim:

  Optional x-axis limits.

- ylim:

  Optional y-axis limits.

- outline:

  Logical, if `TRUE`, draw the outline of the projection.

- graticule:

  Vector of length at most 3, to plot meridians with spacing
  `graticule[1]` degrees and parallels with spacing `graticule[2]`
  degrees. `graticule[3]` optionally specifies the spacing above and
  below the first and last parallel. When `graticule[1]==0` no meridians
  are drawn, and when `graticule[2]==0` no parallels are drawn. Use
  `graticule=NULL` to skip drawing a graticule.

- tissot:

  Vector of length at most 3, to plot Tissot's indicatrices with spacing
  `tissot[1]` degrees and parallels with spacing `tissot[2]` degrees.
  `tissot[3]` specifices a scaling factor. Use `tissot=NULL` to skip
  drawing a Tissot's indicatrices.

- asp:

  The aspect ratio for the plot, default 1.

- add:

  If `TRUE`, add the projecton plot to an existing plot.

- eps:

  Clipping tolerance for rudimentary boundary clipping

- ...:

  Additional arguments passed on to the internal calls to `plot` and
  `lines`.

- by:

  The spacing between `(long, lat, long_at_poles)`
  graticules/indicatrices, see the `graticule` and `tissot` arguments.

- do.plot:

  logical; If TRUE, do plotting

- diff.eps:

  Pre-scaling

## Value

`NULL`, invisibly

## Functions

- `fm_crs_graticule()`: **\[experimental\]** Constructs graticule
  information for a given `CRS` or
  [`fm_crs()`](https://inlabru-org.github.io/fmesher/reference/fm_crs.md)
  and optionally plots the graticules. Returns a list with two elements,
  `meridians` and `parallels`, which are `SpatialLines` objects.

- `fm_crs_tissot()`: **\[experimental\]** Constructs Tissot indicatrix
  information for a given `CRS` or
  [`fm_crs()`](https://inlabru-org.github.io/fmesher/reference/fm_crs.md)
  and optionally plots the indicatrices. Returns a list with one
  element, `tissot`, which is a `SpatialLines` object.

## See also

[`fm_crs()`](https://inlabru-org.github.io/fmesher/reference/fm_crs.md)

## Author

Finn Lindgren <Finn.Lindgren@gmail.com>

## Examples

``` r
# \donttest{
if (require("sf") && require("sp")) {
  for (projtype in c(
    "longlat_norm",
    "lambert_norm",
    "mollweide_norm",
    "hammer_norm"
  )) {
    fm_crs_plot(fm_crs(projtype), main = projtype)
  }
}
#> Loading required package: sf
#> Linking to GEOS 3.12.1, GDAL 3.8.4, PROJ 9.4.0; sf_use_s2() is TRUE
#> Loading required package: sp





if (require("sf") && require("sp")) {
  oblique <- c(0, 45, 45, 0)
  for (projtype in c(
    "longlat_norm",
    "lambert_norm",
    "mollweide_norm",
    "hammer_norm"
  )) {
    fm_crs_plot(
      fm_crs(projtype, oblique = oblique),
      main = paste("oblique", projtype)
    )
  }
}




# }
```
