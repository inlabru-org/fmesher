# Object coordinate transformation

Handle transformation of various inla objects according to coordinate
reference systems of `crs` (from
[`sf::st_crs()`](https://r-spatial.github.io/sf/reference/st_crs.html)),
`fm_crs`,
[`sp::CRS`](https://edzer.github.io/sp/reference/CRS-class.html),
`fm_CRS`, or `INLA::inla.CRS` class.

## Usage

``` r
fm_transform(x, crs, ...)

# Default S3 method
fm_transform(x, crs, ..., crs0 = NULL)

## S3 method for class 'NULL'
fm_transform(x, crs, ...)

# S3 method for class 'matrix'
fm_transform(x, crs, ..., passthrough = FALSE, crs0 = NULL)

# S3 method for class 'sf'
fm_transform(x, crs, ..., passthrough = FALSE)

# S3 method for class 'sfc'
fm_transform(x, crs, ..., passthrough = FALSE)

# S3 method for class 'sfg'
fm_transform(x, crs, ..., passthrough = FALSE)

# S3 method for class 'Spatial'
fm_transform(x, crs, ..., passthrough = FALSE)

# S3 method for class 'fm_mesh_2d'
fm_transform(x, crs = fm_crs(x), ..., passthrough = FALSE, crs0 = fm_crs(x))

# S3 method for class 'fm_collect'
fm_transform(x, crs = fm_crs(x), ..., passthrough = FALSE, crs0 = NULL)

# S3 method for class 'fm_lattice_2d'
fm_transform(x, crs = fm_crs(x), ..., passthrough = FALSE, crs0 = fm_crs(x))

# S3 method for class 'fm_segm'
fm_transform(x, crs = fm_crs(x), ..., passthrough = FALSE, crs0 = fm_crs(x))

# S3 method for class 'fm_list'
fm_transform(x, crs, ...)
```

## Arguments

- x:

  The object that should be transformed from it's current CRS to a new
  CRS

- crs:

  The target crs object

- ...:

  Potential additional arguments

- crs0:

  The source crs object for spatial classes without crs information

- passthrough:

  Default is FALSE. Setting to TRUE allows objects with no CRS
  information to be passed through without transformation. Use with
  care!

## Value

A transformed object, normally of the same class as the input object.

## See also

[`fm_CRS()`](https://inlabru-org.github.io/fmesher/reference/fm_CRS_sp.md)

## Examples

``` r
fm_transform(
  rbind(c(0, 0), c(0, 90), c(0, 91)),
  crs = fm_crs("sphere"),
  crs0 = fm_crs("longlat_norm")
)
#>              [,1] [,2] [,3]
#> [1,] 1.000000e+00    0    0
#> [2,] 6.123234e-17    0    1
#> [3,]           NA   NA   NA
```
