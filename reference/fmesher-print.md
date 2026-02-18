# Print objects

Print objects

## Usage

``` r
# S3 method for class 'fm_segm'
print(x, ..., digits = NULL, verbose = TRUE, newline = TRUE)

# S3 method for class 'fm_segm_list'
print(x, ..., digits = NULL, verbose = FALSE, newline = TRUE)

# S3 method for class 'fm_list'
print(x, ..., digits = NULL, verbose = FALSE, newline = TRUE)

# S3 method for class 'fm_mesh_2d'
print(x, ..., digits = NULL, verbose = FALSE)

# S3 method for class 'fm_mesh_3d'
print(x, ..., digits = NULL, verbose = FALSE)

# S3 method for class 'fm_mesh_1d'
print(x, ..., digits = NULL, verbose = FALSE)

# S3 method for class 'fm_bbox'
print(x, ..., digits = NULL, verbose = TRUE, newline = TRUE)

# S3 method for class 'fm_tensor'
print(x, ..., digits = NULL, verbose = FALSE)

# S3 method for class 'fm_collect'
print(x, ..., digits = NULL, verbose = FALSE)

# S3 method for class 'fm_lattice_2d'
print(x, ..., digits = NULL, verbose = FALSE)

# S3 method for class 'fm_lattice_Nd'
print(x, ..., digits = NULL, verbose = FALSE)

# S3 method for class 'fm_crs'
print(x, ...)

# S3 method for class 'fm_CRS'
print(x, ...)
```

## Arguments

- x:

  an object used to select a method.

- ...:

  further arguments passed to or from other methods.

- digits:

  a positive integer indicating how many significant digits are to be
  used for numeric and complex x. The default, NULL, uses
  `getOption("digits")`.

- verbose:

  logical

- newline:

  logical; if `TRUE` (default), end the printing with `\n`

## Value

The input object `x`

## Examples

``` r
fm_bbox(matrix(1:6, 3, 2))
#> Bounding box: (1,3) x (4,6)
print(fm_bbox(matrix(1:6, 3, 2)), verbose = FALSE)
#> (1,3) x (4,6)

print(fmexample$mesh)
#> fm_mesh_2d object:
#>   Manifold:  R2
#>   V / E / T: 292 / 838 / 547
#>   Euler char.:   1
#>   Constraints:   Boundary: 35 boundary edges (1 group: 1), Interior: 55 interior edges (1 group: 1)
#>   Bounding box: (-5.345477, 4.083580) x (-3.997839, 5.415519)
#>   Basis d.o.f.:  292
print(fmexample$boundary_fm)
#> 43 boundary edges (1 group: 1)
#> 34 boundary edges (1 group: 1)

print(fm_mesh_1d(c(1, 2, 3, 5, 7), degree = 2))
#> fm_mesh_1d object:
#>   Manifold:  R1
#>   #{knots}:  5
#>   Interval:  (1, 7)
#>   Boundary:  (neumann, neumann)
#>   B-spline degree:   2
#>   Basis d.o.f.:  4
```
