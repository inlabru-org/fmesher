# Compute barycentric coordinates

Identify knot intervals or triangles and compute barycentric coordinates

## Usage

``` r
fm_bary(...)

# S3 method for class 'fm_bary'
fm_bary(bary, ..., extra_class = NULL)

# S3 method for class 'list'
fm_bary(bary, ..., extra_class = NULL)

# S3 method for class 'tbl_df'
fm_bary(bary, ..., extra_class = NULL)

# S3 method for class 'fm_mesh_1d'
fm_bary(mesh, loc, method = c("linear", "nearest"), restricted = FALSE, ...)

# S3 method for class 'fm_mesh_2d'
fm_bary(mesh, loc, crs = NULL, ..., max_batch_size = NULL)

# S3 method for class 'fm_mesh_3d'
fm_bary(mesh, loc, ..., max_batch_size = NULL)

# S3 method for class 'fm_lattice_2d'
fm_bary(mesh, loc, crs = NULL, ...)

# S3 method for class 'fm_lattice_Nd'
fm_bary(mesh, loc, ...)
```

## Arguments

- ...:

  Arguments forwarded to sub-methods.

- bary:

  An `fm_bary` object, or an object that can be converted to `fm_bary`.

- extra_class:

  character; If non-`NULL` and not already in the class vector of
  `bary`, add it to the front of the class vector.

- mesh:

  `fm_mesh_1d` or `fm_mesh_2d` object

- loc:

  Points for which to identify the containing interval/triangle, and
  corresponding barycentric coordinates. May be a vector (for 1d) or a
  matrix of raw coordinates, `sf`, or `sp` point information (for 2d).

- method:

  character; method for defining the barycentric coordinates, "linear"
  (default) or "nearest"

- restricted:

  logical, used for `method="linear"`. If `FALSE` (default), points
  outside the mesh interval will be given barycentric weights less than
  0 and greater than 1, according to linear extrapolation. If `TRUE`,
  the barycentric weights are clamped to the (0, 1) interval.

- crs:

  Optional crs information for `loc`

- max_batch_size:

  integer; maximum number of points to process in a single batch. This
  speeds up calculations by avoiding repeated large internal memory
  allocations and data copies. The default, `NULL`, uses
  `max_batch_size = 2e5L`, chosen based on empirical time measurements
  to give an approximately optimal runtime.

## Value

A `fm_bary` object, a `tibble` with columns `index`; either

- vector of triangle indices (triangle meshes),

- vector of knot indices (1D meshes, either for edges or individual
  knots), or

- vector of lower left box indices (2D lattices),

and `where`, a matrix of barycentric coordinates.

## Methods (by class)

- `fm_bary(fm_bary)`: Returns the `bary` input unchanged

- `fm_bary(list)`: Converts a `list` `bary` to `fm_bary`. In the list
  elements are unnamed, the names `index` and `where` are assumed.

- `fm_bary(tbl_df)`: Converts a
  [`tibble::tibble()`](https://tibble.tidyverse.org/reference/tibble.html)
  `bary` to `fm_bary`

- `fm_bary(fm_mesh_1d)`: Return an `fm_bary` object with elements
  `index` (edge index vector pointing to the first knot of each edge)
  and `where` (barycentric coordinates, 2-column matrices). Use
  [`fm_bary_simplex()`](https://inlabru-org.github.io/fmesher/reference/fm_bary_simplex.md)
  to obtain the corresponding endpoint knot indices.

  For `method = "nearest"`, `index` contains the index of the nearest
  mesh knot, and `where` is a single-column all-ones matrix.

- `fm_bary(fm_mesh_2d)`: An `fm_bary` object with columns `index`
  (vector of triangle indices) and `where` (3-column matrix of
  barycentric coordinates). Points that were not found give `NA` entries
  in `index` and `where`.

- `fm_bary(fm_mesh_3d)`: An `fm_bary` object with columns `index`
  (vector of triangle indices) and `where` (4-column matrix of
  barycentric coordinates). Points that were not found give `NA` entries
  in `index` and `where`.

- `fm_bary(fm_lattice_2d)`: An `fm_bary` object with columns `index`
  (vector of lattice cell indices) and `where` (4-column matrix of
  barycentric coordinates). Points that are outside the lattice are
  given `NA` entries in `index` and `where`.

- `fm_bary(fm_lattice_Nd)`: An `fm_bary` object with columns `index`
  (vector of lattice cell indices) and `where` `2^d`-column matrix of
  barycentric coordinates). Points that are outside the lattice are
  given `NA` entries in `index` and `where`.

## See also

[`fm_bary_simplex()`](https://inlabru-org.github.io/fmesher/reference/fm_bary_simplex.md),
[`fm_bary_loc()`](https://inlabru-org.github.io/fmesher/reference/fm_bary_loc.md)

## Examples

``` r
bary <- fm_bary(fm_mesh_1d(1:4), seq(0, 5, by = 0.5))
bary
#> # A tibble: 11 × 2
#>    index where[,1]  [,2]
#>  * <int>     <dbl> <dbl>
#>  1     1       2    -1  
#>  2     1       1.5  -0.5
#>  3     1       1     0  
#>  4     1       0.5   0.5
#>  5     2       1     0  
#>  6     2       0.5   0.5
#>  7     3       1     0  
#>  8     3       0.5   0.5
#>  9     3       0     1  
#> 10     3      -0.5   1.5
#> 11     3      -1     2  
str(fm_bary(fmexample$mesh, fmexample$loc_sf))
#> fm_bary [10 × 2] (S3: fm_bary/tbl_df/tbl/data.frame)
#>  $ index: int [1:10] 358 301 413 337 369 221 363 329 329 142
#>  $ where: num [1:10, 1:3] 0.0699 0.1367 0.4852 0.0831 0.6245 ...
m <- fm_mesh_3d(
  rbind(
    c(1, 0, 0),
    c(0, 1, 0),
    c(0, 0, 1),
    c(0, 0, 0)
  ),
  matrix(c(1, 2, 3, 4), 1, 4)
)
b <- fm_bary(m, matrix(c(1, 1, 1) / 4, 1, 3))
str(fm_bary(fmexample$mesh, fmexample$loc_sf))
#> fm_bary [10 × 2] (S3: fm_bary/tbl_df/tbl/data.frame)
#>  $ index: int [1:10] 358 301 413 337 369 221 363 329 329 142
#>  $ where: num [1:10, 1:3] 0.0699 0.1367 0.4852 0.0831 0.6245 ...
```
