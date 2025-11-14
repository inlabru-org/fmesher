# Bounding box class

Simple class for handling bounding box information

## Usage

``` r
fm_bbox(...)

# S3 method for class 'list'
fm_bbox(x, ...)

## S3 method for class 'NULL'
fm_bbox(...)

# S3 method for class 'numeric'
fm_bbox(x, ...)

# S3 method for class 'matrix'
fm_bbox(x, ...)

# S3 method for class 'Matrix'
fm_bbox(x, ...)

# S3 method for class 'fm_bbox'
fm_bbox(x, ...)

# S3 method for class 'fm_mesh_1d'
fm_bbox(x, ...)

# S3 method for class 'fm_mesh_2d'
fm_bbox(x, ...)

# S3 method for class 'fm_mesh_3d'
fm_bbox(x, ...)

# S3 method for class 'fm_segm'
fm_bbox(x, ...)

# S3 method for class 'fm_lattice_2d'
fm_bbox(x, ...)

# S3 method for class 'fm_lattice_Nd'
fm_bbox(x, ...)

# S3 method for class 'fm_tensor'
fm_bbox(x, ...)

# S3 method for class 'fm_collect'
fm_bbox(x, ...)

# S3 method for class 'sf'
fm_bbox(x, ...)

# S3 method for class 'sfg'
fm_bbox(x, ...)

# S3 method for class 'sfc'
fm_bbox(x, ...)

# S3 method for class 'bbox'
fm_bbox(x, ...)

fm_as_bbox(x, ...)

# S3 method for class 'fm_bbox'
x[i]

# S3 method for class 'fm_bbox'
c(..., .join = FALSE)

fm_as_bbox_list(x, ...)
```

## Arguments

- ...:

  Passed on to sub-methods

- x:

  `fm_bbox` object from which to extract element(s)

- i:

  indices specifying elements to extract

- .join:

  logical; if `TRUE`, concatenate the bounding boxes into a single
  multi-dimensional bounding box. Default is `FALSE`.

## Value

For `c.fm_bbox()`, a `fm_bbox_list` object if `.join = FALSE` (the
default) or an `fm_bbox` object if `.join = TRUE`.

## Methods (by class)

- `fm_bbox(list)`: Construct a bounding box from precomputed interval
  information, stored as a list of 2-vector ranges,
  `list(xlim, ylim, ...)`.

## Methods (by generic)

- `[`: Extract sub-list

- `c(fm_bbox)`: The `...` arguments should be `fm_bbox` objects, or
  coercible with `fm_as_bbox(list(...))`.

## Functions

- `fm_as_bbox_list()`: Convert a list to a `fm_bbox_list` object, with
  each element converted to an `fm_bbox` object.

## Examples

``` r
fm_bbox(matrix(1:6, 3, 2))
#> Bounding box: (1,3) x (4,6)
m <- c(A = fm_bbox(cbind(1, 2)), B = fm_bbox(cbind(3, 4)))
str(m)
#> List of 2
#>  $ A:List of 2
#>   ..$ : num [1:2] 1 1
#>   ..$ : num [1:2] 2 2
#>   ..- attr(*, "class")= chr "fm_bbox"
#>  $ B:List of 2
#>   ..$ : num [1:2] 3 3
#>   ..$ : num [1:2] 4 4
#>   ..- attr(*, "class")= chr "fm_bbox"
#>  - attr(*, "class")= chr [1:3] "fm_bbox_list" "fm_list" "list"
str(m[2])
#> List of 1
#>  $ B:List of 2
#>   ..$ : num [1:2] 3 3
#>   ..$ : num [1:2] 4 4
#>   ..- attr(*, "class")= chr "fm_bbox"
#>  - attr(*, "class")= chr [1:3] "fm_bbox_list" "fm_list" "list"
m <- fm_as_bbox_list(list(
  A = fm_bbox(cbind(1, 2)),
  B = fm_bbox(cbind(3, 4))
))
str(fm_as_bbox_list(m))
#> List of 2
#>  $ A:List of 2
#>   ..$ : num [1:2] 1 1
#>   ..$ : num [1:2] 2 2
#>   ..- attr(*, "class")= chr "fm_bbox"
#>  $ B:List of 2
#>   ..$ : num [1:2] 3 3
#>   ..$ : num [1:2] 4 4
#>   ..- attr(*, "class")= chr "fm_bbox"
#>  - attr(*, "class")= chr [1:3] "fm_bbox_list" "fm_list" "list"
```
