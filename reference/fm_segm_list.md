# Methods for fm_segm lists

`fm_segm` lists can be combined into `fm_segm_list` list objects.

## Usage

``` r
# S3 method for class 'fm_segm'
c(...)

# S3 method for class 'fm_segm_list'
c(...)

# S3 method for class 'fm_segm_list'
x[i]
```

## Arguments

- ...:

  Objects to be combined.

- x:

  `fm_segm_list` object from which to extract element(s)

- i:

  indices specifying elements to extract

## Value

A `fm_segm_list` object

## Methods (by generic)

- `c(fm_segm_list)`: The `...` arguments should be coercible to
  `fm_segm_list` objects.

- `[`: Extract sub-list

## Functions

- `c(fm_segm)`: The `...` arguments should be `fm_segm` objects, or
  coercible with `fm_as_segm_list(list(...))`.

## See also

[`fm_as_segm_list()`](https://inlabru-org.github.io/fmesher/reference/fm_as_segm.md)

## Examples

``` r
m <- c(A = fm_segm(1:2), B = fm_segm(3:4))
str(m)
#> List of 2
#>  $ A:List of 5
#>   ..$ loc   : num [1, 1:3] 1 2 0
#>   ..$ idx   : int [1, 1:2] 1 1
#>   ..$ grp   : NULL
#>   ..$ is.bnd: logi TRUE
#>   ..$ crs   : NULL
#>   ..- attr(*, "class")= chr [1:2] "fm_segm" "inla.mesh.segment"
#>  $ B:List of 5
#>   ..$ loc   : num [1, 1:3] 3 4 0
#>   ..$ idx   : int [1, 1:2] 1 1
#>   ..$ grp   : NULL
#>   ..$ is.bnd: logi TRUE
#>   ..$ crs   : NULL
#>   ..- attr(*, "class")= chr [1:2] "fm_segm" "inla.mesh.segment"
#>  - attr(*, "class")= chr [1:3] "fm_segm_list" "fm_list" "list"
str(m[2])
#> List of 1
#>  $ B:List of 5
#>   ..$ loc   : num [1, 1:3] 3 4 0
#>   ..$ idx   : int [1, 1:2] 1 1
#>   ..$ grp   : NULL
#>   ..$ is.bnd: logi TRUE
#>   ..$ crs   : NULL
#>   ..- attr(*, "class")= chr [1:2] "fm_segm" "inla.mesh.segment"
#>  - attr(*, "class")= chr [1:3] "fm_segm_list" "fm_list" "list"
```
