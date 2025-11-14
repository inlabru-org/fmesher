# Conversion between sparse matrix types

Conversion between sparse matrix types

## Usage

``` r
fm_as_dgCMatrix(x)

fm_as_dgTMatrix(x, unique = TRUE, ...)

fm_as_unpackedMatrix(x)

fm_as_fmesher_sparse(x)

# Default S3 method
fm_as_dgCMatrix(x)

# S3 method for class 'fmesher_sparse'
fm_as_dgCMatrix(x)

# Default S3 method
fm_as_dgTMatrix(x, unique = TRUE, ...)

# Default S3 method
fm_as_unpackedMatrix(x)

# S3 method for class 'fmesher_sparse'
fm_as_unpackedMatrix(x)

# S3 method for class 'fmesher_sparse'
fm_as_dgTMatrix(x, unique = TRUE, ...)
```

## Arguments

- x:

  Object to be converted

- unique:

  logical; if `TRUE`, ensures that the sparse triplet representation has
  a single entry for each non-zero matrix element.

## Value

`fm_as_dgCMatrix` returns a
[Matrix::dgCMatrix](https://rdrr.io/pkg/Matrix/man/dgCMatrix-class.html)
object.

`fm_as_dgTMatrix` returns a
[Matrix::dgTMatrix](https://rdrr.io/pkg/Matrix/man/dgTMatrix-class.html)
object.

`fm_as_unpackedMatrix` returns an object of virtual class
[Matrix::unpackedMatrix](https://rdrr.io/pkg/Matrix/man/unpackedMatrix-class.html).

`fm_as_fmesher_sparse` returns an `fmesher_sparse` object.

## Examples

``` r
library(Matrix)
str(A <- fm_as_dgCMatrix(matrix(c(1, 2, 0, 0, 0, 3, 4, 0, 5), 3, 3)))
#> Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
#>   ..@ i       : int [1:5] 0 1 2 0 2
#>   ..@ p       : int [1:4] 0 2 3 5
#>   ..@ Dim     : int [1:2] 3 3
#>   ..@ Dimnames:List of 2
#>   .. ..$ : NULL
#>   .. ..$ : NULL
#>   ..@ x       : num [1:5] 1 2 3 4 5
#>   ..@ factors : list()
str(fm_as_dgTMatrix(A))
#> Formal class 'dgTMatrix' [package "Matrix"] with 6 slots
#>   ..@ i       : int [1:5] 0 1 2 0 2
#>   ..@ j       : int [1:5] 0 0 1 2 2
#>   ..@ Dim     : int [1:2] 3 3
#>   ..@ Dimnames:List of 2
#>   .. ..$ : NULL
#>   .. ..$ : NULL
#>   ..@ x       : num [1:5] 1 2 3 4 5
#>   ..@ factors : list()
str(fm_as_unpackedMatrix(A))
#> Formal class 'dgeMatrix' [package "Matrix"] with 4 slots
#>   ..@ Dim     : int [1:2] 3 3
#>   ..@ Dimnames:List of 2
#>   .. ..$ : NULL
#>   .. ..$ : NULL
#>   ..@ x       : num [1:9] 1 2 0 0 0 3 4 0 5
#>   ..@ factors : list()
str(fm_as_fmesher_sparse(A))
#> List of 4
#>  $ i   : int [1:5] 0 1 2 0 2
#>  $ j   : int [1:5] 0 0 1 2 2
#>  $ x   : num [1:5] 1 2 3 4 5
#>  $ dims: int [1:2] 3 3
#>  - attr(*, "class")= chr "fmesher_sparse"
```
