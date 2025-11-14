# Sparse partial inverse

Compute sparse partial matrix inverse. As of `0.2.0.9010`, an R
implementation of the Takahashi recursion method, unless a special build
of the `fmesher` package is used.

## Usage

``` r
fm_qinv(A)
```

## Arguments

- A:

  A sparse symmetric positive definite matrix

## Value

A sparse symmetric matrix, with the elements of the inverse of `A` for
the non-zero pattern of `A` plus potential Cholesky in-fill locations.

## Examples

``` r
A <- Matrix::Matrix(
  c(2, -1, 0, 0, -1, 2, -1, 0, 0, -1, 2, -1, 0, 0, -1, 2),
  4,
  4
)
# Partial inverse:
(S <- fm_qinv(A))
#> 4 x 4 sparse Matrix of class "dgCMatrix"
#>                     
#> [1,] 0.8 0.6 0.4 .  
#> [2,] 0.6 1.2 0.8 .  
#> [3,] 0.4 0.8 1.2 0.6
#> [4,] .   .   0.6 0.8
# Full inverse (not guaranteed to be symmetric):
(S2 <- solve(A))
#> 4 x 4 Matrix of class "dsyMatrix"
#>      [,1] [,2] [,3] [,4]
#> [1,]  0.8  0.6  0.4  0.2
#> [2,]  0.6  1.2  0.8  0.4
#> [3,]  0.4  0.8  1.2  0.6
#> [4,]  0.2  0.4  0.6  0.8
# Matrix symmetry:
c(sum((S - Matrix::t(S))^2), sum((S2 - Matrix::t(S2))^2))
#> [1] 0 0
# Accuracy (not that S2 is non-symmetric, and S may be more accurate):
sum((S - S2)[S != 0]^2)
#> [1] 7.395571e-32
```
