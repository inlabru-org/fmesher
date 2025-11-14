# Row-wise Kronecker products

Takes two Matrices and computes the row-wise Kronecker product.
Optionally applies row-wise weights and/or applies an additional 0/1
row-wise Kronecker matrix product.

## Usage

``` r
fm_row_kron(M1, M2, repl = NULL, n.repl = NULL, weights = NULL)
```

## Arguments

- M1:

  A matrix that can be transformed into a sparse Matrix.

- M2:

  A matrix that can be transformed into a sparse Matrix.

- repl:

  An optional index vector. For each entry, specifies which replicate
  the row belongs to, in the sense used in `INLA::inla.spde.make.A`

- n.repl:

  The maximum replicate index, in the sense used in
  `INLA::inla.spde.make.A()`.

- weights:

  Optional scaling weights to be applied row-wise to the resulting
  matrix.

## Value

A
[`Matrix::sparseMatrix`](https://rdrr.io/pkg/Matrix/man/sparseMatrix.html)
object.

## Author

Finn Lindgren <Finn.Lindgren@gmail.com>

## Examples

``` r
fm_row_kron(rbind(c(1, 1, 0), c(0, 1, 1)), rbind(c(1, 2), c(3, 4)))
#> 2 x 6 sparse Matrix of class "dgCMatrix"
#>                 
#> [1,] 1 2 1 2 . .
#> [2,] . . 3 4 3 4
```
