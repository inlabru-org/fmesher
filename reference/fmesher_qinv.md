# Compute sparse matrix inverse

Internal C++ method.

Requires RcppEigen which is not compiled in by default. Enable with
`PKG_CPPFLAGS=-DFMESHER_WITH_EIGEN` in `src/Makevars` and add
`RcppEigen` to the `DESCRIPTION` `LinkingTo` field.

## Usage

``` r
fmesher_qinv(AA)
```

## Arguments

- AA:

  A sparse matrix
