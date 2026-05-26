# Custom mesh classes

``` r

library(fmesher)
```

## Minimal interface

Users and package developers can add `fmesher` support to their own
classes. A minimal interface needs to define
[`fm_dof()`](https://inlabru-org.github.io/fmesher/reference/fm_dof.md)
and
[`fm_basis()`](https://inlabru-org.github.io/fmesher/reference/fm_basis.md)
methods. Assuming the class is called `custom`, the methods should be
named `fm_dof.custom()` and `fm_basis.custom()`.

- The `fm_dof.custom(x)` method *must* take an object `x` of class
  `custom` and return the number of degrees of freedom of the function
  space.

- The `fm_basis.custom(x, loc, ..., full = FALSE)` method *must* take an
  object `x` of class `custom` and return a `sparseMatrix` or `Matrix`
  matrix with each column containing the basis function evaluated at the
  locations determined by `loc`.

  The type of `loc` may be any type (or types) that is supported by the
  custom class.

  The `...` part can include further named arguments specific to the
  custom class. These must be optional arguments so that
  `fm_basis(x, loc)` works.

  When `full = TRUE`, a full `fm_basis` object *must* be returned, which
  is a list containing *at least* the basis matrix as `A`, and a logical
  vector, `ok`, indicating which `loc` values were valid evaluation
  points. The `A` matrix *must* be all-zero for invalid `loc`.

- With the above requirements fulfilled, the default
  [`fm_evaluator()`](https://inlabru-org.github.io/fmesher/reference/fm_evaluate.md)
  and
  [`fm_evaluate()`](https://inlabru-org.github.io/fmesher/reference/fm_evaluate.md)
  methods can be used to evaluate functions at any location, without the
  need for the user to define any further methods.

  Special `fm_evaluator.custom()` and `fm_evaluate.custom()` methods may
  be defined if needed, e.g. to support semi-automated output
  reformatting.

### Example: Harmonic function space of order n

``` r

# Custom class for harmonic functions up to order `n`
create_custom <- function(n) {
  stopifnot(n >= 0)
  structure(
    list(n = n),
    class = "custom"
  )
}
fm_dof.custom <- function(x) {
  # Return the number of degrees of freedom
  1L + 2L * x[["n"]]
}
fm_basis.custom <- function(x, loc, ..., full = FALSE) {
  # Return the evaluated basis functions
  A <- Matrix::Matrix(0.0, NROW(loc), fm_dof(x))
  ok <- !is.na(loc)
  A[ok, 1L] <- 1.0
  for (k in seq_len(x[["n"]])) {
    A[ok, 2 * k] <- cos(2 * pi * k * loc[ok])
    A[ok, 2 * k + 1L] <- sin(2 * pi * k * loc[ok])
  }
  result <- structure(
    list(
      A = A,
      ok = ok, # Required prior to version 0.2.0.9003
      loc = loc
    ),
    class = "fm_basis"
  )

  # Use the fm_basis method to extract the A matrix if full is FALSE:
  fm_basis(result, full = full)
}
```

Note: From version `0.2.0.9004`, the `fm_basis.matrix`,
`fm_basis.Matrix`, and `fm_basis.list` methods provide an easier way to
construct the `fm_basis` object, by creating the object and optionally
extracting `A` in a single call:

``` r

# 'matrix' and 'Matrix' methods:
fm_basis(
  A = A,
  ok = ok, # If missing or NULL, inferred to be all TRUE
  loc = loc, # Optional additional content
  full = full
)
# 'list' method:
fm_basis(
  list(
    A = A,
    ok = ok, # If missing or NULL, inferred to be all TRUE
    loc = loc
  ),
  full = full
)
```

#### Registering the methods

These S3 methods must be registered with the `S3method()` function in
scripts, and with special NAMESPACE tags in packages. In a script, one
should use

``` r

.S3method("fm_dof", "custom", "fm_dof.custom")
.S3method("fm_basis", "custom", "fm_basis.custom")
```

In a package, if R is version 3.6 or newer, one can use roxygen2 tags

``` r

#' @rawNamespace S3method(fmesher::fm_dof, custom)
#' @rawNamespace S3method(fmesher::fm_basis, custom)
```

or before each method, use `@exportS3Method`, like this:

``` r

#' @title Degrees of freedom for custom mesh
#' @description the number of degrees of freedom
#' # The rest of the documentation goes here
#' @exportS3method fmesher::fm_dof
fm_dof.custom <- function(x) {
  1L + 2L * x[["n"]]
}
```

which semi-automates it.

We can the use the new methods with

``` r

m <- create_custom(2)

# How many latent variables are needed?
fm_dof(m)
#> [1] 5

# Evaluate the basis functions at some locations:
fm_basis(m, seq(0, 1, length.out = 6))
#> 6 x 5 sparse Matrix of class "dgCMatrix"
#>                                                       
#> [1,] 1  1.000000  .             1.000000  .           
#> [2,] 1  0.309017  9.510565e-01 -0.809017  5.877853e-01
#> [3,] 1 -0.809017  5.877853e-01  0.309017 -9.510565e-01
#> [4,] 1 -0.809017 -5.877853e-01  0.309017  9.510565e-01
#> [5,] 1  0.309017 -9.510565e-01 -0.809017 -5.877853e-01
#> [6,] 1  1.000000 -2.449294e-16  1.000000 -4.898587e-16
fm_basis(m, seq(0, 1, length.out = 6), full = TRUE)
#> fm_basis object
#>   Projection matrix (A): 6-by-5
#>   Valid evaluations (ok): 6 out of 6
#>   Additional information: loc

# Check if missing values are handled correctly:
fm_basis(m, c(0.1, NA, 0.2))
#> 3 x 5 sparse Matrix of class "dgCMatrix"
#>                                              
#> [1,] 1 0.809017 0.5877853  0.309017 0.9510565
#> [2,] . .        .          .        .        
#> [3,] 1 0.309017 0.9510565 -0.809017 0.5877853
fm_basis(m, c(0.1, NA, 0.2), full = TRUE)
#> fm_basis object
#>   Projection matrix (A): 3-by-5
#>   Valid evaluations (ok): 2 out of 3
#>   Additional information: loc
```

## Expanded implementations

The main additional method that can be defined is the
[`fm_int()`](https://inlabru-org.github.io/fmesher/reference/fm_int.md)
integration scheme method. This must have the call structure
`fm_int.custom(domain, samplers = NULL, name = "x", ...)`.

- The `domain` argument is the `custom` class object over which to
  integrate.
- The `samplers` argument is any object, typically an `sf` or `tibble`
  specifying one or more subsets of the domain, e.g. polygons. When
  `NULL`, the entire domain should be integrated.
- The `name` argument is a character string specifying the name of the
  integration point variable.
- The `...` arguments can be augmented with further optional arguments,
  e.g. options controlling the integration scheme construction and/or
  the output format.

`out <- fm_int(domain, samplers, name)` should return a `data.frame`,
`tibble`, or `sf` object with integration points in a column with the
name indicated by , and additional columns `weight` with corresponding
integration weights, and a `.block` column.

- The column format should be compatible with
  `fm_basis(domain, out[[name]])`.
- The `.block` column should be an integer vector indicating which
  subdomain each integration point belongs to, usable by
  [`fm_block_eval()`](https://inlabru-org.github.io/fmesher/reference/fm_block.md):

&nbsp;

    #>      x weight .block
    #> 1  0.0   0.05      1
    #> 2  0.1   0.10      1
    #> 3  0.2   0.10      1
    #> 4  0.3   0.05      1
    #> 5  0.3   0.05      2
    #> 6  0.4   0.10      2
    #> 7  0.5   0.05      2
    #> 8  0.5   0.05      3
    #> 9  0.6   0.10      3
    #> 10 0.7   0.10      3
    #> 11 0.8   0.10      3
    #> 12 0.9   0.10      3
    #> 13 1.0   0.05      3

``` r

values <- fm_evaluate(
  m,
  field = c(1, 1, 0, 0, 0),
  loc = out[["x"]]
)

# Blockwise aggregation:
fm_block_eval(
  block = out$.block,
  weights = out$weight,
  values = values
)
#> [1] 0.44635255 0.05364745 0.50000000

# Exact integrals:
c(0.3, 0.2, 0.5) +
  c(
    sin(2 * pi * 0.3),
    sin(2 * pi * 0.5) - sin(2 * pi * 0.3),
    sin(2 * pi) - sin(2 * pi * 0.5)
  ) / (2 * pi)
#> [1] 0.45136535 0.04863465 0.50000000
```
