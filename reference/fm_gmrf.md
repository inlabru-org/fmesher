# SPDE, GMRF, and Matérn process methods

**\[experimental\]** Methods for SPDEs and GMRFs.

## Usage

``` r
fm_matern_precision(x, alpha, rho, sigma)

fm_matern_sample(x, alpha = 2, rho, sigma, n = 1, loc = NULL)

fm_covariance(Q, A1 = NULL, A2 = NULL, partial = FALSE)

fm_sample(n, Q, mu = 0, constr = NULL)
```

## Arguments

- x:

  A mesh object, e.g. from
  [`fm_mesh_1d()`](https://inlabru-org.github.io/fmesher/reference/fm_mesh_1d.md),
  [`fm_mesh_2d()`](https://inlabru-org.github.io/fmesher/reference/fm_mesh_2d.md),
  or other object with supporting
  [`fm_fem()`](https://inlabru-org.github.io/fmesher/reference/fm_fem.md)
  and
  [`fm_manifold_dim()`](https://inlabru-org.github.io/fmesher/reference/fm_manifold.md)
  methods.

- alpha:

  The SPDE operator order. The resulting smoothness index is
  `nu = alpha - dim / 2`. Supports integers 1, 2, 3, etc. that give
  `nu > 0`.

- rho:

  The Matérn range parameter (scale parameter
  `kappa = sqrt(8 * nu) / rho`)

- sigma:

  The nominal Matérn std.dev. parameter

- n:

  The number of samples to generate

- loc:

  locations to evaluate the random field, compatible with
  `fm_evaluate(x, loc = loc, field = ...)`

- Q:

  A precision matrix

- A1, A2:

  Matrices, typically obtained from
  [`fm_basis()`](https://inlabru-org.github.io/fmesher/reference/fm_basis.md)
  and/or
  [`fm_block()`](https://inlabru-org.github.io/fmesher/reference/fm_block.md).

- partial:

  **\[experimental\]** If `TRUE`, compute the partial inverse of `Q`,
  i.e. the elements of the inverse corresponding to the non-zero pattern
  of `Q`. (Note: This can be done efficiently with the Takahashi
  recursion method, but to avoid an RcppEigen dependency this is
  currently disabled, and a slower method is used until the efficient
  method is reimplemented.)

- mu:

  Optional mean vector

- constr:

  Optional list of constraint information, with elements `A` and `e`.
  Should only be used for a small number of exact constraints.

## Value

`fm_matern_sample()` returns a matrix, where each column is a sampled
field. If `loc` is `NULL`, the `fm_dof(mesh)` basis weights are given.
Otherwise, the evaluated field at the `nrow(loc)` locations `loc` are
given (from version `0.1.4.9001`)

## Functions

- `fm_matern_precision()`: Construct the (sparse) precision matrix for
  the basis weights for Whittle-Matérn SPDE models. The boundary
  behaviour is determined by the provided mesh function space.

- `fm_matern_sample()`: Simulate a Matérn field given a mesh and
  covariance function parameters, and optionally evaluate at given
  locations.

- `fm_covariance()`: Compute the covariance between "A1 x" and "A2 x",
  when x is a basis vector with precision matrix `Q`.

- `fm_sample()`: Generate `n` samples based on a sparse precision matrix
  `Q`

## Examples

``` r
library(Matrix)
mesh <- fm_mesh_1d(-20:120, degree = 2)
Q <- fm_matern_precision(mesh, alpha = 2, rho = 15, sigma = 1)
x <- seq(0, 100, length.out = 601)
A <- fm_basis(mesh, x)
plot(x,
  as.vector(Matrix::diag(fm_covariance(Q, A))),
  type = "l",
  ylab = "marginal variances"
)


plot(x,
  fm_evaluate(mesh, loc = x, field = fm_sample(1, Q)[, 1]),
  type = "l",
  ylab = "process sample"
)

```
