# Compute mapping matrix between mesh function space and points

Computes the basis mapping matrix between a function space on a mesh,
and locations.

## Usage

``` r
fm_basis(x, ..., full = FALSE)

# Default S3 method
fm_basis(x, ..., full = FALSE)

# S3 method for class 'fm_mesh_1d'
fm_basis(x, loc, weights = NULL, derivatives = NULL, ..., full = FALSE)

# S3 method for class 'fm_mesh_2d'
fm_basis(x, loc, weights = NULL, derivatives = NULL, ..., full = FALSE)

# S3 method for class 'fm_mesh_3d'
fm_basis(x, loc, weights = NULL, ..., full = FALSE)

# S3 method for class 'fm_lattice_2d'
fm_basis(x, loc, weights = NULL, ..., full = FALSE)

# S3 method for class 'fm_lattice_Nd'
fm_basis(x, loc, weights = NULL, ..., full = FALSE)

# S3 method for class 'fm_tensor'
fm_basis(x, loc, weights = NULL, ..., full = FALSE)

# S3 method for class 'fm_collect'
fm_basis(x, loc, weights = NULL, ..., full = FALSE)

# S3 method for class 'matrix'
fm_basis(x, ok = NULL, weights = NULL, ..., full = FALSE)

# S3 method for class 'Matrix'
fm_basis(x, ok = NULL, weights = NULL, ..., full = FALSE)

# S3 method for class 'list'
fm_basis(x, weights = NULL, ..., full = FALSE)

# S3 method for class 'fm_basis'
fm_basis(x, ..., full = FALSE)

# S3 method for class 'fm_evaluator'
fm_basis(x, ..., full = FALSE)
```

## Arguments

- x:

  An function space object, or other supported object (`matrix`,
  `Matrix`, `list`)

- ...:

  Passed on to submethods

- full:

  logical; if `TRUE`, return a `fm_basis` object, containing at least a
  projection matrix `A` and logical vector `ok` indicating which
  evaluations are valid. If `FALSE`, return only the projection matrix
  `A`. Default is `FALSE`.

- loc:

  A location/value information object (`numeric`, `matrix`, `sf`,
  `fm_bary`, etc, depending on the class of `x`)

- weights:

  Optional weight vector to apply (from the left, one weight for each
  row of the basis matrix)

- derivatives:

  If non-NULL and logical, include derivative matrices in the output.
  Forces `full = TRUE`.

- ok:

  numerical of length `NROW(x)`, indicating which rows of `x` are
  valid/successful basis evaluations. If `NULL`, inferred as
  `rep(TRUE, NROW(x))`.

## Value

A `sparseMatrix` object (if `full = FALSE`), or a `fm_basis` object (if
`full = TRUE` or `isTRUE(derivatives)`). The `fm_basis` object contains
at least the projection matrix `A` and logical vector `ok`; If `x_j`
denotes the latent basis coefficient for basis function `j`, the field
is defined as `u(loc_i)=sum_j A_ij x_j` for all `i` where `ok[i]` is
`TRUE`, and `u(loc_i)=0.0` where `ok[i]` is `FALSE`.

## Methods (by class)

- `fm_basis(fm_mesh_1d)`: If `derivatives=TRUE`, the `fm_basis` object
  contains additional derivative weight matrices, `d1A` and `d2A`,
  `du/dx(loc_i)=sum_j dx_ij w_i`.

- `fm_basis(fm_mesh_2d)`: If `derivatives=TRUE`, additional derivative
  weight matrices are included in the `full=TRUE` output: Derivative
  weight matrices `dx`, `dy`, `dz`; `du/dx(loc_i)=sum_j dx_ij w_i`, etc.

- `fm_basis(fm_mesh_3d)`: `fm_mesh_3d` basis functions.

- `fm_basis(fm_lattice_2d)`: `fm_lattice_2d` bilinear basis functions.

- `fm_basis(fm_lattice_Nd)`: `fm_lattice_Nd` multilinear basis
  functions.

- `fm_basis(fm_tensor)`: Evaluates a basis matrix for a `fm_tensor`
  function space.

- `fm_basis(fm_collect)`: Evaluates a basis matrix for a `fm_collect`
  function space. The `loc` argument must be a `list` or `tibble` with
  elements `loc` (the locations) and `index` (the indices into the
  function space collection).

- `fm_basis(matrix)`: Creates a new `fm_basis` object with elements `A`
  and `ok`, from a pre-evaluated basis matrix, including optional
  additional elements in the `...` arguments. If a `ok` is `NULL`, it is
  inferred as `rep(TRUE, NROW(x))`, indicating that all rows correspond
  to successful basis evaluations. If `full = FALSE`, returns the matrix
  unchanged.

- `fm_basis(Matrix)`: Creates a new `fm_basis` object with elements `A`
  and `ok`, from a pre-evaluated basis matrix, including optional
  additional elements in the `...` arguments. If a `ok` is `NULL`, it is
  inferred as `rep(TRUE, NROW(x))`, indicating that all rows correspond
  to successful basis evaluations. If `full = FALSE`, returns the matrix
  unchanged.

- `fm_basis(list)`: Creates a new `fm_basis` object from a plain list
  containing at least an element `A`. If an `ok` element is missing, it
  is inferred as `rep(TRUE, NROW(x$A))`. If `full = FALSE`, extracts the
  `A` matrix.

- `fm_basis(fm_basis)`: If `full` is `TRUE`, returns `x` unchanged,
  otherwise returns the `A` matrix contained in `x`.

- `fm_basis(fm_evaluator)`: Extract `fm_basis` information from an
  `fm_evaluator` object. If `full = FALSE`, returns the `A` matrix
  contained in the `fm_basis` object.

## See also

[`fm_raw_basis()`](https://inlabru-org.github.io/fmesher/reference/fm_raw_basis.md)

## Examples

``` r
# Compute basis mapping matrix
dim(fm_basis(fmexample$mesh, fmexample$loc))
#> [1]  10 279
print(fm_basis(fmexample$mesh, fmexample$loc, full = TRUE))
#> fm_basis object
#>   Projection matrix (A): 10-by-279
#>   Valid evaluations (ok): 10 out of 10
#>   Additional information: bary

# From precomputed `fm_bary` information:
bary <- fm_bary(fmexample$mesh, fmexample$loc)
print(fm_basis(fmexample$mesh, bary, full = TRUE))
#> fm_basis object
#>   Projection matrix (A): 10-by-279
#>   Valid evaluations (ok): 10 out of 10
#>   Additional information: bary
```
