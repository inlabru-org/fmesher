# Query the mesh manifold type

Extract a manifold definition string, or a logical for matching manifold
type

## Usage

``` r
fm_manifold(x, type = NULL)

fm_manifold_get(x)

# Default S3 method
fm_manifold_get(x)

# S3 method for class 'character'
fm_manifold_get(x)

# S3 method for class 'fm_lattice_2d'
fm_manifold_get(x)

# S3 method for class 'fm_lattice_Nd'
fm_manifold_get(x)

fm_manifold_type(x)

fm_manifold_dim(x)
```

## Arguments

- x:

  An object with `manifold` information, or a character string

- type:

  `character`; if `NULL` (the default), returns the manifold definition
  string by calling `fm_manifold_get(x)`. If `character`, returns `TRUE`
  if the manifold type of `x` matches at least one of the character
  vector elements.

## Value

`fm_manifold()`: Either logical (matching manifold type yes/no), or
character (the stored manifold, when `is.null(type)` is `TRUE`)

`fm_manifold_get()`: `character` or `NULL`

`fm_manifold_type()`: character or NULL; "M" (curved manifold), "R"
(flat space), "S" (generalised spherical space), "T" (general tensor
product space), or "G" (metric graph)

`fm_manifold_dim()`: integer or NULL

## Functions

- `fm_manifold_get()`: Method for obtaining a text representation of the
  manifold characteristics, e.g. "R1", "R2", "M2", or "T3". The default
  method assumes that the manifold is stored as a `character` string in
  a "manifold" element of the object, so it can be extracted with
  `x[["manifold"]]`. Object classes that do not store the information in
  this way need to implement their own method.

## Examples

``` r
fm_manifold_get(fmexample$mesh)
#> [1] "R2"
fm_manifold(fmexample$mesh)
#> [1] "R2"
fm_manifold(fmexample$mesh, "R2")
#> [1] TRUE
fm_manifold_type(fmexample$mesh)
#> [1] "R"
fm_manifold_dim(fmexample$mesh)
#> [1] 2
```
