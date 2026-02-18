# Function spece degrees of freedom

Obtain the degrees of freedom of a function space, i.e. the number of
basis functions it uses.

## Usage

``` r
fm_dof(x)

# S3 method for class 'fm_mesh_1d'
fm_dof(x)

# S3 method for class 'fm_mesh_2d'
fm_dof(x)

# S3 method for class 'fm_mesh_3d'
fm_dof(x)

# S3 method for class 'fm_tensor'
fm_dof(x)

# S3 method for class 'fm_collect'
fm_dof(x)

# S3 method for class 'fm_lattice_2d'
fm_dof(x)

# S3 method for class 'fm_lattice_Nd'
fm_dof(x)
```

## Arguments

- x:

  A function space object, such as
  [`fm_mesh_1d()`](https://inlabru-org.github.io/fmesher/reference/fm_mesh_1d.md)
  or
  [`fm_mesh_2d()`](https://inlabru-org.github.io/fmesher/reference/fm_mesh_2d.md)

## Value

An integer

## Examples

``` r
fm_dof(fmexample$mesh)
#> [1] 292
```
